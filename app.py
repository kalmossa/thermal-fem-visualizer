"""API REST du visualiseur FEM thermo-mécanique.

Endpoints :

    GET  /api/mesh              → maillage courant (nœuds, connectivité, bords)
    POST /api/run-simulation    → lance un calcul, renvoie l'identifiant et les résultats
    GET  /api/results/<id>      → relit un calcul passé
    GET  /api/results/<id>.csv  → export CSV des déplacements de surface
    GET  /api/history           → liste des calculs de la session
    GET  /api/materials         → bibliothèque de matériaux disponibles
    POST /api/compare-materials → superpose les profils de plusieurs matériaux
    POST /api/twin-rail         → flexion de la traverse sous les deux rails (Winkler)
    POST /api/profile           → compare des profils d'assise multicouches
    GET  /api/health            → état du service

    POST /api/simulate          → alias historique de /api/run-simulation

Lancement en développement :

    py -m pip install -r requirements.txt
    py app.py

En production, `debug` reste à False et le service passe derrière un serveur
WSGI (waitress sous Windows, gunicorn ailleurs) — voir README.
"""
from __future__ import annotations

import json
import logging
import os
import threading
import time
import uuid
from collections import deque
from datetime import datetime, timezone

import numpy as np
from flask import Flask, Response, jsonify, request, send_from_directory

import materials as MAT
import track_beam as TB
from fem_ballast_beton import (
    BALLAST_PARAMS,
    Layer,
    SolverDivergence,
    CONCRETE_PARAMS,
    CrushableCap,
    FEMSolver,
    LinearElastic,
    Mesh,
)

# ─────────────────────────────────────────────────────────────────────────────
# Configuration
# ─────────────────────────────────────────────────────────────────────────────
ALLOWED_ORIGIN = os.environ.get("FEM_ALLOWED_ORIGIN", "*")
MAX_RUNS_PER_MINUTE = int(os.environ.get("FEM_RATE_LIMIT", "20"))
HISTORY_SIZE = 50
SOLVE_TIMEOUT_S = 30.0

app = Flask(__name__, static_folder=".")

# Un corps de requête ne dépasse jamais quelques kilo-octets ici : plafonner à
# 1 Mo coupe court à un envoi destiné à saturer la mémoire. Flask répond 413.
app.config["MAX_CONTENT_LENGTH"] = 1 * 1024 * 1024

# Journalisation : sans trace, un plantage en production ne laisse rien à
# analyser. Le fichier reste local et n'est jamais servi.
if not app.debug:
    _handler = logging.FileHandler(
        os.environ.get("FEM_LOG_FILE", "fem_server.log"), encoding="utf-8")
    _handler.setFormatter(logging.Formatter(
        "%(asctime)s %(levelname)-8s %(message)s"))
    _handler.setLevel(logging.INFO)
    app.logger.addHandler(_handler)
    app.logger.setLevel(logging.INFO)

# Un calcul mobilise un cœur pendant ~0,7 s. Sans garde-fou, une poignée de
# requêtes simultanées suffit à saturer l'instance : on sérialise les calculs et
# on limite la cadence par client.
_solve_lock = threading.Lock()
_rate_hits: dict[str, deque] = {}
_rate_lock = threading.Lock()

_history: deque = deque(maxlen=HISTORY_SIZE)
_results: dict[str, dict] = {}
_cache: dict[tuple, dict] = {}
_store_lock = threading.Lock()


# ─────────────────────────────────────────────────────────────────────────────
# Validation des entrées
# ─────────────────────────────────────────────────────────────────────────────
class ParamError(ValueError):
    """Paramètre hors domaine ou non numérique."""


MSG_DIVERGENCE = (
    "le modèle ne couvre pas cette combinaison de paramètres : le calcul a "
    "divergé avant convergence. Adoucis le chargement, rapproche le module du "
    "ballast de sa plage nominale, ou réduis l'écart de température."
)


def verifier_fini(charge, chemin="réponse"):
    """Garde-fou : refuse de sérialiser un NaN ou un infini.

    `jsonify` écrit `NaN` tel quel, ce qui n'est pas du JSON valide : le client
    reçoit un 200 et un corps que `response.json()` ne sait pas lire. Mieux vaut
    un 422 explicite qu'un succès mensonger.
    """
    pile = [(chemin, charge)]
    while pile:
        ou, v = pile.pop()
        if isinstance(v, float):
            if not np.isfinite(v):
                raise SolverDivergence(f"valeur non finie dans {ou}", etape="sérialisation")
        elif isinstance(v, dict):
            pile.extend((f"{ou}.{k}", x) for k, x in v.items())
        elif isinstance(v, (list, tuple)):
            pile.extend((f"{ou}[{i}]", x) for i, x in enumerate(v))
    return charge


# nom → (défaut, min, max, unité)
BOUNDS = {
    "p":       (300.0,   0.0,   2000.0,  "kN/m"),
    "delta_T": (40.0,   -60.0,  80.0,    "°C"),
    "E_b":     (150.0,   10.0,  2000.0,  "MPa"),
    "E_c":     (30000.0, 1000.0, 100000.0, "MPa"),
    "alpha_b": (1.2e-5,  0.0,   5.0e-5,  "1/°C"),
    "alpha_c": (1.0e-5,  0.0,   5.0e-5,  "1/°C"),
    "span":    (0.25,    0.05,  0.60,    "m"),
    "lx":      (0.80,    0.30,  2.00,    "m"),
    "ly":      (0.35,    0.10,  1.00,    "m"),
    "nx":      (32,      8,     64,      "éléments"),
    "ny":      (14,      4,     32,      "éléments"),
    "nsteps":  (6,       1,     20,      "incréments"),
}


def parse_params(data: dict) -> dict:
    """Lit, type et borne les paramètres. Lève ParamError avec un message clair."""
    if not isinstance(data, dict):
        raise ParamError("le corps de la requête doit être un objet JSON")

    out = {}
    for name, (default, lo, hi, unit) in BOUNDS.items():
        raw = data.get(name, default)
        try:
            val = float(raw)
        except (TypeError, ValueError):
            raise ParamError(f"« {name} » doit être un nombre (reçu : {raw!r})")
        if not np.isfinite(val):
            raise ParamError(f"« {name} » doit être fini")
        if not (lo <= val <= hi):
            raise ParamError(f"« {name} » doit être compris entre {lo} et {hi} {unit} (reçu : {val})")
        out[name] = int(val) if name in ("nx", "ny", "nsteps") else val

    if out["span"] >= out["lx"]:
        raise ParamError("la semelle (span) doit être plus étroite que le domaine (lx)")
    return out


# ─────────────────────────────────────────────────────────────────────────────
# Limitation de cadence
# ─────────────────────────────────────────────────────────────────────────────
def rate_limited(client: str) -> bool:
    now = time.monotonic()
    with _rate_lock:
        hits = _rate_hits.setdefault(client, deque())
        while hits and now - hits[0] > 60.0:
            hits.popleft()
        if len(hits) >= MAX_RUNS_PER_MINUTE:
            return True
        hits.append(now)
        return False


# ─────────────────────────────────────────────────────────────────────────────
# Calcul
# ─────────────────────────────────────────────────────────────────────────────
def run_simulation(p: dict) -> dict:
    """Lance les deux calculs et met en forme la réponse."""
    mesh = Mesh(p["lx"], p["ly"], p["nx"], p["ny"])
    p_service = p["p"] * 1e3            # kN/m → N/m
    t0 = time.perf_counter()

    # H_cap est calé pour E0 = 150 MPa ; le curseur va bien au-delà. Un
    # squelette granulaire plus raide se consolide aussi plus raidement, donc on
    # échelonne l'écrouissage du cap sur le module plutôt que de le figer — sans
    # quoi le retour radial décroche dès qu'on s'éloigne de la valeur nominale.
    ballast = CrushableCap(**dict(BALLAST_PARAMS, E0=p["E_b"] * 1e6,
                                  H_cap=RATIO_H_CAP * p["E_b"] * 1e6))
    s_b = FEMSolver(mesh, ballast, is_plastic=True,
                    delta_T=p["delta_T"], alpha=p["alpha_b"])
    U_b, *_ = s_b.run(p["nsteps"], p_service, p["span"])
    x, uy_b = s_b.top_profile(U_b)

    concrete = LinearElastic(p["E_c"] * 1e6, CONCRETE_PARAMS["nu"])
    s_c = FEMSolver(mesh, concrete, is_plastic=False,
                    delta_T=p["delta_T"], alpha=p["alpha_c"])
    U_c, *_ = s_c.run(p["nsteps"], p_service, p["span"])
    _, uy_c = s_c.top_profile(U_c)

    solve_ms = (time.perf_counter() - t0) * 1e3

    # La cuvette isole la réponse mécanique du soulèvement thermique d'ensemble :
    # c'est elle qui compare réellement la rigidité des deux matériaux.
    bowl_b = float(uy_b[0] - uy_b.min())
    bowl_c = float(uy_c[0] - uy_c.min())

    return {
        "id": f"run_{uuid.uuid4().hex[:8]}",
        "timestamp": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "params": p,
        "mesh": {"lx": p["lx"], "ly": p["ly"], "nx": p["nx"], "ny": p["ny"],
                 "n_nodes": int(mesh.nodes.shape[0]), "n_elements": int(mesh.conn.shape[0])},
        "x_cm": (x * 100).round(4).tolist(),
        "uy_ballast": uy_b.round(4).tolist(),
        "uy_beton": uy_c.round(4).tolist(),
        "summary": {
            "uy_min_ballast": float(uy_b.min()), "uy_max_ballast": float(uy_b.max()),
            "uy_min_beton": float(uy_c.min()), "uy_max_beton": float(uy_c.max()),
            "bowl_ballast": bowl_b,
            "bowl_beton": bowl_c,
            "ratio_bowl": bowl_b / bowl_c if bowl_c else None,
            "evp_max_pct": float(s_b.state_max("evp") * 100),
            "B_max": float(s_b.state_max("B")),
            "pc_max_kpa": float(s_b.state_max("pc") / 1e3),
            "solve_ms": round(solve_ms, 1),
            "iterations_ballast": s_b.iterations,
            "iterations_beton": s_c.iterations,
        },
        "fields": {
            # par élément : cartes de contraintes et d'état plastique
            "vm_ballast": s_b.stress_field(U_b).round(3).tolist(),
            "vm_beton": s_c.stress_field(U_c).round(3).tolist(),
            "evp_ballast": (s_b.gauss_field("evp") * 100).round(5).tolist(),
            "pc_ballast": (s_b.gauss_field("pc") / 1e3).round(3).tolist(),
            # par nœud : maillage déformé
            "uy_nodal_ballast": (U_b[1::2] * 1e6).round(3).tolist(),
            "uy_nodal_beton": (U_c[1::2] * 1e6).round(3).tolist(),
            "ux_nodal_ballast": (U_b[0::2] * 1e6).round(3).tolist(),
            "ux_nodal_beton": (U_c[0::2] * 1e6).round(3).tolist(),
        },
    }


# ─────────────────────────────────────────────────────────────────────────────
# Routes
# ─────────────────────────────────────────────────────────────────────────────
@app.after_request
def add_headers(response):
    response.headers["Access-Control-Allow-Origin"] = ALLOWED_ORIGIN
    response.headers["Access-Control-Allow-Methods"] = "POST, GET, OPTIONS"
    response.headers["Access-Control-Allow-Headers"] = "Content-Type"
    # Empêche le navigateur de deviner le type MIME, qui ouvre la porte à des
    # injections par contenu mal typé.
    response.headers["X-Content-Type-Options"] = "nosniff"
    # Interdit l'affichage dans une iframe : protection contre le clickjacking.
    response.headers["X-Frame-Options"] = "DENY"
    response.headers["Referrer-Policy"] = "strict-origin-when-cross-origin"
    if request.is_secure:
        response.headers["Strict-Transport-Security"] = "max-age=31536000; includeSubDomains"
    # Un résultat de simulation est calculé à la demande : il ne doit pas être
    # resservi depuis le cache du navigateur.
    if request.path.startswith("/api/"):
        response.headers["Cache-Control"] = "no-store"
    return response


@app.errorhandler(413)
def trop_volumineux(_):
    return jsonify({"error": "corps de requête trop volumineux (limite : 1 Mo)"}), 413


@app.errorhandler(429)
def trop_de_requetes(_):
    return jsonify({"error": "trop de requêtes, réessaie dans une minute"}), 429


@app.route("/")
def index():
    return send_from_directory(".", "fem_visualizer.html")


@app.route("/api/health")
def health():
    return jsonify({"status": "ok", "runs": len(_history), "limite_par_minute": MAX_RUNS_PER_MINUTE})


@app.route("/api/mesh")
def api_mesh():
    """Maillage pour l'affichage — les paramètres géométriques passent en query string."""
    try:
        p = parse_params(request.args.to_dict())
    except ParamError as e:
        return jsonify({"error": str(e)}), 400

    mesh = Mesh(p["lx"], p["ly"], p["nx"], p["ny"])
    return jsonify({
        "lx": p["lx"], "ly": p["ly"], "nx": p["nx"], "ny": p["ny"],
        "nodes": mesh.nodes.round(6).tolist(),
        "conn": mesh.conn.tolist(),
        "boundaries": {k: v.tolist() for k, v in mesh.bnds.items()},
        "n_nodes": int(mesh.nodes.shape[0]),
        "n_elements": int(mesh.conn.shape[0]),
    })


@app.route("/api/run-simulation", methods=["POST", "OPTIONS"])
@app.route("/api/simulate", methods=["POST", "OPTIONS"])      # alias historique
def api_run():
    if request.method == "OPTIONS":
        return "", 204

    client = request.headers.get("X-Forwarded-For", request.remote_addr or "?").split(",")[0].strip()
    if rate_limited(client):
        return jsonify({"error": f"trop de calculs — maximum {MAX_RUNS_PER_MINUTE} par minute"}), 429

    try:
        p = parse_params(request.get_json(silent=True) or {})
    except ParamError as e:
        return jsonify({"error": str(e)}), 400

    key = tuple(sorted(p.items()))
    with _store_lock:
        if key in _cache:
            cached = dict(_cache[key], cached=True)
            return jsonify(cached)

    if not _solve_lock.acquire(timeout=SOLVE_TIMEOUT_S):
        return jsonify({"error": "service occupé, réessayez dans quelques secondes"}), 503
    try:
        result = verifier_fini(run_simulation(p))
    except SolverDivergence as e:
        # 422 : la requête est bien formée, c'est le modèle qui ne couvre pas ce
        # jeu de paramètres. Le message est destiné à l'utilisateur.
        return jsonify({"error": MSG_DIVERGENCE, "detail": str(e)}), 422
    except Exception:
        # Le détail part dans les logs du serveur, pas dans la réponse : un
        # message d'exception peut révéler des chemins ou la structure du code.
        app.logger.exception("échec du calcul pour %s", p)
        return jsonify({"error": "le calcul a échoué pour ces paramètres"}), 500
    finally:
        _solve_lock.release()

    with _store_lock:
        _results[result["id"]] = result
        _cache[key] = result
        _history.append({
            "id": result["id"],
            "timestamp": result["timestamp"],
            "params": result["params"],
            "summary": {k: result["summary"][k] for k in
                        ("uy_min_ballast", "uy_min_beton", "bowl_ballast",
                         "bowl_beton", "ratio_bowl", "solve_ms")},
        })
    return jsonify(result)


MAX_MATERIAUX = 6

# Bornes de la flèche de rail. Modèle analytique, donc pas de limitation de
# cadence : une évaluation coûte quelques microsecondes.
BOUNDS_RAIL = {
    "f_left_kn":     (64.0,  0.0,  250.0, "kN"),
    "f_right_kn":    (64.0,  0.0,  250.0, "kN"),
    "k_dff_kn_mm":   (27.9,  2.0,  500.0, "kN/mm"),
    "gauge":         (1.435, 0.60, 1.80,  "m"),
    "ei_MNm2":       (5.83,  0.20, 60.0,  "MN·m²"),
    "k_ballast_MPa": (40.0,  2.0,  500.0, "MPa"),
    "longueur":      (2.50,  1.00, 12.0,  "m"),
}


def parse_rail(data: dict) -> dict:
    out = {}
    for name, (default, lo, hi, unit) in BOUNDS_RAIL.items():
        try:
            val = float(data.get(name, default))
        except (TypeError, ValueError):
            raise ParamError(f"« {name} » doit être un nombre (reçu : {data.get(name)!r})")
        if not np.isfinite(val) or not (lo <= val <= hi):
            raise ParamError(f"« {name} » doit être compris entre {lo} et {hi} {unit} (reçu : {val})")
        out[name] = val
    if out["f_left_kn"] + out["f_right_kn"] <= 0:
        raise ParamError("au moins une des deux files doit être chargée")
    if out["gauge"] >= out["longueur"]:
        raise ParamError("l'écartement doit être inférieur à la longueur de la traverse")
    return out


RATIO_H_CAP = BALLAST_PARAMS["H_cap"] / BALLAST_PARAMS["E0"]

MAX_COUCHES = 4
MAX_PROFILS = 4


def construire_couches(spec, ly):
    """Transforme une liste [{id, epaisseur}] en couches du solveur.

    Les épaisseurs sont normalisées pour totaliser exactement la hauteur du
    domaine : le client raisonne en proportions, le solveur exige une partition
    exacte.
    """
    if not isinstance(spec, list) or not spec:
        raise ParamError("un profil doit être une liste non vide de couches")
    if len(spec) > MAX_COUCHES:
        raise ParamError(f"au maximum {MAX_COUCHES} couches par profil")

    epaisseurs = []
    for couche in spec:
        if not isinstance(couche, dict) or couche.get("id") not in MAT.MATERIALS:
            raise ParamError(f"couche invalide : {couche!r}")
        try:
            e = float(couche.get("epaisseur", 0.0))
        except (TypeError, ValueError):
            raise ParamError("« epaisseur » doit être un nombre")
        if not (e > 0) or not np.isfinite(e):
            raise ParamError("chaque épaisseur doit être strictement positive")
        epaisseurs.append(e)

    total = sum(epaisseurs)
    couches = []
    for couche, e in zip(spec, epaisseurs):
        meta = MAT.MATERIALS[couche["id"]]
        loi, plastique = MAT.constitutive(couche["id"])
        couches.append(Layer(loi, e / total * ly, meta["alpha"], plastique, couche["id"]))
    return couches


@app.route("/api/profile", methods=["POST", "OPTIONS"])
def api_profile():
    """Compare plusieurs empilements sous la même géométrie et le même chargement.

    C'est la question que se pose un projeteur : à hauteur d'assise donnée,
    quelle répartition entre ballast, sous-couche et sol support minimise le
    tassement.
    """
    if request.method == "OPTIONS":
        return "", 204

    client = request.headers.get("X-Forwarded-For", request.remote_addr or "?").split(",")[0].strip()
    if rate_limited(client):
        return jsonify({"error": f"trop de calculs — maximum {MAX_RUNS_PER_MINUTE} par minute"}), 429

    corps = request.get_json(silent=True) or {}
    profils = corps.get("profiles")
    if not isinstance(profils, list) or not profils:
        return jsonify({"error": "« profiles » doit être une liste non vide"}), 400
    if len(profils) > MAX_PROFILS:
        return jsonify({"error": f"au maximum {MAX_PROFILS} profils par comparaison"}), 400

    try:
        p = parse_params(corps)
        specs = [construire_couches(prof.get("layers"), p["ly"]) for prof in profils]
    except (ParamError, AttributeError) as e:
        return jsonify({"error": str(e) if isinstance(e, ParamError) else "profil mal formé"}), 400

    key = ("prof", json.dumps(profils, sort_keys=True)) + tuple(sorted(p.items()))
    with _store_lock:
        if key in _cache:
            return jsonify(dict(_cache[key], cached=True))

    if not _solve_lock.acquire(timeout=SOLVE_TIMEOUT_S):
        return jsonify({"error": "service occupé, réessayez dans quelques secondes"}), 503
    try:
        mesh = Mesh(p["lx"], p["ly"], p["nx"], p["ny"])
        t0 = time.perf_counter()
        series = []
        for prof, couches in zip(profils, specs):
            s = FEMSolver(mesh, delta_T=p["delta_T"], layers=couches)
            U, *_ = s.run(p["nsteps"], p["p"] * 1e3, p["span"])
            x, uy = s.top_profile(U)
            series.append({
                "nom": prof.get("nom") or " / ".join(c.nom for c in couches),
                "layers": [{"id": c.nom, "epaisseur": round(c.epaisseur, 4),
                            "name": MAT.MATERIALS[c.nom]["name"],
                            "color": MAT.MATERIALS[c.nom]["color"]} for c in couches],
                "uy": uy.round(3).tolist(),
                "uy_min": float(uy.min()),
                "bowl": float(uy[0] - uy.min()),
                "evp_max_pct": s.state_max("evp") * 100,
                "vm_max_kpa": float(s.stress_field(U).max()),
                "layer_map": s.layer_of_element().tolist(),
            })
        x_cm = (x * 100).round(4).tolist()
        solve_ms = round((time.perf_counter() - t0) * 1e3, 1)
        verifier_fini(series, "series")
    except SolverDivergence as e:
        return jsonify({"error": MSG_DIVERGENCE, "detail": str(e)}), 422
    except ValueError as e:
        return jsonify({"error": str(e)}), 400
    except Exception:
        app.logger.exception("échec de la comparaison de profils")
        return jsonify({"error": "le calcul a échoué pour ces paramètres"}), 500
    finally:
        _solve_lock.release()

    resultat = {
        "id": f"prof_{uuid.uuid4().hex[:8]}",
        "timestamp": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "params": p,
        "mesh": {"lx": p["lx"], "ly": p["ly"], "nx": p["nx"], "ny": p["ny"]},
        "x_cm": x_cm,
        "series": series,
        "ranking": [{"nom": s["nom"], "bowl": s["bowl"]}
                    for s in sorted(series, key=lambda s: s["bowl"])],
        "solve_ms": solve_ms,
    }
    with _store_lock:
        _cache[key] = resultat
        _results[resultat["id"]] = resultat
    return jsonify(resultat)


@app.route("/api/twin-rail", methods=["POST", "OPTIONS"])
def api_twin_rail():
    """Flexion de la traverse sous les deux files de rail (poutre de Winkler)."""
    if request.method == "OPTIONS":
        return "", 204
    try:
        p = parse_rail(request.get_json(silent=True) or {})
    except ParamError as e:
        return jsonify({"error": str(e)}), 400

    try:
        commun = dict(k_dff_kn_mm=p["k_dff_kn_mm"], gauge=p["gauge"],
                      ei=p["ei_MNm2"] * 1e6, k_ballast=p["k_ballast_MPa"] * 1e6,
                      longueur=p["longueur"])
        resultat = TB.twin_rail(f_left_kn=p["f_left_kn"], f_right_kn=p["f_right_kn"], **commun)
        # Référence : la même charge totale reportée sur une seule file, cas du
        # dévers ou du report de charge en courbe.
        mono = TB.twin_rail(f_left_kn=p["f_left_kn"] + p["f_right_kn"], f_right_kn=0.0,
                            n=len(resultat["x_cm"]), **commun)
        resultat["w_mono"] = mono["w_traverse"]
        resultat["moment_mono_kNm"] = mono["moment_kNm"]
        resultat["summary"]["mono_w_max_um"] = min(mono["w_traverse"])
        resultat["summary"]["mono_moment_max_kNm"] = mono["summary"]["moment_max_kNm"]
    except (ValueError, ZeroDivisionError, np.linalg.LinAlgError):
        app.logger.exception("échec du calcul de flexion pour %s", p)
        return jsonify({"error": "le calcul a échoué pour ces paramètres"}), 500

    return jsonify(resultat)


@app.route("/api/materials")
def api_materials():
    return jsonify({"materials": MAT.catalogue(),
                    "default": MAT.DEFAUT, "max": MAX_MATERIAUX})


@app.route("/api/compare-materials", methods=["POST", "OPTIONS"])
def api_compare_materials():
    """Calcule le profil de surface de plusieurs matériaux sous le même chargement.

    C'est la comparaison qui intéresse le concepteur : à géométrie et charge
    identiques, comment se classent les matériaux d'assise disponibles.
    """
    if request.method == "OPTIONS":
        return "", 204

    client = request.headers.get("X-Forwarded-For", request.remote_addr or "?").split(",")[0].strip()
    if rate_limited(client):
        return jsonify({"error": f"trop de calculs — maximum {MAX_RUNS_PER_MINUTE} par minute"}), 429

    corps = request.get_json(silent=True) or {}
    # `or MAT.DEFAUT` confondrait une liste vide avec un champ absent : la
    # première est une erreur du client, la seconde une requête par défaut.
    choisis = corps["materials"] if "materials" in corps else MAT.DEFAUT
    if not isinstance(choisis, list) or not choisis:
        return jsonify({"error": "« materials » doit être une liste non vide"}), 400
    if len(choisis) > MAX_MATERIAUX:
        return jsonify({"error": f"au maximum {MAX_MATERIAUX} matériaux par comparaison"}), 400
    inconnus = [m for m in choisis if m not in MAT.MATERIALS]
    if inconnus:
        return jsonify({"error": f"matériau inconnu : {', '.join(map(str, inconnus))}"}), 400

    try:
        p = parse_params(corps)
    except ParamError as e:
        return jsonify({"error": str(e)}), 400

    key = ("cmp", tuple(choisis)) + tuple(sorted(p.items()))
    with _store_lock:
        if key in _cache:
            return jsonify(dict(_cache[key], cached=True))

    if not _solve_lock.acquire(timeout=SOLVE_TIMEOUT_S):
        return jsonify({"error": "service occupé, réessayez dans quelques secondes"}), 503
    try:
        mesh = Mesh(p["lx"], p["ly"], p["nx"], p["ny"])
        t0 = time.perf_counter()
        series = []
        for mid in choisis:
            loi, plastique = MAT.constitutive(mid)
            meta = MAT.MATERIALS[mid]
            s = FEMSolver(mesh, loi, is_plastic=plastique,
                          delta_T=p["delta_T"], alpha=meta["alpha"])
            U, *_ = s.run(p["nsteps"], p["p"] * 1e3, p["span"])
            x, uy = s.top_profile(U)
            series.append({
                "id": mid,
                "name": meta["name"],
                "color": meta["color"],
                "E_MPa": meta["E"] / 1e6,
                "nu": meta["nu"],
                "alpha": meta["alpha"],
                "model": meta["model"],
                "uy": uy.round(3).tolist(),
                "uy_min": float(uy.min()),
                "uy_max": float(uy.max()),
                "bowl": float(uy[0] - uy.min()),
                "evp_max_pct": float(s.state_max("evp") * 100) if plastique else 0.0,
                "vm_max_kpa": float(s.stress_field(U).max()),
            })
        x_cm = (x * 100).round(4).tolist()
        solve_ms = round((time.perf_counter() - t0) * 1e3, 1)
        verifier_fini(series, "series")
    except SolverDivergence as e:
        return jsonify({"error": MSG_DIVERGENCE, "detail": str(e)}), 422
    except Exception:
        app.logger.exception("échec de la comparaison %s", choisis)
        return jsonify({"error": "le calcul a échoué pour ces paramètres"}), 500
    finally:
        _solve_lock.release()

    # Le classement se fait sur la cuvette : elle isole la réponse mécanique du
    # soulèvement thermique, qui dépend surtout de alpha et masque les écarts
    # de rigidité si on compare les déplacements bruts.
    ordre = sorted(series, key=lambda s: s["bowl"])
    resultat = {
        "id": f"cmp_{uuid.uuid4().hex[:8]}",
        "timestamp": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "params": p,
        "mesh": {"lx": p["lx"], "ly": p["ly"], "nx": p["nx"], "ny": p["ny"]},
        "x_cm": x_cm,
        "series": series,
        "ranking": [{"id": s["id"], "bowl": s["bowl"]} for s in ordre],
        "solve_ms": solve_ms,
    }
    with _store_lock:
        _cache[key] = resultat
        _results[resultat["id"]] = resultat
    return jsonify(resultat)


@app.route("/api/results/<run_id>")
def api_result(run_id):
    with _store_lock:
        res = _results.get(run_id)
    if res is None:
        return jsonify({"error": "calcul inconnu ou expiré"}), 404
    return jsonify(res)


@app.route("/api/results/<run_id>.csv")
def api_result_csv(run_id):
    with _store_lock:
        res = _results.get(run_id)
    if res is None:
        return jsonify({"error": "calcul inconnu ou expiré"}), 404

    lignes = ["x_cm,uy_ballast_um,uy_beton_um"]
    lignes += [f"{x},{b},{c}" for x, b, c in
               zip(res["x_cm"], res["uy_ballast"], res["uy_beton"])]
    return Response(
        "\n".join(lignes) + "\n",
        mimetype="text/csv",
        headers={"Content-Disposition": f'attachment; filename="{run_id}.csv"'},
    )


@app.route("/api/history")
def api_history():
    with _store_lock:
        return jsonify({"runs": list(reversed(_history))})


@app.route("/<path:filename>")
def static_files(filename):
    # Une route /api/ inconnue doit répondre en JSON. Sans ce garde-fou, elle
    # tombe dans le service de fichiers statiques et renvoie la page HTML : le
    # client reçoit un 200 et une page à la place de son erreur, ce qui donne
    # un « Unexpected token < » au lieu d'un message exploitable.
    if filename.startswith("api/"):
        return jsonify({"error": f"route inconnue : /{filename}"}), 404
    return send_from_directory(".", filename)


if __name__ == "__main__":
    port = int(os.environ.get("PORT", "5000"))
    print("=" * 58)
    print("  FEM Visualizer — API Flask")
    print(f"  http://localhost:{port}")
    print("  GET  /api/mesh  |  POST /api/run-simulation  |  GET /api/history")
    print("=" * 58)
    # debug=False : le débogueur Werkzeug expose une console d'exécution de code.
    app.run(host="127.0.0.1", port=port, debug=False, threaded=True)
