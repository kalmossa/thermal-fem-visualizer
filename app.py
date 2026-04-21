"""
app.py — Serveur Flask sécurisé
================================
Ce fichier fait le lien entre le frontend HTML et le solveur FEM Python.
Il reçoit les paramètres de simulation via une requête HTTP POST,
les valide, lance le calcul, et renvoie les résultats en JSON.

Pour lancer : py app.py
Puis ouvrir  : http://localhost:5000
"""

# Flask : micro-framework web Python pour créer des API REST
from flask import Flask, request, jsonify, send_from_directory

# importlib : permet d'importer un fichier .py sans lancer son bloc __main__
import importlib.util

# numpy : calcul matriciel (utilisé pour le post-traitement des résultats FEM)
import numpy as np

# os : accès au système de fichiers (pour trouver fem_ballast_beton.py)
import os

# logging : journalisation des événements dans un fichier texte
import logging

# time : horodatage des requêtes (utilisé pour le rate limiting)
import time

# defaultdict : dictionnaire avec valeur par défaut (compteur de requêtes par IP)
from collections import defaultdict


# ══════════════════════════════════════════════════════════════════════════════
# JOURNALISATION (LOGGING)
# ══════════════════════════════════════════════════════════════════════════════
# Toutes les requêtes, erreurs et avertissements sont écrits dans fem_server.log
# Le client ne voit jamais les erreurs internes brutes — sécurité + débogage facile

logging.basicConfig(
    filename='fem_server.log',    # fichier de log créé automatiquement
    level=logging.INFO,           # niveau minimum : INFO (en dessous = ignoré)
    format='%(asctime)s [%(levelname)s] %(message)s',  # format de chaque ligne
    datefmt='%Y-%m-%d %H:%M:%S'  # format de la date
)
log = logging.getLogger(__name__)  # logger associé à ce fichier


# ══════════════════════════════════════════════════════════════════════════════
# CHARGEMENT DU SOLVEUR FEM
# ══════════════════════════════════════════════════════════════════════════════
# On importe fem_ballast_beton.py comme un module Python normal,
# MAIS sans déclencher son bloc "if __name__ == '__main__'" qui lancerait
# un calcul complet au démarrage du serveur — ce qu'on ne veut pas.

def load_fem_module():
    """Charge fem_ballast_beton.py sans exécuter son bloc principal."""
    path = os.path.join(os.path.dirname(__file__), "fem_ballast_beton.py")
    spec = importlib.util.spec_from_file_location("fem_solver", path)
    mod  = importlib.util.module_from_spec(spec)
    mod.__name__ = "fem_solver"   # nom différent de "__main__" → bloque l'exécution auto
    spec.loader.exec_module(mod)  # charge le code sans lancer __main__
    return mod

# On récupère les 4 classes dont Flask a besoin pour lancer les calculs
fem           = load_fem_module()
Mesh          = fem.Mesh           # génère le maillage rectangulaire Q4
LinearElastic = fem.LinearElastic  # loi de Hooke pour le béton
CrushableCap  = fem.CrushableCap   # loi élasto-plastique pour le ballast
FEMSolver     = fem.FEMSolver      # assemblage + résolution du système K·U = F

# Maillage créé UNE SEULE FOIS au démarrage (32×14 éléments, 0.8m × 0.35m)
# Le recréer à chaque requête serait inutile et plus lent
MESH = Mesh(0.8, 0.35, 32, 14)


# ══════════════════════════════════════════════════════════════════════════════
# INITIALISATION FLASK
# ══════════════════════════════════════════════════════════════════════════════
app = Flask(__name__, static_folder=".")

# Limite la taille maximale d'une requête entrante à 1 Mo
# Empêche un attaquant d'envoyer un JSON énorme pour saturer la mémoire
app.config['MAX_CONTENT_LENGTH'] = 1 * 1024 * 1024  # 1 Mo en octets


# ══════════════════════════════════════════════════════════════════════════════
# RATE LIMITING — LIMITATION DU NOMBRE DE REQUÊTES
# ══════════════════════════════════════════════════════════════════════════════
# On autorise au maximum 10 requêtes par minute par adresse IP.
# Au-delà, le serveur refuse et retourne une erreur 429.
# Protège contre les boucles accidentelles ou les abus.

RATE_LIMIT     = 10   # nombre max de requêtes autorisées
RATE_WINDOW    = 60   # fenêtre de temps en secondes (1 minute)

# Dictionnaire : IP → liste des timestamps des requêtes récentes
# defaultdict(list) crée automatiquement une liste vide pour chaque nouvelle IP
request_counts = defaultdict(list)

def is_rate_limited(ip):
    """
    Vérifie si l'adresse IP a dépassé la limite de requêtes.
    Retourne True si bloquée, False si autorisée.
    """
    now = time.time()  # timestamp actuel en secondes

    # On supprime les timestamps plus vieux que la fenêtre de temps
    request_counts[ip] = [t for t in request_counts[ip] if now - t < RATE_WINDOW]

    # Si le compteur est plein → bloqué
    if len(request_counts[ip]) >= RATE_LIMIT:
        return True

    # Sinon on enregistre cette requête et on laisse passer
    request_counts[ip].append(now)
    return False


# ══════════════════════════════════════════════════════════════════════════════
# VALIDATION DES PARAMÈTRES
# ══════════════════════════════════════════════════════════════════════════════
# Chaque paramètre reçu est vérifié :
#   1. C'est bien un nombre (pas du texte ou du code malveillant)
#   2. Il est dans une plage physiquement réaliste
# Si un paramètre est absent, on utilise la valeur par défaut.

# Format : "nom": (valeur_min, valeur_max, valeur_par_défaut)
PARAM_RANGES = {
    "p":       (10,    1000,   300),     # pression de service en kN/m
    "delta_T": (-30,   100,    40.0),    # variation thermique en °C
    "E_b":     (10,    1000,   150),     # module de Young du ballast en MPa
    "alpha_b": (0.1,   5.0,    1.2),     # coefficient de dilatation ballast (×10⁻⁵)
    "E_c":     (1000,  100000, 30000),   # module de Young du béton en MPa
    "alpha_c": (0.1,   5.0,    1.0),     # coefficient de dilatation béton (×10⁻⁵)
}

def validate_params(data):
    """
    Valide et retourne les paramètres nettoyés.
    Retourne (dict_valide, None) si OK, ou (None, message_erreur) si problème.
    """
    validated = {}
    for key, (min_val, max_val, default) in PARAM_RANGES.items():
        raw = data.get(key, default)  # récupère la valeur ou utilise le défaut

        # Vérification 1 : est-ce un nombre ?
        try:
            val = float(raw)
        except (TypeError, ValueError):
            return None, f"Parametre '{key}' invalide : '{raw}' n'est pas un nombre."

        # Vérification 2 : est-ce dans la plage autorisée ?
        if not (min_val <= val <= max_val):
            return None, (f"Parametre '{key}' hors plage : {val} "
                          f"(valeur attendue entre {min_val} et {max_val}).")

        validated[key] = val  # paramètre validé → on le garde

    return validated, None  # tout OK, pas d'erreur


# ══════════════════════════════════════════════════════════════════════════════
# HEADERS DE SÉCURITÉ HTTP
# ══════════════════════════════════════════════════════════════════════════════
# Ces headers sont ajoutés automatiquement à TOUTES les réponses du serveur.
# Ils indiquent au navigateur comment traiter les données reçues.

@app.after_request
def add_security_headers(response):
    """Ajoute des headers de sécurité à chaque réponse HTTP."""

    # CORS : autorise le HTML ouvert localement à appeler cette API
    # Sans ça, le navigateur bloquerait les requêtes cross-origin
    response.headers["Access-Control-Allow-Origin"]  = "*"
    response.headers["Access-Control-Allow-Methods"] = "POST, GET, OPTIONS"
    response.headers["Access-Control-Allow-Headers"] = "Content-Type"

    # Empêche le navigateur de deviner le type MIME d'un fichier
    # Protection contre certaines attaques XSS par injection de contenu
    response.headers["X-Content-Type-Options"] = "nosniff"

    # Interdit l'affichage du site dans une iframe
    # Protection contre le clickjacking (piège visuel pour tromper l'utilisateur)
    response.headers["X-Frame-Options"] = "DENY"

    # Demande au navigateur d'utiliser HTTPS si disponible (pendant 1 an)
    response.headers["Strict-Transport-Security"] = "max-age=31536000"

    # Sur les routes API : désactive le cache navigateur
    # Les résultats de simulation ne doivent pas être mis en cache
    if request.path.startswith("/api/"):
        response.headers["Cache-Control"] = "no-store, no-cache, must-revalidate"

    return response


# ══════════════════════════════════════════════════════════════════════════════
# GESTION DES ERREURS HTTP
# ══════════════════════════════════════════════════════════════════════════════

@app.errorhandler(413)
def request_too_large(e):
    """Erreur 413 : requête trop volumineuse (dépassement des 1 Mo)."""
    log.warning(f"Requete trop grande depuis {request.remote_addr}")
    return jsonify({"error": "Requete trop volumineuse (maximum 1 Mo)."}), 413


# ══════════════════════════════════════════════════════════════════════════════
# ROUTES
# ══════════════════════════════════════════════════════════════════════════════

@app.route("/api/simulate", methods=["OPTIONS"])
def options():
    """Réponse aux pré-requêtes CORS (obligatoire pour les navigateurs)."""
    return "", 204  # 204 = No Content, tout va bien


@app.route("/")
def index():
    """Sert la page principale — ouvre fem_visualizer.html dans le navigateur."""
    return send_from_directory(".", "fem_visualizer.html")


@app.route("/api/simulate", methods=["POST"])
def simulate():
    """
    Route principale de l'API.
    Reçoit les paramètres en JSON, valide, calcule, renvoie les résultats.
    """
    ip = request.remote_addr  # adresse IP du client (pour logs et rate limiting)

    # ── Étape 1 : Rate limiting ───────────────────────────────────────────────
    # On refuse si l'IP a fait trop de requêtes récemment
    if is_rate_limited(ip):
        log.warning(f"Rate limit depasse pour {ip}")
        return jsonify({"error": "Trop de requetes. Patientez une minute."}), 429

    # ── Étape 2 : Parsing du JSON ─────────────────────────────────────────────
    # On tente de lire le corps de la requête comme du JSON
    # silent=True évite une exception si le JSON est mal formé
    try:
        data = request.get_json(force=True, silent=True)
        if data is None:
            # Le corps n'est pas du JSON valide
            return jsonify({"error": "Corps de requete JSON invalide."}), 400
    except Exception:
        return jsonify({"error": "Requete mal formee."}), 400

    # ── Étape 3 : Validation des paramètres ───────────────────────────────────
    # On vérifie que chaque valeur est un nombre dans la plage autorisée
    validated, err = validate_params(data)
    if err:
        log.warning(f"Parametre invalide depuis {ip} : {err}")
        return jsonify({"error": err}), 422  # 422 = Unprocessable Entity

    # ── Étape 4 : Conversion des unités ───────────────────────────────────────
    # Le frontend envoie en MPa et kN/m — le solveur FEM travaille en Pa et N/m
    p_kn_m  = validated["p"]
    delta_T = validated["delta_T"]
    E_b     = validated["E_b"]     * 1e6    # MPa → Pa
    alpha_b = validated["alpha_b"] * 1e-5   # ×10⁻⁵ → valeur réelle
    E_c     = validated["E_c"]     * 1e6    # MPa → Pa
    alpha_c = validated["alpha_c"] * 1e-5
    p_service    = p_kn_m * 1e3             # kN/m → N/m
    sleeper_span = 0.25                      # largeur traverse en m (fixe)

    try:
        # Log de la requête pour traçabilité
        log.info(
            f"Calcul FEM depuis {ip} — "
            f"p={p_kn_m}kN/m dT={delta_T}C "
            f"E_b={validated['E_b']}MPa E_c={validated['E_c']}MPa"
        )

        # ── Calcul Ballast ────────────────────────────────────────────────────
        # Loi crushable-cap avec les 8 paramètres matériau du ballast
        ballast  = CrushableCap(E_b, 0.25, 1.2, 50e3, 80e3, 1e5, 8.0, 2.5)
        solver_b = FEMSolver(MESH, ballast, is_plastic=True,
                             delta_T=delta_T, alpha=alpha_b)
        U_b, y_b, x_b, top_nodes = solver_b.run(6, p_service, sleeper_span)

        # ── Calcul Béton ──────────────────────────────────────────────────────
        # Loi de Hooke : E=E_c, nu=0.2, plane strain
        concrete = LinearElastic(E_c, 0.2)
        solver_c = FEMSolver(MESH, concrete, is_plastic=False,
                             delta_T=delta_T, alpha=alpha_c)
        U_c, y_c, x_c, _ = solver_c.run(6, p_service, sleeper_span)

        # ── Post-traitement ───────────────────────────────────────────────────
        # Tri des nœuds de surface par position x croissante
        order  = np.argsort(x_b)
        x_plot = x_b[order]

        # Extraction des déplacements verticaux (composante y = indice 2n+1)
        # Conversion m → µm (× 1e6) pour l'affichage
        uy_b = np.array([U_b[2*n+1] for n in np.array(top_nodes)[order]]) * 1e6
        uy_c = np.array([U_c[2*n+1] for n in np.array(top_nodes)[order]]) * 1e6

        max_b = float(np.max(np.abs(uy_b)))  # déplacement max ballast en µm
        max_c = float(np.max(np.abs(uy_c)))  # déplacement max béton en µm

        log.info(f"Calcul termine — ballast={max_b:.2f}um beton={max_c:.2f}um")

        # ── Réponse JSON ──────────────────────────────────────────────────────
        return jsonify({
            "x_m":         x_plot.tolist(),   # positions x en mètres
            "uy_ballast":  uy_b.tolist(),      # déplacements ballast en µm
            "uy_beton":    uy_c.tolist(),      # déplacements béton en µm
            "max_ballast": max_b,              # valeur max ballast
            "max_beton":   max_c,              # valeur max béton
            "ratio":       round(max_b / max_c, 3) if max_c else 0,  # ratio B/C
            "params": {                        # paramètres utilisés (pour l'historique)
                "p": p_kn_m, "delta_T": delta_T,
                "E_b": validated["E_b"], "alpha_b": validated["alpha_b"],
                "E_c": validated["E_c"], "alpha_c": validated["alpha_c"],
            }
        })

    except Exception as e:
        # On log l'erreur complète côté serveur (avec traceback)
        # Mais on renvoie un message générique au client — pas d'infos internes exposées
        log.error(f"Erreur solveur FEM depuis {ip} : {e}", exc_info=True)
        return jsonify({"error": "Erreur interne du serveur. Verifiez les parametres."}), 500


# ══════════════════════════════════════════════════════════════════════════════
# POINT D'ENTRÉE
# ══════════════════════════════════════════════════════════════════════════════
if __name__ == "__main__":
    print("=" * 55)
    print("  FEM Visualizer — API Flask (securise)")
    print("  http://localhost:5000")
    print("  POST /api/simulate  -> calcul FEM")
    print("  Logs : fem_server.log")
    print("=" * 55)
    # debug=False : en mode debug, Flask expose des infos internes si le code plante
    # On le désactive pour ne pas révéler l'architecture au client
    app.run(debug=False, port=5000)