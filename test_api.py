"""Tests de l'API REST et de la bibliothèque de matériaux.

Complète `test_fem.py`, qui valide la physique du solveur : ici on vérifie le
contrat de l'API — codes de retour, forme des réponses, bornes des paramètres —
et la cohérence de la bibliothèque de matériaux.

    py -m pytest test_api.py -v
"""
import json

import pytest

import app as API
import materials as MAT
from fem_ballast_beton import CrushableCap, LinearElastic


@pytest.fixture
def client():
    """Client de test, avec les compteurs de cadence remis à zéro.

    La limitation est un état de processus : sans ce nettoyage, une suite un peu
    longue finirait par se faire répondre 429 par sa propre protection.
    """
    API._rate_hits.clear()
    API._cache.clear()
    API._results.clear()
    API._history.clear()
    return API.app.test_client()


# ─────────────────────────────────────────────────────────────────────────────
# Bibliothèque de matériaux
# ─────────────────────────────────────────────────────────────────────────────
def test_catalogue_est_bilingue():
    """Le front bascule FR/EN sans redemander le catalogue : les deux langues
    doivent donc arriver ensemble."""
    for entree in MAT.catalogue():
        assert set(entree["name"]) >= {"en", "fr"}
        assert set(entree["note"]) >= {"en", "fr"}
        assert entree["name"]["en"] and entree["name"]["fr"]


def test_catalogue_parametres_physiques_plausibles():
    for entree in MAT.catalogue():
        assert entree["E_MPa"] > 0
        assert 0.0 <= entree["nu"] < 0.5, "nu = 0,5 rend la loi de Hooke singulière"
        assert 0.0 <= entree["alpha"] <= 1e-4
        assert entree["model"] in ("elastic", "cap")
        assert entree["color"].startswith("#")


def test_defauts_presents_dans_la_bibliotheque():
    assert set(MAT.DEFAUT) <= set(MAT.MATERIALS)


def test_constitutive_instancie_la_bonne_loi():
    loi, plastique = MAT.constitutive("ballast")
    assert isinstance(loi, CrushableCap) and plastique is True

    loi, plastique = MAT.constitutive("beton")
    assert isinstance(loi, LinearElastic) and plastique is False


def test_constitutive_respecte_le_module_impose():
    loi, _ = MAT.constitutive("beton", E_override=12e9)
    assert loi.E == pytest.approx(12e9)


# ─────────────────────────────────────────────────────────────────────────────
# Endpoints de base
# ─────────────────────────────────────────────────────────────────────────────
def test_health(client):
    d = client.get("/api/health").get_json()
    assert d["status"] == "ok"


def test_mesh(client):
    d = client.get("/api/mesh").get_json()
    assert d["n_nodes"] == (d["nx"] + 1) * (d["ny"] + 1)
    assert d["n_elements"] == d["nx"] * d["ny"]
    assert len(d["nodes"]) == d["n_nodes"]
    assert len(d["conn"]) == d["n_elements"]
    assert set(d["boundaries"]) == {"left", "right", "top", "bottom"}


def test_mesh_refuse_une_geometrie_hors_bornes(client):
    assert client.get("/api/mesh?nx=9999").status_code == 400


# ─────────────────────────────────────────────────────────────────────────────
# Simulation
# ─────────────────────────────────────────────────────────────────────────────
def test_run_simulation_nominal(client):
    d = client.post("/api/run-simulation", json={}).get_json()
    s = d["summary"]

    assert len(d["x_cm"]) == len(d["uy_ballast"]) == len(d["uy_beton"])
    assert s["bowl_ballast"] > 100 * s["bowl_beton"], "le ballast s'enfonce bien davantage"
    assert s["uy_min_beton"] > 0, "le béton reste soulevé par la dilatation"
    assert s["uy_min_ballast"] < 0, "le ballast s'enfonce sous la traverse"
    assert s["evp_max_pct"] > 0, "la plasticité doit s'activer"

    n_elem = d["mesh"]["n_elements"]
    assert len(d["fields"]["vm_ballast"]) == n_elem
    assert len(d["fields"]["uy_nodal_ballast"]) == d["mesh"]["n_nodes"]


@pytest.mark.parametrize("charge_utile,attendu", [
    ({"E_b": 0}, "E_b"),
    ({"E_b": "abc"}, "E_b"),
    ({"delta_T": 500}, "delta_T"),
    ({"nx": 9999}, "nx"),
    ({"span": 0.59, "lx": 0.3}, "span"),
])
def test_parametres_hors_bornes(client, charge_utile, attendu):
    r = client.post("/api/run-simulation", json=charge_utile)
    assert r.status_code == 400
    assert attendu in r.get_json()["error"]


def test_erreur_ne_fuite_pas_la_stack(client):
    """Un message d'exception peut révéler des chemins ou la structure du code."""
    r = client.post("/api/run-simulation", json={"E_b": "abc"})
    message = r.get_json()["error"]
    assert "Traceback" not in message and ".py" not in message


def test_mise_en_cache(client):
    premier = client.post("/api/run-simulation", json={"p": 250}).get_json()
    second = client.post("/api/run-simulation", json={"p": 250}).get_json()
    assert second.get("cached") is True
    assert second["uy_ballast"] == premier["uy_ballast"]


def test_relecture_et_export_csv(client):
    run_id = client.post("/api/run-simulation", json={}).get_json()["id"]
    assert client.get(f"/api/results/{run_id}").status_code == 200

    csv = client.get(f"/api/results/{run_id}.csv")
    assert csv.status_code == 200
    lignes = csv.get_data(as_text=True).strip().splitlines()
    assert lignes[0] == "x_cm,uy_ballast_um,uy_beton_um"
    assert len(lignes) == 34                      # en-tête + 33 nœuds de surface


def test_resultat_inconnu(client):
    assert client.get("/api/results/inexistant").status_code == 404


def test_historique(client):
    client.post("/api/run-simulation", json={"p": 100})
    client.post("/api/run-simulation", json={"p": 200})
    runs = client.get("/api/history").get_json()["runs"]
    assert len(runs) == 2
    assert runs[0]["params"]["p"] == 200, "le plus récent d'abord"


# ─────────────────────────────────────────────────────────────────────────────
# Comparaison de matériaux
# ─────────────────────────────────────────────────────────────────────────────
def test_materials_catalogue(client):
    d = client.get("/api/materials").get_json()
    assert len(d["materials"]) == len(MAT.MATERIALS)
    assert d["default"] == MAT.DEFAUT
    assert d["max"] >= 2


def test_compare_materials(client):
    choisis = ["beton", "ballast", "sol_support"]
    d = client.post("/api/compare-materials", json={"materials": choisis}).get_json()

    assert [s["id"] for s in d["series"]] == choisis
    for s in d["series"]:
        assert len(s["uy"]) == len(d["x_cm"])
        assert s["bowl"] >= 0

    cuvettes = {s["id"]: s["bowl"] for s in d["series"]}
    assert cuvettes["beton"] < cuvettes["sol_support"] < cuvettes["ballast"], (
        "le classement doit suivre la souplesse, la plasticité du ballast "
        "l'emportant sur la faible raideur du sol"
    )


def test_compare_materials_classement_trie(client):
    d = client.post("/api/compare-materials",
                    json={"materials": ["ballast", "beton", "grave_bitume"]}).get_json()
    cuvettes = [r["bowl"] for r in d["ranking"]]
    assert cuvettes == sorted(cuvettes)


def test_compare_materials_seul_le_materiau_change(client):
    """Deux matériaux identiques doivent donner exactement le même profil : c'est
    le garde-fou contre un état qui fuirait d'un calcul au suivant."""
    d = client.post("/api/compare-materials",
                    json={"materials": ["beton", "beton"]}).get_json()
    a, b = d["series"]
    assert a["uy"] == b["uy"]


@pytest.mark.parametrize("charge_utile", [
    {"materials": []},
    {"materials": "ballast"},
    {"materials": ["materiau_qui_nexiste_pas"]},
    {"materials": ["ballast"] * 9},
    {"materials": ["ballast"], "delta_T": 999},
])
def test_compare_materials_entrees_invalides(client, charge_utile):
    assert client.post("/api/compare-materials", json=charge_utile).status_code == 400


def test_compare_materials_sans_liste_utilise_les_defauts(client):
    d = client.post("/api/compare-materials", json={}).get_json()
    assert [s["id"] for s in d["series"]] == MAT.DEFAUT


# ─────────────────────────────────────────────────────────────────────────────
# Limitation de cadence
# ─────────────────────────────────────────────────────────────────────────────
def test_limitation_de_cadence(client, monkeypatch):
    """Chaque calcul coûte ~0,7 s de CPU : sans plafond, une boucle curl suffit
    à coucher l'instance."""
    monkeypatch.setattr(API, "MAX_RUNS_PER_MINUTE", 3)
    API._rate_hits.clear()
    for _ in range(3):
        assert not API.rate_limited("1.2.3.4")
    assert API.rate_limited("1.2.3.4")
    assert not API.rate_limited("5.6.7.8"), "le plafond est par client"


# ─────────────────────────────────────────────────────────────────────────────
# Flexion de la traverse
# ─────────────────────────────────────────────────────────────────────────────
def test_twin_rail_nominal(client):
    d = client.post("/api/twin-rail", json={}).get_json()
    s = d["summary"]

    n = len(d["x_cm"])
    for cle in ("w_traverse", "w_left", "w_right", "reaction_kN_m", "moment_kNm", "w_mono"):
        assert len(d[cle]) == n

    assert s["reaction_totale_kN"] == pytest.approx(128.0, rel=2e-3), "équilibre global"
    assert s["moment_siege_kNm"] > 0 > s["moment_centre_kNm"]
    assert s["w_traverse_gauche_um"] < 0
    assert abs(s["w_rail_gauche_um"]) > abs(s["w_traverse_gauche_um"]), (
        "le rail descend de la flèche de la traverse plus l'écrasement de la semelle"
    )


def test_twin_rail_charge_dissymetrique(client):
    d = client.post("/api/twin-rail", json={"f_left_kn": 100, "f_right_kn": 20}).get_json()
    s = d["summary"]
    assert abs(s["w_traverse_gauche_um"]) > abs(s["w_traverse_droit_um"])
    assert s["reaction_totale_kN"] == pytest.approx(120.0, rel=2e-3)


def test_twin_rail_report_sur_une_file_est_penalisant(client):
    """Concentrer la même charge sur une seule file creuse davantage — c'est ce
    que montre la comparaison mono / bi-rail."""
    d = client.post("/api/twin-rail", json={}).get_json()["summary"]
    assert abs(d["mono_w_max_um"]) > abs(d["w_traverse_gauche_um"])
    assert d["mono_moment_max_kNm"] > d["moment_max_kNm"]


@pytest.mark.parametrize("charge_utile,attendu", [
    ({"f_left_kn": -1}, "f_left_kn"),
    ({"k_dff_kn_mm": 0}, "k_dff_kn_mm"),
    ({"gauge": 5}, "gauge"),
    ({"longueur": 0.2}, "longueur"),
    ({"ei_MNm2": "abc"}, "ei_MNm2"),
    ({"f_left_kn": 0, "f_right_kn": 0}, "chargée"),
    ({"gauge": 1.7, "longueur": 1.2}, "écartement"),
])
def test_twin_rail_entrees_invalides(client, charge_utile, attendu):
    r = client.post("/api/twin-rail", json=charge_utile)
    assert r.status_code == 400
    assert attendu in r.get_json()["error"]


# ─────────────────────────────────────────────────────────────────────────────
# Profils multicouches
# ─────────────────────────────────────────────────────────────────────────────
def profil(nom, *couches):
    return {"nom": nom, "layers": [{"id": m, "epaisseur": e} for m, e in couches]}


def test_profile_nominal(client):
    d = client.post("/api/profile", json={"profiles": [
        profil("ballast seul", ("ballast", 0.35)),
        profil("ballast + GB", ("ballast", 0.25), ("grave_bitume", 0.10)),
    ]}).get_json()

    assert len(d["series"]) == 2
    for s in d["series"]:
        assert len(s["uy"]) == len(d["x_cm"])
        assert len(s["layer_map"]) == d["mesh"]["nx"] * d["mesh"]["ny"]
        assert s["bowl"] > 0


def test_profile_sous_couche_raide_gagne(client):
    """Le classement doit refléter le gain d'une sous-couche rigide."""
    d = client.post("/api/profile", json={"profiles": [
        profil("seul", ("ballast", 0.35)),
        profil("avec GB", ("ballast", 0.25), ("grave_bitume", 0.10)),
    ]}).get_json()
    cuvettes = {s["nom"]: s["bowl"] for s in d["series"]}
    assert cuvettes["avec GB"] < cuvettes["seul"]
    assert d["ranking"][0]["nom"] == "avec GB"


def test_profile_epaisseurs_normalisees(client):
    """Le client raisonne en proportions ; le solveur exige une partition
    exacte du domaine. Doubler toutes les épaisseurs ne doit rien changer."""
    a = client.post("/api/profile", json={"profiles": [
        profil("a", ("ballast", 0.25), ("grave_bitume", 0.10))]}).get_json()
    b = client.post("/api/profile", json={"profiles": [
        profil("a", ("ballast", 0.50), ("grave_bitume", 0.20))]}).get_json()
    assert a["series"][0]["uy"] == b["series"][0]["uy"]


def test_profile_carte_des_couches(client):
    d = client.post("/api/profile", json={"profiles": [
        profil("deux", ("ballast", 0.25), ("grave_bitume", 0.10))]}).get_json()
    carte = d["series"][0]["layer_map"]
    assert set(carte) == {0, 1}


@pytest.mark.parametrize("charge_utile", [
    {"profiles": []},
    {"profiles": "ballast"},
    {"profiles": [{"layers": []}]},
    {"profiles": [{"layers": [{"id": "inconnu", "epaisseur": 1}]}]},
    {"profiles": [{"layers": [{"id": "ballast", "epaisseur": 0}]}]},
    {"profiles": [{"layers": [{"id": "ballast", "epaisseur": "x"}]}]},
    {"profiles": [{"layers": [{"id": "ballast", "epaisseur": 1}]}] * 9},
    {"profiles": [{"layers": [{"id": "ballast", "epaisseur": 1}] * 9}]},
])
def test_profile_entrees_invalides(client, charge_utile):
    assert client.post("/api/profile", json=charge_utile).status_code == 400


def test_profile_couche_trop_fine_renvoie_400(client):
    """Une couche sans élément est une erreur d'entrée, pas une erreur serveur."""
    r = client.post("/api/profile", json={
        "ny": 4,
        "profiles": [profil("fine", ("ballast", 0.99), ("grave_bitume", 0.01))]})
    assert r.status_code == 400
    assert "élément" in r.get_json()["error"]


def test_route_api_inconnue_repond_en_json(client):
    """Sans garde-fou, une route /api/ inconnue tombe dans le service de
    fichiers statiques et renvoie la page HTML avec un code 200."""
    r = client.get("/api/route-qui-nexiste-pas")
    assert r.status_code == 404
    assert r.is_json
    assert "inconnue" in r.get_json()["error"]


# ─────────────────────────────────────────────────────────────────────────────
# Robustesse de la sérialisation
# ─────────────────────────────────────────────────────────────────────────────
@pytest.mark.parametrize("charge_utile", [
    {},
    {"E_b": 2000, "p": 2000, "delta_T": 80},
    {"E_b": 10, "p": 2000, "delta_T": -60, "ly": 0.6},
    {"E_b": 430, "alpha_b": 2.6e-5, "E_c": 50500, "alpha_c": 1.7e-5,
     "p": 225, "delta_T": 40, "span": 0.39, "ly": 0.60, "lx": 1.40},
])
def test_reponse_toujours_json_valide(client, charge_utile):
    """Une réponse ne doit jamais contenir NaN ni Infinity.

    `jsonify` les écrit tels quels, ce qui n'est pas du JSON : le client reçoit
    un 200 et un corps que `response.json()` refuse, avec un « Unexpected token
    N » pour tout diagnostic. Une combinaison hors domaine doit donner un 422
    explicite, pas un succès illisible.
    """
    r = client.post("/api/run-simulation", json=charge_utile)
    corps = r.get_data(as_text=True)
    assert "NaN" not in corps and "Infinity" not in corps
    assert r.status_code in (200, 422)
    json.loads(corps)                       # doit être analysable sans erreur
    if r.status_code == 422:
        assert r.get_json()["error"]


def test_divergence_est_signalee(client):
    """Le solveur doit refuser une décharge complète plutôt que diverger."""
    from fem_ballast_beton import BALLAST_PARAMS, CrushableCap, FEMSolver, Mesh, SolverDivergence

    s = FEMSolver(Mesh(0.8, 0.35, 12, 6), CrushableCap(**BALLAST_PARAMS), is_plastic=True)
    with pytest.raises((SolverDivergence, ValueError)):
        s.run_cyclic(n_cycles=1, lam_min=0.0)


def test_valeur_non_finie_refusee_a_la_serialisation():
    """Le garde-fou doit voir un NaN enfoui dans une structure imbriquée."""
    import app as API
    from fem_ballast_beton import SolverDivergence

    API.verifier_fini({"a": [1.0, 2.0], "b": {"c": 3.0}})        # passe
    with pytest.raises(SolverDivergence):
        API.verifier_fini({"a": [1.0, float("nan")]})
    with pytest.raises(SolverDivergence):
        API.verifier_fini({"a": {"b": [{"c": float("inf")}]}})
