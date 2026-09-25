"""Validation du modèle de traverse sur appui élastique.

Deux niveaux de contrôle :

- la solution analytique de Hetényi (poutre infinie) se vérifie exactement —
  équilibre, position des zéros, superposition ;
- le modèle aux éléments finis se vérifie contre elle dans la limite des
  poutres longues, puis sur ses propres invariants (équilibre, symétrie,
  convergence en maillage).

C'est précisément ce qu'aucune de ces vérifications n'aurait validé dans
l'approximation par cloches gaussiennes qu'il remplace.

    py -m pytest test_track_beam.py -v
"""
import math

import numpy as np
import pytest

import track_beam as TB

K_TEST = 20e6          # N/m²
EI_TEST = TB.EI_TRAVERSE
P_TEST = 64e3          # N


# ─────────────────────────────────────────────────────────────────────────────
# 1. Solution analytique : poutre infinie
# ─────────────────────────────────────────────────────────────────────────────
def test_beta_suit_sa_definition():
    assert TB.beta(K_TEST, EI_TEST) == pytest.approx((K_TEST / (4 * EI_TEST)) ** 0.25)


@pytest.mark.parametrize("k,ei", [(0, EI_TEST), (-1, EI_TEST), (K_TEST, 0)])
def test_beta_refuse_les_valeurs_non_physiques(k, ei):
    with pytest.raises(ValueError):
        TB.beta(k, ei)


def test_fleche_maximale_sous_la_charge():
    """w(0) = −P·β / 2k, et c'est bien le minimum du profil."""
    b = TB.beta(K_TEST, EI_TEST)
    attendu = -P_TEST * b / (2 * K_TEST)
    assert float(TB.deflection(0.0, P_TEST, K_TEST, EI_TEST)) == pytest.approx(attendu)

    xs = np.linspace(-8, 8, 4001)
    w = TB.deflection(xs, P_TEST, K_TEST, EI_TEST)
    assert w.min() == pytest.approx(attendu, rel=1e-9)
    assert abs(xs[np.argmin(w)]) < 1e-9


def test_equilibre_de_la_poutre_infinie():
    """∫ k·w dx = P : la réaction répartie reprend exactement la charge."""
    xs = np.linspace(-60, 60, 400001)
    reaction = -K_TEST * TB.deflection(xs, P_TEST, K_TEST, EI_TEST)
    assert np.trapezoid(reaction, xs) == pytest.approx(P_TEST, rel=1e-6)


def test_poutre_infinie_change_de_signe():
    b = TB.beta(K_TEST, EI_TEST)
    w = TB.deflection(np.linspace(0, 10 / b, 20001), P_TEST, K_TEST, EI_TEST)
    assert w.max() > 0 and w.min() < 0


def test_premier_zero_en_trois_quarts_de_pi():
    """e^(−βx)(cos βx + sin βx) s'annule en βx = 3π/4."""
    b = TB.beta(K_TEST, EI_TEST)
    xs = np.linspace(1e-6, 4 / b, 200001)
    w = TB.deflection(xs, P_TEST, K_TEST, EI_TEST)
    assert b * xs[np.argmax(w > 0)] == pytest.approx(3 * math.pi / 4, rel=1e-3)


def test_soulevement_maximal_vers_pi():
    b = TB.beta(K_TEST, EI_TEST)
    xs = np.linspace(0, 5 / b, 200001)
    w = TB.deflection(xs, P_TEST, K_TEST, EI_TEST)
    assert b * xs[np.argmax(w)] == pytest.approx(math.pi, rel=1e-3)


def test_moment_maximal_sous_la_charge():
    b = TB.beta(K_TEST, EI_TEST)
    assert float(TB.moment_infini(0.0, P_TEST, K_TEST, EI_TEST)) == pytest.approx(P_TEST / (4 * b))


# ─────────────────────────────────────────────────────────────────────────────
# 2. Éléments finis de poutre
# ─────────────────────────────────────────────────────────────────────────────
def test_poutre_longue_rejoint_hetenyi():
    """Le modèle numérique doit retrouver la solution analytique quand la poutre
    devient longue devant β⁻¹ — c'est le seul cas où les deux décrivent le même
    problème."""
    b = TB.beta(K_TEST, EI_TEST)
    longueur = 30.0 / b                          # βL = 30 : très au-delà du régime court
    poutre = TB.SleeperBeam(longueur=longueur, ei=EI_TEST, k=K_TEST, n_elem=1200)
    w, _, _ = poutre.solve([(0.0, P_TEST)])
    exact = TB.deflection(poutre.x, P_TEST, K_TEST, EI_TEST)

    assert np.max(np.abs(w - exact)) < 0.005 * np.max(np.abs(exact))


def test_equilibre_elements_finis():
    """∫ k·w dx = ΣP, à la précision de la quadrature."""
    r = TB.twin_rail(f_left_kn=50, f_right_kn=80)
    assert r["summary"]["reaction_totale_kN"] == pytest.approx(130.0, rel=2e-3)


def test_convergence_en_maillage():
    """Le raffinement doit faire converger la flèche."""
    vals = []
    for n_elem in (20, 40, 80, 160):
        poutre = TB.SleeperBeam(n_elem=n_elem)
        w, theta, _ = poutre.solve([(-0.7175, 64e3), (0.7175, 64e3)])
        vals.append(poutre.deflection_at(w, theta, -0.7175))
    e1, e2 = abs(vals[1] - vals[0]), abs(vals[3] - vals[2])
    assert e2 < e1
    assert abs(vals[3] - vals[2]) < 1e-3 * abs(vals[3])


def test_matrice_non_singuliere_sans_appui():
    """Aucun déplacement n'est bloqué : c'est la fondation qui équilibre la
    poutre. Le système doit néanmoins être inversible."""
    poutre = TB.SleeperBeam(n_elem=40)
    assert np.linalg.matrix_rank(poutre.K) == poutre.K.shape[0]


def test_charge_hors_noeud_est_repartie():
    """Une charge entre deux nœuds ne doit pas être perdue ni dupliquée."""
    poutre = TB.SleeperBeam(n_elem=37)           # maillage impair : rien ne tombe juste
    w, _, _ = poutre.solve([(0.123456, P_TEST)])
    reaction = np.trapezoid(-poutre.k * w, poutre.x)
    assert reaction == pytest.approx(P_TEST, rel=5e-3)


def test_linearite():
    poutre = TB.SleeperBeam(n_elem=60)
    w1, _, _ = poutre.solve([(0.0, P_TEST)])
    w2, _, _ = poutre.solve([(0.0, 2 * P_TEST)])
    assert np.allclose(w2, 2 * w1)


# ─────────────────────────────────────────────────────────────────────────────
# 3. Configuration traverse sous deux files
# ─────────────────────────────────────────────────────────────────────────────
def test_symetrie_sous_charges_egales():
    r = TB.twin_rail(f_left_kn=64, f_right_kn=64)
    s = r["summary"]
    assert s["w_rail_gauche_um"] == pytest.approx(s["w_rail_droit_um"], rel=1e-9)
    w = np.array(r["w_traverse"])
    assert np.allclose(w, w[::-1], atol=1e-6 * abs(w).max())


def test_superposition_des_deux_files():
    commun = dict(k_dff_kn_mm=27.9, gauge=1.435, k_ballast=40e6)
    deux = np.array(TB.twin_rail(f_left_kn=50, f_right_kn=80, **commun)["w_traverse"])
    gauche = np.array(TB.twin_rail(f_left_kn=50, f_right_kn=0, **commun)["w_left"])
    droite = np.array(TB.twin_rail(f_left_kn=0, f_right_kn=80, **commun)["w_right"])
    # les profils sont arrondis au millième de µm à la sérialisation
    assert np.allclose(deux, gauche + droite, atol=3e-3)


def test_moment_positif_au_siege_negatif_au_centre():
    """Signature de flexion d'une traverse : elle porte sur le ballast entre les
    rails, d'où l'inversion du moment au centre."""
    s = TB.twin_rail()["summary"]
    assert s["moment_siege_kNm"] > 0
    assert s["moment_centre_kNm"] < 0


def test_traverse_reelle_ne_se_souleve_pas():
    """Avec βL ≈ 2,9, la traverse s'enfonce partout : le centre moins que les
    sièges, mais jamais au-dessus de zéro. Les lobes de soulèvement de la
    solution infinie demandent une poutre bien plus longue."""
    s = TB.twin_rail()["summary"]
    assert s["beta_L"] < 4.0
    assert s["souleve"] is False
    assert s["w_haut_um"] < 0
    assert s["differentiel_siege_centre_um"] < 0, "le centre doit être moins enfoncé"


def test_le_soulevement_apparait_pour_une_poutre_longue():
    """Transition entre les deux régimes — le modèle doit produire les deux."""
    court = TB.twin_rail(longueur=2.5, n=401)["summary"]
    long_ = TB.twin_rail(longueur=8.0, n=401)["summary"]
    assert court["souleve"] is False
    assert long_["souleve"] is True


def test_raideur_du_ballast_change_la_forme():
    """β dépend de k : un ballast plus raide resserre le bassin de flexion.

    L'approximation gaussienne remplacée fixait la largeur à 0,55 fois la
    longueur du siège ; ni k ni EI n'y agissaient sur la forme.
    """
    mou = TB.twin_rail(k_ballast=10e6)["summary"]
    dur = TB.twin_rail(k_ballast=200e6)["summary"]
    assert dur["longueur_caracteristique_m"] < mou["longueur_caracteristique_m"]
    assert abs(dur["w_traverse_gauche_um"]) < abs(mou["w_traverse_gauche_um"])


def test_rigidite_de_la_traverse_elargit_le_bassin():
    souple = TB.twin_rail(ei=TB.EI_TRAVERSE / 8)["summary"]
    raide = TB.twin_rail(ei=TB.EI_TRAVERSE * 8)["summary"]
    assert raide["longueur_caracteristique_m"] > souple["longueur_caracteristique_m"]


def test_semelle_agit_sur_le_rail_pas_sur_la_traverse():
    """La semelle est un ressort ponctuel entre rail et traverse : elle change
    la flèche du rail sans toucher celle de la traverse."""
    souple = TB.twin_rail(k_dff_kn_mm=10)["summary"]
    raide = TB.twin_rail(k_dff_kn_mm=300)["summary"]
    assert souple["w_traverse_gauche_um"] == pytest.approx(raide["w_traverse_gauche_um"])
    assert abs(souple["w_rail_gauche_um"]) > abs(raide["w_rail_gauche_um"])
    assert abs(souple["ecrasement_gauche_um"]) > abs(raide["ecrasement_gauche_um"])


def test_ecartement_trop_grand_refuse():
    with pytest.raises(ValueError):
        TB.twin_rail(gauge=3.0, longueur=2.5)


def test_ordres_de_grandeur_ferroviaires():
    """Garde-fou : une charge de roue courante doit donner des valeurs de voie."""
    s = TB.twin_rail(f_left_kn=64, f_right_kn=64)["summary"]
    assert 0.4 < s["longueur_caracteristique_m"] < 2.0
    assert 0.2 < abs(s["w_traverse_gauche_um"]) / 1000 < 5.0       # mm
    assert 1.0 < s["moment_max_kNm"] < 40.0                        # kN·m


def test_profils_de_meme_longueur():
    r = TB.twin_rail(n=121)
    n = len(r["x_cm"])
    for cle in ("w_traverse", "w_left", "w_right", "reaction_kN_m", "moment_kNm"):
        assert len(r[cle]) == n
