"""Tests du profil multicouche et du chargement cyclique.

Le multicouche est la brique qui fait passer le comparateur d'un classement de
matériaux à un outil de dimensionnement : une plateforme réelle est stratifiée,
et c'est l'épaisseur de chaque couche qui constitue le levier de conception.

Le chargement cyclique est livré avec son domaine de validité explicite — le
modèle de cap ne supporte pas une décharge complète, et les tests le vérifient
plutôt que de le passer sous silence.

    py -m pytest test_profil.py -v
"""
import numpy as np
import pytest

import materials as MAT
from fem_ballast_beton import (
    BALLAST_PARAMS, CrushableCap, FEMSolver, Layer, LinearElastic, Mesh,
)


def couches(specs):
    """specs : [(identifiant de matériau, épaisseur)] du haut vers le bas."""
    out = []
    for mid, ep in specs:
        loi, plastique = MAT.constitutive(mid)
        out.append(Layer(loi, ep, MAT.MATERIALS[mid]["alpha"], plastique, mid))
    return out


def resoudre(specs, ly=0.35, nx=24, ny=10, dT=40.0):
    s = FEMSolver(Mesh(0.8, ly, nx, ny), delta_T=dT, layers=couches(specs))
    U, *_ = s.run(6, 300e3, 0.25)
    return s, U, s.top_profile(U)[1]


# ─────────────────────────────────────────────────────────────────────────────
# Cohérence du découpage
# ─────────────────────────────────────────────────────────────────────────────
def test_une_seule_couche_egale_le_cas_homogene():
    """Non-régression : le chemin multicouche ne doit rien changer quand il n'y
    a qu'une couche."""
    mesh = Mesh(0.8, 0.35, 24, 10)
    homogene = FEMSolver(mesh, CrushableCap(**BALLAST_PARAMS), is_plastic=True,
                         delta_T=40.0, alpha=1.2e-5)
    U1, *_ = homogene.run(6, 300e3, 0.25)

    stratifie = FEMSolver(mesh, delta_T=40.0, layers=couches([("ballast", 0.35)]))
    U2, *_ = stratifie.run(6, 300e3, 0.25)

    assert np.allclose(U1, U2, atol=1e-14)


def test_les_couches_partitionnent_le_maillage():
    s = FEMSolver(Mesh(0.8, 0.35, 24, 10), delta_T=0.0,
                  layers=couches([("ballast", 0.20), ("grave_bitume", 0.15)]))
    tous = np.concatenate([c.elements for c in s.layers])
    assert np.array_equal(np.sort(tous), np.arange(s.nelem))


def test_epaisseurs_incoherentes_refusees():
    with pytest.raises(ValueError, match="totalisent"):
        FEMSolver(Mesh(0.8, 0.35, 24, 10), delta_T=0.0,
                  layers=couches([("ballast", 0.20), ("grave_bitume", 0.20)]))


def test_couche_trop_fine_pour_le_maillage_refusee():
    """Une couche sans aucun élément est une erreur de maillage, pas un cas à
    traiter silencieusement."""
    with pytest.raises(ValueError, match="aucun élément"):
        FEMSolver(Mesh(0.8, 0.35, 24, 4), delta_T=0.0,      # éléments de 8,75 cm
                  layers=couches([("ballast", 0.34), ("grave_bitume", 0.01)]))


def test_epaisseur_negative_refusee():
    with pytest.raises(ValueError):
        Layer(LinearElastic(30e9, 0.2), -0.1)


def test_solveur_sans_loi_ni_couches_refuse():
    with pytest.raises(ValueError):
        FEMSolver(Mesh(0.8, 0.35, 8, 4))


def test_affectation_par_ordonnee():
    """Les couches sont données du haut vers le bas."""
    s = FEMSolver(Mesh(0.8, 0.30, 16, 6), delta_T=0.0,      # éléments de 5 cm
                  layers=couches([("ballast", 0.10), ("sol_support", 0.20)]))
    yc = s.nodes[s.conn, 1].mean(axis=1)
    assert yc[s.layers[0].elements].min() > 0.20 - 1e-9, "la 1re couche est en haut"
    assert yc[s.layers[1].elements].max() < 0.20 + 1e-9


# ─────────────────────────────────────────────────────────────────────────────
# Comportement mécanique du profil
# ─────────────────────────────────────────────────────────────────────────────
def test_une_sous_couche_raide_reduit_le_tassement():
    _, _, seul = resoudre([("ballast", 0.35)])
    _, _, avec = resoudre([("ballast", 0.25), ("grave_bitume", 0.10)])
    assert (avec[0] - avec.min()) < 0.9 * (seul[0] - seul.min())


def test_une_sous_couche_molle_n_apporte_rien():
    _, _, seul = resoudre([("ballast", 0.35)])
    _, _, mou = resoudre([("ballast", 0.25), ("sol_support", 0.10)])
    assert (mou[0] - mou.min()) > 0.9 * (seul[0] - seul.min())


def test_amincir_le_ballast_reduit_le_tassement():
    """Contre-intuitif mais conforme : le ballast est le maillon plastique.
    Remplacer une partie de son épaisseur par une grave-bitume réduit le
    tassement — c'est la justification des sous-couches bitumineuses.
    """
    _, _, epais = resoudre([("ballast", 0.25), ("grave_bitume", 0.10)])
    _, _, mince = resoudre([("ballast", 0.15), ("grave_bitume", 0.20)])
    assert (mince[0] - mince.min()) < (epais[0] - epais.min())


def test_chaque_couche_a_sa_dilatation():
    """La grave-bitume dilate 2,5 fois plus que la grave non traitée : à ΔT
    identique, le soulèvement dépend de l'empilement."""
    _, _, a = resoudre([("ballast", 0.25), ("grave_bitume", 0.10)])
    _, _, b = resoudre([("ballast", 0.25), ("grave_non_traitee", 0.10)])
    assert a[0] > b[0], "la couche la plus dilatante soulève davantage"


def test_champs_definis_sur_tout_le_maillage():
    s, U, _ = resoudre([("ballast", 0.25), ("grave_bitume", 0.10)])
    for champ in (s.stress_field(U), s.gauss_field("evp")):
        assert champ.shape == (s.nelem,)
        assert np.all(np.isfinite(champ))

    evp = s.gauss_field("evp")
    assert np.all(evp[s.layers[1].elements] == 0.0), "couche élastique : pas de plasticité"
    assert evp[s.layers[0].elements].max() > 0.0


def test_layer_of_element():
    s, _, _ = resoudre([("ballast", 0.25), ("grave_bitume", 0.10)])
    idx = s.layer_of_element()
    assert set(np.unique(idx)) == {0, 1}
    assert np.all(idx[s.layers[0].elements] == 0)


# ─────────────────────────────────────────────────────────────────────────────
# Chargement cyclique
# ─────────────────────────────────────────────────────────────────────────────
def test_adaptation_apres_le_premier_cycle():
    """Le cap ne s'écrouit que dans le sens de la compaction et ne se dégrade
    jamais : il finit par envelopper l'état atteint, et les cycles suivants sont
    élastiques. C'est l'adaptation, ou shakedown.
    """
    s = FEMSolver(Mesh(0.8, 0.35, 20, 8), CrushableCap(**BALLAST_PARAMS),
                  is_plastic=True, delta_T=40.0, alpha=1.2e-5)
    _, hist = s.run_cyclic(n_cycles=4, p_total=300e3)

    maxima = [h["uy_centre_um"] for h in hist if h["phase"] == "charge"]
    assert abs(maxima[-1] - maxima[1]) < 0.01 * abs(maxima[1]), f"pas d'adaptation : {maxima}"

    evp = [h["evp_max_pct"] for h in hist]
    assert max(evp) - min(evp) < 1e-9, "plus aucune plasticité après le 1er cycle"


def test_la_decharge_laisse_un_tassement_residuel():
    s = FEMSolver(Mesh(0.8, 0.35, 20, 8), CrushableCap(**BALLAST_PARAMS),
                  is_plastic=True, delta_T=40.0, alpha=1.2e-5)
    _, hist = s.run_cyclic(n_cycles=2, p_total=300e3)
    en_charge = next(h for h in hist if h["phase"] == "charge")
    decharge = next(h for h in hist if h["phase"] == "decharge")
    assert decharge["uy_centre_um"] < 0, "la surface reste enfoncée après décharge"
    assert abs(decharge["uy_centre_um"]) < abs(en_charge["uy_centre_um"])


@pytest.mark.parametrize("lam_min", [0.0, 0.1, 0.34, -0.2, 1.0])
def test_decharge_complete_refusee(lam_min):
    """Hors du domaine de validité, le modèle lève plutôt que de produire un
    résultat : l'ellipse est fermée, une décharge profonde la franchit par le
    bas, et le retour radial fait résister le matériau à la décharge.
    """
    s = FEMSolver(Mesh(0.8, 0.35, 12, 6), CrushableCap(**BALLAST_PARAMS),
                  is_plastic=True)
    with pytest.raises(ValueError):
        s.run_cyclic(n_cycles=1, lam_min=lam_min)


def test_ecrouissage_du_cap_est_monotone():
    """La compaction des grains est irréversible : pc ne doit jamais décroître."""
    mat = CrushableCap(**BALLAST_PARAMS)
    state = mat.init_state(1)
    pcs = []
    for eps in (-1e-3, -2e-3, -3e-3, -2e-3, -1e-3, -5e-4):
        _, _, state = mat.stress_update(np.array([[eps, eps, 0.0]]),
                                        np.zeros((1, 3)), state)
        pcs.append(float(state["pc"][0]))
    assert all(b >= a - 1e-9 for a, b in zip(pcs, pcs[1:])), f"pc décroît : {pcs}"
