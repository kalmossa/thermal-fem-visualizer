"""Tests de validation du solveur FEM.

Trois familles :

1. **Vérification** — le code résout-il correctement les équations ? Patch test,
   intégration exacte, comparaison de l'assemblage vectorisé à une boucle naïve.
2. **Validation** — les équations sont-elles les bonnes ? Confrontation à des
   solutions analytiques (dilatation thermique contrainte, linéarité).
3. **Robustesse numérique** — la réponse est-elle indépendante de la
   discrétisation ? Indépendance au nombre d'incréments, convergence en maillage.

    py -m pytest test_fem.py -v
"""
import numpy as np
import pytest
from scipy.sparse.linalg import spsolve

from fem_ballast_beton import (
    BALLAST_PARAMS,
    CrushableCap,
    FEMSolver,
    LinearElastic,
    Mesh,
    M_VOL,
    B_matrix,
    gauss_points,
    shape_Q4,
)

E_BETON, NU_BETON = 30e9, 0.2
ALPHA_BETON = 1.0e-5


def moduli(E, nu):
    return E / (3 * (1 - 2 * nu)), E / (2 * (1 + nu))


# ─────────────────────────────────────────────────────────────────────────────
# 1. Vérification
# ─────────────────────────────────────────────────────────────────────────────
def test_maillage_compte_noeuds_et_elements():
    mesh = Mesh(0.8, 0.35, 32, 14)
    assert mesh.nodes.shape == (33 * 15, 2)
    assert mesh.conn.shape == (32 * 14, 4)
    assert len(mesh.bnds["top"]) == 33
    assert len(mesh.bnds["left"]) == 15


def test_quadrature_integre_exactement_l_aire():
    """Σ detJ·w sur tous les points de Gauss doit rendre l'aire du domaine."""
    mesh = Mesh(0.8, 0.35, 32, 14)
    s = FEMSolver(mesh, LinearElastic(E_BETON, NU_BETON))
    assert s.dV.sum() == pytest.approx(0.8 * 0.35, rel=1e-14)


def test_fonctions_de_forme():
    """Partition de l'unité, et dérivées de somme nulle."""
    for xi, eta, _ in gauss_points():
        N, dN = shape_Q4(xi, eta)
        assert N.sum() == pytest.approx(1.0)
        assert dN.sum(axis=0) == pytest.approx(np.zeros(2))


def test_assemblage_vectorise_egale_boucle_naive():
    """La version einsum doit rendre exactement la même matrice que l'assemblage
    élément par élément — c'est le garde-fou de la vectorisation."""
    mesh = Mesh(0.2, 0.1, 5, 3)
    mat = LinearElastic(210e9, 0.3)
    s = FEMSolver(mesh, mat)
    K_vec, _, _ = s._assemble(np.zeros(s.ndof), 0.0, s._init_state())

    D = mat.D()
    K_ref = np.zeros((s.ndof, s.ndof))
    for e, enodes in enumerate(mesh.conn):
        Xe = mesh.nodes[enodes, :]
        Ke = np.zeros((8, 8))
        for xi, eta, w in gauss_points():
            _, dN_dxi = shape_Q4(xi, eta)
            J = dN_dxi.T @ Xe
            detJ = np.linalg.det(J)
            dN_dx = dN_dxi @ np.linalg.inv(J).T
            Bm = B_matrix(dN_dx)
            Ke += (Bm.T @ D @ Bm) * detJ * w
        edofs = np.empty(8, dtype=int)
        edofs[0::2] = 2 * enodes
        edofs[1::2] = 2 * enodes + 1
        K_ref[np.ix_(edofs, edofs)] += Ke

    # Les termes de K valent ici ~5e11 : la tolérance se juge en relatif, au
    # niveau de l'epsilon machine, pas en absolu.
    ecart = np.abs(K_vec.toarray() - K_ref).max()
    assert ecart <= 1e-13 * np.abs(K_ref).max(), f"écart relatif {ecart / np.abs(K_ref).max():.2e}"


def test_patch_test():
    """Patch test : un champ de déplacement linéaire imposé au bord doit être
    reproduit exactement à l'intérieur. C'est la condition de convergence d'un
    élément fini — si elle échoue, B, le jacobien ou l'assemblage sont faux."""
    mesh = Mesh(1.0, 1.0, 4, 4)
    s = FEMSolver(mesh, LinearElastic(210e9, 0.3))
    K, _, _ = s._assemble(np.zeros(s.ndof), 0.0, s._init_state())

    a, b, c, d = 1e-4, 2e-4, -3e-4, 5e-4
    x, y = s.nodes[:, 0], s.nodes[:, 1]
    U_exact = np.empty(s.ndof)
    U_exact[0::2] = a * x + b * y
    U_exact[1::2] = c * x + d * y

    bnd = np.unique(np.concatenate([mesh.bnds[k] for k in ("left", "right", "top", "bottom")]))
    bdofs = np.sort(np.concatenate([2 * bnd, 2 * bnd + 1]))
    free = np.setdiff1d(np.arange(s.ndof), bdofs)

    U = U_exact.copy()
    rhs = -np.asarray(K[free][:, bdofs] @ U_exact[bdofs]).ravel()
    U[free] = spsolve(K[free][:, free].tocsc(), rhs)

    assert np.allclose(U, U_exact, atol=1e-14)


def test_jacobien_negatif_detecte():
    """Un maillage retourné doit lever, pas produire des résultats silencieux."""
    mesh = Mesh(1.0, 1.0, 2, 2)
    mesh.conn = mesh.conn[:, ::-1]          # inversion du sens de parcours
    with pytest.raises(RuntimeError, match="jacobien"):
        FEMSolver(mesh, LinearElastic(210e9, 0.3))


def test_von_mises_cisaillement_pur():
    """En cisaillement pur σ_xy = τ, l'équivalent de von Mises vaut √3·τ."""
    tau = 1e5
    sig4 = np.array([[0.0, 0.0, tau, 0.0]])
    p, q = CrushableCap._invariants(sig4)
    assert p == pytest.approx(0.0)
    assert q == pytest.approx(np.sqrt(3.0) * tau)


@pytest.mark.parametrize("nx", [16, 24, 32, 40, 48, 64, 57])
def test_resultante_de_charge(nx):
    """Les forces nodales doivent totaliser exactement p × span, quel que soit
    le maillage.

    La version précédente sélectionnait les nœuds situés dans la zone chargée
    puis les pondérait par la règle du trapèze : la résultante valait alors
    p × (abscisse du dernier nœud − abscisse du premier), soit 80 % à 100 % de
    la charge prescrite selon la finesse. Ce test ne passait qu'en 32×14, par
    coïncidence d'alignement — d'où le paramétrage sur plusieurs maillages,
    dont un nombre impair d'éléments où rien ne tombe juste.
    """
    s = FEMSolver(Mesh(0.8, 0.35, nx, 14), LinearElastic(E_BETON, NU_BETON))
    p_lin, span = 300e3, 0.25
    F = s._load_vector(p_lin, span)
    assert F.sum() == pytest.approx(-p_lin * span, rel=1e-12)
    assert np.all(F[0::2] == 0.0)           # charge purement verticale


def test_charge_repartie_symetriquement():
    """La semelle étant centrée, les forces nodales doivent l'être aussi."""
    mesh = Mesh(0.8, 0.35, 37, 14)          # maillage impair : bords dans un élément
    s = FEMSolver(mesh, LinearElastic(E_BETON, NU_BETON))
    F = s._load_vector(300e3, 0.25)
    top = np.asarray(mesh.bnds["top"])
    x = mesh.nodes[top, 0]
    ordre = np.argsort(x)
    fy = F[2 * top[ordre] + 1]
    assert np.allclose(fy, fy[::-1], atol=1e-9 * abs(fy).max())


# ─────────────────────────────────────────────────────────────────────────────
# 2. Validation contre solutions analytiques
# ─────────────────────────────────────────────────────────────────────────────
def test_dilatation_libre():
    """Couche encastrée en base, bords à glissement, sommet libre, chauffée de ΔT.

    Seule ε_yy est libre et σ_yy = 0, d'où ε_yy = 3K·αΔT/(K + 4G/3) et un
    soulèvement uniforme u_y = ε_yy · L. C'est le test qui attrape le facteur
    (1+ν) de la réduction en déformation plane.
    """
    Ly, dT = 0.35, 40.0
    mesh = Mesh(0.8, Ly, 32, 14)
    s = FEMSolver(mesh, LinearElastic(E_BETON, NU_BETON),
                  delta_T=dT, alpha=ALPHA_BETON)
    U, *_ = s.run(nsteps=3, p_total=0.0, sleeper_span=0.25)
    _, uy = s.top_profile(U)

    K, G = moduli(E_BETON, NU_BETON)
    attendu = 3 * K * ALPHA_BETON * dT / (K + 4 * G / 3) * Ly * 1e6   # µm

    assert uy.max() - uy.min() < 1e-6, "le soulèvement doit être uniforme"
    assert uy.mean() == pytest.approx(attendu, rel=1e-10)


def test_signe_thermique():
    """Un échauffement soulève la surface, un refroidissement l'abaisse.

    C'est le bug qui rendait la conclusion du rapport initial fausse : la
    contrainte thermique entrait dans le résidu avec le signe opposé.
    """
    mesh = Mesh(0.8, 0.35, 16, 8)
    haut = FEMSolver(mesh, LinearElastic(E_BETON, NU_BETON),
                     delta_T=+40.0, alpha=ALPHA_BETON)
    bas = FEMSolver(mesh, LinearElastic(E_BETON, NU_BETON),
                    delta_T=-40.0, alpha=ALPHA_BETON)
    _, uy_chaud = haut.top_profile(haut.run(nsteps=2, p_total=0.0)[0])
    _, uy_froid = bas.top_profile(bas.run(nsteps=2, p_total=0.0)[0])

    assert uy_chaud.mean() > 0
    assert uy_froid.mean() < 0
    assert uy_chaud.mean() == pytest.approx(-uy_froid.mean(), rel=1e-12)


def test_elasticite_lineaire_est_lineaire():
    """Doubler la charge doit doubler le déplacement (sans plasticité ni ΔT)."""
    mesh = Mesh(0.8, 0.35, 16, 8)
    res = []
    for p in (150e3, 300e3):
        s = FEMSolver(mesh, LinearElastic(E_BETON, NU_BETON))
        U, *_ = s.run(nsteps=2, p_total=p, sleeper_span=0.25)
        res.append(s.top_profile(U)[1])
    assert np.allclose(res[1], 2.0 * res[0], rtol=1e-10)


def test_superposition_thermique_et_mecanique():
    """En élasticité linéaire, charge et thermique se superposent exactement."""
    mesh = Mesh(0.8, 0.35, 16, 8)

    def profil(p, dT):
        s = FEMSolver(mesh, LinearElastic(E_BETON, NU_BETON),
                      delta_T=dT, alpha=ALPHA_BETON)
        return s.top_profile(s.run(nsteps=2, p_total=p, sleeper_span=0.25)[0])[1]

    assert np.allclose(profil(300e3, 40.0), profil(300e3, 0.0) + profil(0.0, 40.0),
                       rtol=1e-10)


def test_symetrie_du_profil():
    """Géométrie, matériau et chargement symétriques ⇒ profil symétrique."""
    mesh = Mesh(0.8, 0.35, 32, 14)
    s = FEMSolver(mesh, CrushableCap(**BALLAST_PARAMS), is_plastic=True,
                  delta_T=40.0, alpha=1.2e-5)
    U, *_ = s.run(nsteps=4, p_total=300e3, sleeper_span=0.25)
    _, uy = s.top_profile(U)
    # tolérance rapportée à l'amplitude du profil, pas à chaque valeur (certaines
    # passent par zéro, où une tolérance relative n'a pas de sens)
    amplitude = uy.max() - uy.min()
    assert np.abs(uy - uy[::-1]).max() <= 1e-9 * amplitude


def test_ballast_s_enfonce_beaucoup_plus_que_le_beton():
    """Le béton est 200 fois plus rigide et reste élastique : sa cuvette doit
    être inférieure de plusieurs ordres de grandeur.

    Le rapport initial concluait l'inverse (rapport 0,89) parce que les deux
    courbes étaient dominées par un décalage thermique uniforme.
    """
    from fem_ballast_beton import compare
    r = compare(nx=24, ny=10, nsteps=4)
    assert r["bowl_ballast"] > 100 * r["bowl_beton"]
    assert r["uy_beton"].min() > 0, "le béton, très rigide, reste soulevé par ΔT"
    assert r["uy_ballast"].min() < 0, "le ballast s'enfonce sous la traverse"


# ─────────────────────────────────────────────────────────────────────────────
# 3. Robustesse numérique
# ─────────────────────────────────────────────────────────────────────────────
def test_independance_au_nombre_d_increments():
    """La réponse plastique ne doit pas dépendre du découpage du chargement.

    Deux défauts distincts se manifestaient ici. La version initiale divisait la
    charge par nsteps sans jamais l'accumuler : le résultat était inversement
    proportionnel au nombre d'incréments. Et reconstruire σ_zz par ν(σ_xx+σ_yy)
    après un retour plastique faisait dériver la réponse de 17 % entre 4 et 24
    incréments. Aucune valeur de référence codée en dur ici : les résultats sont
    comparés entre eux, ce qui est précisément la propriété recherchée.
    """
    mesh = Mesh(0.8, 0.35, 24, 10)
    tassements = []
    for nsteps in (4, 8, 16):
        s = FEMSolver(mesh, CrushableCap(**BALLAST_PARAMS), is_plastic=True,
                      delta_T=40.0, alpha=1.2e-5)
        U, *_ = s.run(nsteps=nsteps, p_total=300e3, sleeper_span=0.25)
        tassements.append(s.top_profile(U)[1].min())

    ecart = (max(tassements) - min(tassements)) / abs(np.mean(tassements))
    assert ecart < 0.02, f"dispersion {ecart:.1%} sur {tassements}"


def test_convergence_en_maillage():
    """Le raffinement doit converger : l'écart entre niveaux successifs décroît."""
    vals = []
    for nx, ny in ((16, 7), (32, 14), (64, 28)):
        mesh = Mesh(0.8, 0.35, nx, ny)
        s = FEMSolver(mesh, LinearElastic(E_BETON, NU_BETON),
                      delta_T=40.0, alpha=ALPHA_BETON)
        U, *_ = s.run(nsteps=2, p_total=300e3, sleeper_span=0.25)
        _, uy = s.top_profile(U)
        vals.append(uy[0] - uy.min())

    e1, e2 = abs(vals[1] - vals[0]), abs(vals[2] - vals[1])
    assert e2 < e1, f"pas de convergence : écarts {e1:.3e} puis {e2:.3e}"


def test_convergence_en_maillage_du_cas_plastique():
    """Le cas plastique doit converger de façon monotone.

    Il ne le faisait pas tant que la résultante appliquée dépendait du
    maillage : l'étude comparait alors des calculs sous charges différentes.
    """
    from fem_ballast_beton import compare
    cuvettes = [compare(nx=nx, ny=ny)["bowl_ballast"]
                for nx, ny in ((16, 7), (24, 10), (32, 14), (48, 21))]
    ecarts = np.abs(np.diff(cuvettes))
    assert np.all(np.diff(cuvettes) > 0), f"convergence non monotone : {cuvettes}"
    assert ecarts[-1] < ecarts[0] / 3, f"convergence trop lente : {ecarts}"


def test_newton_converge_en_peu_d_iterations():
    """Avec la tangente consistante, chaque incrément doit tenir en une dizaine
    d'itérations. Une tangente approchée en demandait plus de cent."""
    mesh = Mesh(0.8, 0.35, 24, 10)
    s = FEMSolver(mesh, CrushableCap(**BALLAST_PARAMS), is_plastic=True,
                  delta_T=40.0, alpha=1.2e-5)
    s.run(nsteps=6, p_total=300e3, sleeper_span=0.25, max_iter=60)
    assert max(s.iterations) <= 15, f"itérations par incrément : {s.iterations}"


# ─────────────────────────────────────────────────────────────────────────────
# 4. Loi de comportement
# ─────────────────────────────────────────────────────────────────────────────
def test_cap_reste_elastique_sous_faible_contrainte():
    """Sous le seuil, aucune déformation plastique ne doit apparaître."""
    mat = CrushableCap(**BALLAST_PARAMS)
    state = mat.init_state(1)
    eps = np.array([[1e-7, -1e-7, 0.0]])
    _, _, new = mat.stress_update(eps, np.zeros((1, 3)), state)
    assert new["evp"][0] == pytest.approx(0.0, abs=1e-18)
    assert new["pc"][0] == pytest.approx(BALLAST_PARAMS["pc0"])


def test_cap_se_consolide_en_compression():
    """En compaction, pc doit croître (consolidation) et B augmenter (écrasement)."""
    mat = CrushableCap(**BALLAST_PARAMS)
    state = mat.init_state(1)
    for _ in range(5):
        eps = state["eps"][:, :3] - np.array([[2e-3, 2e-3, 0.0]])
        _, _, state = mat.stress_update(eps, np.zeros((1, 3)), state)

    assert state["pc"][0] > BALLAST_PARAMS["pc0"]
    assert state["evp"][0] > 0.0
    assert 0.0 < state["B"][0] < 1.0


def test_cap_suit_l_elasticite_hors_plasticite():
    """Hors surface de charge, la loi doit coïncider avec Hooke élastique."""
    mat = CrushableCap(**BALLAST_PARAMS)
    state = mat.init_state(1)
    deps = np.array([[3e-8, -1e-8, 2e-8]])
    sig, _, _ = mat.stress_update(deps, np.zeros((1, 3)), state)

    K, G = moduli(mat.E0, mat.nu0)
    lam = K - 2 * G / 3
    tr = deps[0, 0] + deps[0, 1]
    attendu = np.array([lam * tr + 2 * G * deps[0, 0],
                        lam * tr + 2 * G * deps[0, 1],
                        G * deps[0, 2]])
    assert np.allclose(sig[0], attendu, rtol=1e-9)
