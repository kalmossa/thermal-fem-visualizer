"""
FEM 2D — Ballast (crushable-cap) vs Béton (élastique linéaire), avec effets thermiques.

Formulation
-----------
Éléments quadrangulaires Q4 bilinéaires, quadrature de Gauss 2x2, déformation plane.

Le problème est résolu par Newton-Raphson incrémental sur un chargement
*proportionnel* : un facteur de charge lambda croît de 0 à 1 et pilote
simultanément la pression mécanique et la variation de température. À chaque
incrément on résout l'équilibre

    f_int(U) = f_ext        avec   f_int(U) = ∫ Bᵀ σ(U) dΩ
    K_t ΔU   = f_ext − f_int(U)     puis     U ← U + ΔU

La contrainte thermique entre par la loi de comportement, σ = D:(ε − ε_th) avec
ε_th = α ΔT [1, 1, 0]ᵀ. Elle apparaît donc automatiquement dans le résidu, avec
le bon signe : un échauffement (ΔT > 0) soulève la surface libre.

Convention de signe : traction positive. Une compression donne p < 0.

L'assemblage est vectorisé : les matrices B et les jacobiens sont précalculés une
fois pour toutes, les matrices élémentaires sont obtenues par `einsum` et la
matrice globale par un unique `coo_matrix` (qui somme les doublons).
"""
from dataclasses import dataclass, field

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve

M_VOL = np.array([1.0, 1.0, 0.0])  # opérateur "trace" en notation de Voigt 2D


class SolverDivergence(RuntimeError):
    """Le calcul est sorti du domaine de validité du modèle.

    Type distinct plutôt que RuntimeError générique : l'appelant doit pouvoir
    répondre « ce jeu de paramètres dépasse le modèle » sans attraper au passage
    les vraies erreurs de programmation. Une divergence non détectée se propage
    en NaN, que Flask sérialise tel quel — et `NaN` n'est pas du JSON valide,
    si bien que le client reçoit un 200 accompagné d'un corps illisible.
    """

    def __init__(self, message, etape=None):
        super().__init__(message)
        self.etape = etape


# ─────────────────────────────────────────────────────────────────────────────
# Maillage
# ─────────────────────────────────────────────────────────────────────────────
@dataclass
class Mesh:
    """Maillage rectangulaire structuré d'éléments Q4."""

    lx: float
    ly: float
    nx: int
    ny: int
    nodes: np.ndarray = field(default=None, init=False)
    conn: np.ndarray = field(default=None, init=False)
    bnds: dict = field(default=None, init=False)

    def __post_init__(self):
        self.rect_mesh(self.lx, self.ly, self.nx, self.ny)

    def rect_mesh(self, lx, ly, nx, ny):
        x = np.linspace(0.0, lx, nx + 1)
        y = np.linspace(0.0, ly, ny + 1)
        xv, yv = np.meshgrid(x, y, indexing="xy")
        self.nodes = np.column_stack([xv.ravel(), yv.ravel()])

        conn = []
        for j in range(ny):
            for i in range(nx):
                n0 = j * (nx + 1) + i
                n1 = n0 + 1
                n3 = n0 + (nx + 1)
                n2 = n3 + 1
                conn.append([n0, n1, n2, n3])
        self.conn = np.array(conn, dtype=int)

        tol = 1e-12
        self.bnds = {
            "left":   np.where(np.isclose(self.nodes[:, 0], 0.0, atol=tol))[0],
            "right":  np.where(np.isclose(self.nodes[:, 0], lx, atol=tol))[0],
            "bottom": np.where(np.isclose(self.nodes[:, 1], 0.0, atol=tol))[0],
            "top":    np.where(np.isclose(self.nodes[:, 1], ly, atol=tol))[0],
        }

    def generate(self):
        return self.nodes, self.conn, self.bnds


def gauss_points():
    """Points de Gauss 2x2 sur l'élément de référence [-1, 1]²."""
    g = 1.0 / np.sqrt(3.0)
    return [(-g, -g, 1.0), (g, -g, 1.0), (g, g, 1.0), (-g, g, 1.0)]


def shape_Q4(xi, eta):
    """Fonctions de forme bilinéaires et leurs dérivées dans l'élément de référence."""
    N = 0.25 * np.array([
        (1 - xi) * (1 - eta),
        (1 + xi) * (1 - eta),
        (1 + xi) * (1 + eta),
        (1 - xi) * (1 + eta),
    ])
    dN_dxi = 0.25 * np.array([
        [-(1 - eta), -(1 - xi)],
        [(1 - eta),  -(1 + xi)],
        [(1 + eta),   (1 + xi)],
        [-(1 + eta),  (1 - xi)],
    ])
    return N, dN_dxi


def B_matrix(dN_dx):
    """Matrice déformation-déplacement (3x8) pour un Q4, à partir de dN/dx (4x2)."""
    B = np.zeros((3, 8))
    for a in range(4):
        B[0, 2 * a] = dN_dx[a, 0]
        B[1, 2 * a + 1] = dN_dx[a, 1]
        B[2, 2 * a] = dN_dx[a, 1]
        B[2, 2 * a + 1] = dN_dx[a, 0]
    return B


# ─────────────────────────────────────────────────────────────────────────────
# Lois de comportement
# ─────────────────────────────────────────────────────────────────────────────
class LinearElastic:
    """Hooke isotrope en déformation plane."""

    def __init__(self, E, nu, plane_strain=True):
        self.E = E
        self.nu = nu
        self.plane_strain = plane_strain

    def D(self):
        E, nu = self.E, self.nu
        c = E / ((1 + nu) * (1 - 2 * nu))
        return c * np.array([
            [1 - nu, nu, 0.0],
            [nu, 1 - nu, 0.0],
            [0.0, 0.0, (1 - 2 * nu) / 2],
        ])

    def init_state(self, ngp):
        return None

    def stress_update(self, eps, eps_th, state):
        """σ = D:(ε − (1+ν)·ε_th). Renvoie (σ, D, état inchangé).

        Le facteur (1+ν) n'est pas un ajustement : il tombe de la réduction en
        déformation plane. La dilatation libre est isotrope (α ΔT sur les trois
        directions), mais la condition ε_zz = 0 bloque la composante hors-plan et
        réinjecte σ_zz dans le plan. En repassant par la loi 3D complète on
        trouve σ = D_plan:(ε − (1+ν) α ΔT [1,1,0]ᵀ).

        Contrôle : une couche d'épaisseur L, encastrée en base, bords à
        glissement et sommet libre, chauffée de ΔT, se soulève de
        u_y = 3K α ΔT L / (K + 4G/3) — c'est le test `test_dilatation_libre`.
        """
        sig = (eps - (1.0 + self.nu) * eps_th) @ self.D().T
        return sig, self.D(), state


class CrushableCap:
    """Modèle élasto-plastique à cap pour matériau granulaire écrasable.

    Surface de charge dans le plan (p̄, q), famille Cam-Clay / cap :

        f(p̄, q, pc) = (p̄ − pc)² / X²  +  (q / (M·pc))²  −  1

    C'est une ellipse centrée en p̄ = pc, de demi-axes X en pression et M·pc en
    déviateur.

    ⚠ p̄ = −p est la pression moyenne en convention *compression positive*.
    Le reste du code travaille en traction positive ; la conversion se fait ici,
    à l'entrée et à la sortie de `stress_update`. Cette distinction n'est pas
    cosmétique : un cap doit s'ouvrir du côté comprimé, sinon l'écrouissage
    (pc croissant) déplace le plafond de compression dans le mauvais sens et
    *réduit* la capacité portante au lieu de l'augmenter.

    Écoulement associé            : dε^p = λ ∂f/∂σ
    Écrouissage (consolidation)   : dpc = −H_cap · dε_v^p, **monotone** — la
                                    compaction ouvre le cap et rien ne le
                                    referme : l'écrasement des grains ne se
                                    défait pas à la décharge
    Indice d'écrasement des grains: dB  = kB · |dε_v^p|, qui ramollit le module
                                    volumique via K = K0 / (1 + m_soft·B).
                                    L'écrasement des grains est irréversible,
                                    d'où la valeur absolue ici.

    Le retour sur la surface est résolu par bisection sur le multiplicateur
    plastique λ — moins élégant qu'un Newton local, mais inconditionnellement
    robuste et vectorisable sur tous les points de Gauss à la fois.
    """

    def __init__(self, E0, nu0, M, pc0, X, H_cap, m_soft, kB):
        self.E0 = E0
        self.nu0 = nu0
        self.M = M
        self.pc0 = pc0
        self.X = X
        self.H_cap = H_cap
        self.m_soft = m_soft
        self.kB = kB

    # ── état interne ─────────────────────────────────────────────────────────
    # Les champs `eps` et `sig` ont 4 composantes : [xx, yy, xy, zz]. Suivre σ_zz
    # explicitement est indispensable dès qu'il y a plasticité — la relation
    # élastique σ_zz = ν(σ_xx + σ_yy) cesse d'être vraie après un retour sur la
    # surface de charge, et la reconstruire à chaque pas perd la part plastique
    # hors-plan (l'erreur s'accumule alors avec le nombre d'incréments).
    def init_state(self, ngp):
        return {
            "eps": np.zeros((ngp, 4)),   # déformation mécanique convergée
            "sig": np.zeros((ngp, 4)),   # contrainte convergée
            "pc":  np.full(ngp, self.pc0),
            "B":   np.zeros(ngp),        # indice d'écrasement ∈ [0, 1)
            "evp": np.zeros(ngp),        # déformation volumique plastique cumulée
        }

    @staticmethod
    def _hooke4(deps4, K, G):
        """Incrément de contrainte [xx, yy, xy, zz] par la loi de Hooke 3D isotrope."""
        lam = K - 2.0 * G / 3.0
        tr = deps4[:, 0] + deps4[:, 1] + deps4[:, 3]
        return np.column_stack([
            lam * tr + 2.0 * G * deps4[:, 0],
            lam * tr + 2.0 * G * deps4[:, 1],
            G * deps4[:, 2],                     # γ_xy est la distorsion d'ingénieur
            lam * tr + 2.0 * G * deps4[:, 3],
        ])

    # ── élasticité courante (dégradée par l'écrasement) ──────────────────────
    def _moduli(self, Bidx):
        K0 = self.E0 / (3.0 * (1.0 - 2.0 * self.nu0))
        G0 = self.E0 / (2.0 * (1.0 + self.nu0))
        K = K0 / (1.0 + self.m_soft * Bidx)
        E = 9.0 * K * G0 / (3.0 * K + G0)
        nu = (3.0 * K - 2.0 * G0) / (2.0 * (3.0 * K + G0))
        return K, np.full_like(Bidx, G0), E, nu

    @staticmethod
    def _D_from(E, nu):
        """Matrice de Hooke en déformation plane, une par point de Gauss."""
        c = E / ((1.0 + nu) * (1.0 - 2.0 * nu))
        D = np.zeros(E.shape + (3, 3))
        D[..., 0, 0] = c * (1 - nu)
        D[..., 0, 1] = c * nu
        D[..., 1, 0] = c * nu
        D[..., 1, 1] = c * (1 - nu)
        D[..., 2, 2] = c * (1 - 2 * nu) / 2
        return D

    @classmethod
    def _D_from_KG(cls, K, G):
        """Hooke en déformation plane à partir de (K, G).

        ν est borné : la tangente algorithmique peut faire chuter K jusqu'à
        produire un ν non physique, ce qui n'affecterait que le conditionnement
        (la solution converge sur le résidu, pas sur la tangente).
        """
        nu = np.clip((3.0 * K - 2.0 * G) / (2.0 * (3.0 * K + G)), -0.9, 0.49)
        E = 9.0 * K * G / (3.0 * K + G)
        return cls._D_from(E, nu)

    @staticmethod
    def _invariants(sig4):
        """Pression moyenne p et déviateur de von Mises q, depuis [xx, yy, xy, zz].

        q = √(3·J₂) avec J₂ = ½ s_ij s_ij, la double contraction comptant deux
        fois le terme de cisaillement hors-diagonale.
        """
        sxx, syy, sxy, szz = sig4[:, 0], sig4[:, 1], sig4[:, 2], sig4[:, 3]
        p = (sxx + syy + szz) / 3.0
        s = np.stack([sxx - p, syy - p, szz - p], axis=1)
        J2 = 0.5 * (s ** 2).sum(axis=1) + sxy ** 2
        return p, np.sqrt(np.maximum(0.0, 3.0 * J2))

    # ── retour radial ────────────────────────────────────────────────────────
    def _cap_state(self, lam, pb_tr, q_tr, pc_n, K, G):
        """(f, p̄, q, pc, dε_v^p) pour un multiplicateur plastique λ donné.

        Tout est exprimé en compression positive (p̄ = −p).

        Schéma implicite : les gradients sont évalués à l'état mis à jour. Les
        trois inconnues (p̄, q, pc) se découplent alors en formes fermées.

        Pour pc : dε_v^p = −c·(p̄_tr − pc) avec c = 2λ/(X²(1+a)), et l'écrouissage
        dpc = −H_cap·dε_v^p donne

            pc = (pc_n + H_cap·c·p̄_tr) / (1 + H_cap·c)

        soit une combinaison convexe de pc_n et p̄_tr. Résoudre cette équation par
        point fixe au lieu de l'inverser diverge dès que H_cap approche le module
        volumique K (le rapport H_cap/K est le facteur de contraction de la
        récurrence) — d'où la forme explicite.

        L'écrouissage est piloté par dε_v^p *signé* et non par sa valeur absolue :
        la compaction referme le cap, la dilatance le rouvre. Prendre |dε_v^p|
        ferait durcir le matériau y compris quand il se décomprime.
        """
        X2 = self.X ** 2
        a = 2.0 * K * lam / X2
        c = 2.0 * lam / (X2 * (1.0 + a))
        hc = self.H_cap * c
        # L'écrouissage est monotone : le cap ne se referme jamais. La
        # compaction des grains est irréversible, une décharge ne la défait pas.
        pc = np.maximum(pc_n, (pc_n + hc * pb_tr) / (1.0 + hc))
        # dε_v^p = λ ∂f/∂p = −c(p̄ − pc) : négatif en compression (compaction)
        devp = -c * (pb_tr - pc)
        pb = (pb_tr + a * pc) / (1.0 + a)
        Mpc2 = (self.M * np.maximum(pc, 1e-12)) ** 2
        q = q_tr / (1.0 + 6.0 * G * lam / Mpc2)
        f = (pb - pc) ** 2 / X2 + q ** 2 / Mpc2 - 1.0
        return f, pb, q, pc, devp

    def _solve_lambda(self, pb_tr, q_tr, pc_n, K, G):
        """Plus petit λ ≥ 0 ramenant l'état sur la surface, par encadrement + bisection."""
        # λ = X²/(2K) est l'échelle naturelle du problème (elle rend a = 2Kλ/X²
        # égal à 1) ; partir de là évite une trentaine de doublements à vide.
        lo = np.zeros_like(pb_tr)
        hi = 1e-6 * self.X ** 2 / (2.0 * K)
        bracketed = np.zeros_like(pb_tr, dtype=bool)
        for _ in range(60):
            f_hi = self._cap_state(hi, pb_tr, q_tr, pc_n, K, G)[0]
            bracketed |= f_hi <= 0.0
            if np.all(bracketed):
                break
            hi = np.where(bracketed, hi, hi * 8.0)
        if not np.all(bracketed):
            raise SolverDivergence(
                "retour sur le cap impossible : l'écrouissage H_cap est trop faible "
                "pour le niveau de contrainte atteint — la surface de charge ne "
                "rattrape jamais l'état d'essai.", etape="retour radial")
        for _ in range(80):
            mid = 0.5 * (lo + hi)
            f_m = self._cap_state(mid, pb_tr, q_tr, pc_n, K, G)[0]
            outside = f_m > 0.0
            lo = np.where(outside, mid, lo)
            hi = np.where(outside, hi, mid)
        return 0.5 * (lo + hi)

    def _return_map(self, eps_m, state, masque=None):
        """Prédicteur élastique puis retour sur le cap, à déformation mécanique donnée.

        Fonction pure : `state` n'est pas modifié. Renvoie (σ à 4 composantes,
        pc, dε_v^p, masque actif).

        `masque` force l'ensemble des points plastifiants au lieu de le
        redéterminer. La tangente par différences finies s'en sert : sans cela
        la perturbation ferait basculer certains points d'un côté à l'autre du
        critère, et la dérivée mesurerait la discontinuité du basculement plutôt
        que la pente de la loi. Newton reçoit alors une tangente aberrante et
        cesse de converger.
        """
        K, G, _, _ = self._moduli(state["B"])
        sig_tr = state["sig"] + self._hooke4(eps_m - state["eps"], K, G)

        pc_n = state["pc"]
        p_tr, q_tr = self._invariants(sig_tr)
        pb_tr = -p_tr                      # passage en compression positive
        X2 = self.X ** 2
        Mpc2_n = (self.M * np.maximum(pc_n, 1e-12)) ** 2
        f_tr = (pb_tr - pc_n) ** 2 / X2 + q_tr ** 2 / Mpc2_n - 1.0

        # Le cap n'est actif que du côté compaction (p̄ > pc). L'ellipse est une
        # surface fermée : en dessous de p̄ = pc − X elle est de nouveau franchie,
        # mais par le bas — c'est-à-dire en décompression. Y déclencher un retour
        # plastique fait décroître pc, qui poursuit la contrainte vers le bas :
        # la décharge diverge. Un modèle de cap complet borne ce côté par une
        # enveloppe de rupture en cisaillement ; celle-ci n'est pas implémentée,
        # le domaine y est donc traité comme élastique. Le modèle ne doit pas
        # être utilisé sous faible confinement.
        # Le cap n'est actif que sur sa branche de compaction (p̄ > pc).
        # L'ellipse est une surface fermée : sous p̄ = pc − X on la franchit de
        # nouveau, mais par le bas, c'est-à-dire en décompression. Y déclencher
        # un retour plastique fait que le matériau résiste à la décharge, et la
        # décharge diverge. Un modèle de cap complet borne ce côté par une
        # enveloppe de rupture en cisaillement ; elle n'est pas implémentée, le
        # domaine y est donc élastique — le modèle ne vaut pas sous faible
        # confinement.
        act = (f_tr > 0.0) if masque is None else masque
        pb, q, pc = pb_tr.copy(), q_tr.copy(), pc_n.copy()
        devp = np.zeros_like(pb_tr)
        if np.any(act):
            lam = self._solve_lambda(pb_tr[act], q_tr[act], pc_n[act], K[act], G[act])
            _, pb_a, q_a, pc_a, devp_a = self._cap_state(
                lam, pb_tr[act], q_tr[act], pc_n[act], K[act], G[act]
            )
            pb[act], q[act], pc[act], devp[act] = pb_a, q_a, pc_a, devp_a
        p = -pb                            # retour en traction positive

        # reconstruction du tenseur : retour radial sur la partie déviatorique
        scale = np.where(q_tr > 1e-12, q / np.where(q_tr > 1e-12, q_tr, 1.0), 1.0)
        sig = np.column_stack([
            p + scale * (sig_tr[:, 0] - p_tr),
            p + scale * (sig_tr[:, 1] - p_tr),
            scale * sig_tr[:, 2],
            p + scale * (sig_tr[:, 3] - p_tr),
        ])
        return sig, pc, devp, act

    def stress_update(self, eps, eps_th, state):
        """Intégration de la loi sur l'incrément courant.

        `eps` est la déformation *totale*, `eps_th` la déformation thermique
        libre ; l'incrément mécanique se mesure depuis le dernier état convergé.
        La fonction est pure : l'état renvoyé n'est commité par le solveur
        qu'une fois l'incrément convergé.
        """
        # Déformation mécanique à 4 composantes. La dilatation libre est isotrope
        # et agit donc aussi hors-plan ; comme ε_zz = 0 en déformation plane, la
        # composante mécanique hors-plan vaut −α ΔT.
        th = eps_th[0, 0]
        eps_m = np.column_stack([
            eps[:, 0] - th,
            eps[:, 1] - th,
            eps[:, 2],
            np.full(eps.shape[0], -th),
        ])
        sig, pc, devp, act = self._return_map(eps_m, state)

        # Tangente consistante, obtenue par différences finies : on perturbe
        # chacune des trois composantes de déformation dans le plan et on mesure
        # la réponse du retour radial complet, écrouissage compris.
        #
        # Une tangente approchée (élastique, ou « à λ gelé » : K/(1+2Kλ/X²) et
        # G/(1+6Gλ/(M·pc)²)) fait converger Newton linéairement — de l'ordre de
        # cent itérations par incrément sur ce cas. Les trois évaluations
        # supplémentaires du retour radial sont largement rentabilisées.
        h = 1e-7
        D_tan = np.empty((eps.shape[0], 3, 3))
        for j in range(3):
            eps_p = eps_m.copy()
            eps_p[:, j] += h
            sig_p, _, _, _ = self._return_map(eps_p, state, masque=act)
            D_tan[:, :, j] = (sig_p[:, :3] - sig[:, :3]) / h

        new_state = {
            "eps": eps_m,
            "sig": sig,
            "pc": pc,
            "B": np.minimum(0.999, state["B"] + self.kB * np.abs(devp)),
            "evp": state["evp"] + np.abs(devp),
        }
        # Seules les composantes dans le plan alimentent les forces internes.
        return sig[:, :3], D_tan, new_state


# ─────────────────────────────────────────────────────────────────────────────
# Solveur
# ─────────────────────────────────────────────────────────────────────────────
class Layer:
    """Une couche horizontale du profil : loi, dilatation, épaisseur.

    Une plateforme ferroviaire réelle est stratifiée — ballast sur sous-couche
    sur sol support — et c'est l'épaisseur de chaque couche qui constitue le
    levier de conception. Le domaine homogène ne permet de comparer que des
    matériaux pris isolément.
    """

    def __init__(self, law, epaisseur, alpha=0.0, is_plastic=False, nom=""):
        if epaisseur <= 0:
            raise ValueError("l'épaisseur d'une couche doit être strictement positive")
        self.law = law
        self.epaisseur = float(epaisseur)
        self.alpha = float(alpha)
        self.is_plastic = bool(is_plastic)
        self.nom = nom or law.__class__.__name__
        self.elements = None          # rempli par le solveur
        self.gp = None                # indices des points de Gauss
        self.y_min = self.y_max = None

    def __repr__(self):
        return f"Layer({self.nom!r}, e={self.epaisseur:.3f} m, α={self.alpha:.1e})"


class FEMSolver:
    """Solveur Newton-Raphson incrémental, assemblage vectorisé.

    Accepte soit une loi de comportement unique (domaine homogène), soit une
    liste de `Layer` décrivant le profil **du haut vers le bas**, à la manière
    dont un projeteur décrit une structure d'assise.
    """

    def __init__(self, mesh, constitutive=None, is_plastic=False, delta_T=0.0,
                 alpha=0.0, layers=None):
        if constitutive is None and not layers:
            raise ValueError("il faut une loi de comportement ou une liste de couches")
        self.mesh = mesh
        self.constitutive = constitutive
        self.is_plastic = is_plastic
        self.delta_T = delta_T
        self.alpha = alpha
        self.nodes, self.conn, self.bnds = mesh.generate()
        self.nelem = self.conn.shape[0]
        self.ngp = self.nelem * 4
        self.ndof = 2 * self.nodes.shape[0]
        self._precompute()
        self._setup_bcs()
        self._setup_layers(layers)
        self.iterations = []

    # ── précalculs géométriques ──────────────────────────────────────────────
    def _precompute(self):
        """Matrices B et poids d'intégration, une fois pour toutes."""
        Xe = self.nodes[self.conn]                      # (nelem, 4, 2)
        self.B = np.zeros((self.nelem, 4, 3, 8))
        self.dV = np.zeros((self.nelem, 4))

        for g, (xi, eta, w) in enumerate(gauss_points()):
            _, dN_dxi = shape_Q4(xi, eta)               # (4, 2)
            # J[k, i] = ∂x_i/∂ξ_k
            J = np.einsum("ak,eai->eki", dN_dxi, Xe)
            detJ = np.linalg.det(J)
            if np.any(detJ <= 0):
                bad = int(np.argmin(detJ))
                raise RuntimeError(f"jacobien non positif sur l'élément {bad} : detJ={detJ[bad]:.3e}")
            Jinv = np.linalg.inv(J)
            # ∂N_a/∂x_i = Σ_k Jinv[i, k] ∂N_a/∂ξ_k  (noter la transposition)
            dN_dx = np.einsum("ak,eik->eai", dN_dxi, Jinv)

            self.B[:, g, 0, 0::2] = dN_dx[:, :, 0]
            self.B[:, g, 1, 1::2] = dN_dx[:, :, 1]
            self.B[:, g, 2, 0::2] = dN_dx[:, :, 1]
            self.B[:, g, 2, 1::2] = dN_dx[:, :, 0]
            self.dV[:, g] = detJ * w

        # table des degrés de liberté et indices d'assemblage
        self.edofs = np.empty((self.nelem, 8), dtype=int)
        self.edofs[:, 0::2] = 2 * self.conn
        self.edofs[:, 1::2] = 2 * self.conn + 1
        self._rows = np.repeat(self.edofs, 8, axis=1).ravel()
        self._cols = np.tile(self.edofs, (1, 8)).ravel()

    def _setup_layers(self, layers):
        """Affecte chaque élément à sa couche, d'après l'ordonnée de son centre.

        Les couches sont données du haut vers le bas ; leurs épaisseurs doivent
        totaliser la hauteur du domaine. Un élément à cheval sur une interface
        serait ambigu : on exige donc que les interfaces tombent sur des lignes
        du maillage, ce qui est de toute façon la seule discrétisation
        acceptable pour un contraste de rigidité.
        """
        if not layers:
            unique = Layer(self.constitutive, self.mesh.ly, self.alpha,
                           self.is_plastic, "homogène")
            unique.elements = np.arange(self.nelem)
            unique.y_min, unique.y_max = 0.0, self.mesh.ly
            self.layers = [unique]
            self.homogene = True
        else:
            total = sum(c.epaisseur for c in layers)
            if abs(total - self.mesh.ly) > 1e-9:
                raise ValueError(
                    f"les épaisseurs totalisent {total:.4f} m pour un domaine de "
                    f"{self.mesh.ly:.4f} m")

            yc = self.nodes[self.conn, 1].mean(axis=1)      # ordonnée du centre
            haut = self.mesh.ly
            for couche in layers:
                couche.y_max = haut
                couche.y_min = haut - couche.epaisseur
                couche.elements = np.where((yc > couche.y_min) & (yc <= couche.y_max + 1e-12))[0]
                if couche.elements.size == 0:
                    raise ValueError(
                        f"la couche {couche.nom!r} ({couche.epaisseur:.3f} m) ne "
                        f"contient aucun élément : raffine le maillage vertical")
                haut = couche.y_min
            self.layers = list(layers)
            self.homogene = False

        vus = np.concatenate([c.elements for c in self.layers])
        if vus.size != self.nelem or np.unique(vus).size != self.nelem:
            raise ValueError("le découpage en couches ne partitionne pas le maillage")

        for couche in self.layers:
            couche.gp = (4 * couche.elements[:, None] + np.arange(4)).ravel()

    def _setup_bcs(self):
        """Encastrement en base, appuis à glissement (u_x = 0) sur les bords latéraux."""
        fixed = set()
        for nid in self.bnds["bottom"]:
            fixed.add(2 * nid)
            fixed.add(2 * nid + 1)
        for nid in np.concatenate([self.bnds["left"], self.bnds["right"]]):
            fixed.add(2 * nid)
        self.fixed = np.array(sorted(fixed), dtype=int)
        mask = np.ones(self.ndof, dtype=bool)
        mask[self.fixed] = False
        self.free = np.where(mask)[0]

    # ── chargement ───────────────────────────────────────────────────────────
    def _load_vector(self, p_total, sleeper_span):
        """Forces nodales cohérentes pour une charge linéique uniforme sur [xL, xR].

        Les forces s'obtiennent en intégrant les fonctions de forme linéaires du
        bord sur l'**intersection** de chaque segment de surface avec la zone
        chargée. La résultante vaut donc exactement p_total × sleeper_span, quel
        que soit le maillage, y compris quand les bords de la semelle tombent à
        l'intérieur d'un élément.

        Sélectionner les nœuds dont l'abscisse est dans [xL, xR] puis les
        pondérer par la règle du trapèze, comme le faisait la version
        précédente, revient à appliquer p × (x du dernier nœud − x du premier),
        c'est-à-dire une charge qui dépend de la discrétisation : entre 80 % et
        100 % de la valeur prescrite selon la finesse du maillage. Toute étude
        de convergence bâtie là-dessus compare des calculs sous charges
        différentes.
        """
        F = np.zeros(self.ndof)
        top = np.asarray(self.bnds["top"])
        x_top = self.nodes[top, 0]
        order = np.argsort(x_top)
        top, x_top = top[order], x_top[order]

        xmid = 0.5 * self.mesh.lx
        xL, xR = xmid - 0.5 * sleeper_span, xmid + 0.5 * sleeper_span
        charges = []

        for i in range(top.size - 1):
            x1, x2 = x_top[i], x_top[i + 1]
            a, b = max(x1, xL), min(x2, xR)
            if b <= a:
                continue
            h = x2 - x1
            # ∫ₐᵇ (x₂−x)/h dx  et  ∫ₐᵇ (x−x₁)/h dx
            f1 = ((x2 - a) ** 2 - (x2 - b) ** 2) / (2.0 * h)
            f2 = ((b - x1) ** 2 - (a - x1) ** 2) / (2.0 * h)
            F[2 * top[i] + 1] -= p_total * f1
            F[2 * top[i + 1] + 1] -= p_total * f2
            charges.extend((top[i], top[i + 1]))

        self.loaded_nodes = np.unique(charges) if charges else np.array([], dtype=int)
        return F

    # ── assemblage ───────────────────────────────────────────────────────────
    def _assemble(self, U, lam, state):
        """Renvoie (K_tangent, f_interne, nouvel état) pour un champ U donné.

        `lam` est le facteur de charge : chaque couche a sa propre dilatation,
        donc sa propre déformation thermique λ·α·ΔT.
        """
        Ue = U[self.edofs]                                     # (nelem, 8)
        eps = np.einsum("egji,ei->egj", self.B, Ue)            # (nelem, 4, 3)

        Ke = np.empty((self.nelem, 8, 8))
        fe = np.empty((self.nelem, 8))
        nouveaux = []

        for i, couche in enumerate(self.layers):
            idx = couche.elements
            eps_c = eps[idx].reshape(-1, 3)
            eps_th = lam * couche.alpha * self.delta_T * M_VOL

            sig_c, D, etat = couche.law.stress_update(eps_c, eps_th[None, :], state[i])
            nouveaux.append(etat)

            B_c, dV_c = self.B[idx], self.dV[idx]
            sig = sig_c.reshape(idx.size, 4, 3)

            if D.ndim == 2:  # loi linéaire : une seule matrice de Hooke
                Ke[idx] = np.einsum("egji,jk,egkl,eg->eil", B_c, D, B_c, dV_c,
                                    optimize=True)
            else:
                Dg = D.reshape(idx.size, 4, 3, 3)
                Ke[idx] = np.einsum("egji,egjk,egkl,eg->eil", B_c, Dg, B_c, dV_c,
                                    optimize=True)
            fe[idx] = np.einsum("egji,egj,eg->ei", B_c, sig, dV_c, optimize=True)

        K = coo_matrix(
            (Ke.ravel(), (self._rows, self._cols)), shape=(self.ndof, self.ndof)
        ).tocsr()
        f_int = np.bincount(self.edofs.ravel(), weights=fe.ravel(), minlength=self.ndof)
        return K, f_int, nouveaux

    # ── résolution ───────────────────────────────────────────────────────────
    def run(self, nsteps=6, p_total=300e3, sleeper_span=0.25, max_iter=30, tol=1e-9):
        """Applique le chargement thermo-mécanique en nsteps incréments.

        Renvoie (U, y_top, x_top, top_nodes) — même signature que la version
        d'origine, pour rester compatible avec l'API Flask.
        """
        U = np.zeros(self.ndof)
        state = self._init_state()
        F_ref = self._load_vector(p_total, sleeper_span)
        self.iterations = []

        for step in range(1, nsteps + 1):
            lam = step / nsteps                 # chargement proportionnel
            U, state, it = self._increment(U, state, lam * F_ref, lam, max_iter, tol)
            self.iterations.append(it)

        self.state = state
        top = np.asarray(self.bnds["top"])
        return U, self.nodes[top, 1], self.nodes[top, 0], top

    def _init_state(self):
        return [c.law.init_state(c.gp.size) for c in self.layers]

    def _increment(self, U, state, F_ext, lam, max_iter=30, tol=1e-9):
        """Un incrément de Newton-Raphson. Renvoie (U, état convergé, itérations).

        La divergence est détectée à chaque itération plutôt qu'en bout de
        chaîne : une fois le champ pollué par un NaN, il se propage dans tout le
        résultat et l'origine devient impossible à situer.
        """
        ref = max(np.linalg.norm(F_ext[self.free]), 1.0)
        # Un déplacement de l'ordre de la hauteur du domaine n'a plus de sens
        # dans une formulation en petites déformations.
        plafond = 10.0 * max(self.mesh.lx, self.mesh.ly)

        for it in range(1, max_iter + 1):
            K, f_int, _ = self._assemble(U, lam, state)
            res = F_ext - f_int
            if not np.all(np.isfinite(res)):
                raise SolverDivergence(
                    "le résidu d'équilibre est devenu non fini", etape="assemblage")
            if np.linalg.norm(res[self.free]) <= tol * ref:
                break
            Kff = K[self.free][:, self.free].tocsc()
            dU = np.zeros(self.ndof)
            dU[self.free] = spsolve(Kff, res[self.free])
            if not np.all(np.isfinite(dU)):
                raise SolverDivergence(
                    "la résolution du système tangent a rendu un incrément non fini — "
                    "la matrice est singulière ou mal conditionnée", etape="système tangent")
            U = U + dU
            if np.max(np.abs(U)) > plafond:
                raise SolverDivergence(
                    f"le déplacement dépasse {plafond:.1f} m, soit dix fois la taille "
                    "du domaine : le calcul a divergé", etape="Newton")
        # l'état n'est commité qu'une fois l'équilibre atteint
        _, _, state = self._assemble(U, lam, state)
        return U, state, it

    def run_cyclic(self, n_cycles=5, p_total=300e3, sleeper_span=0.25,
                   lam_min=0.4, nsteps=3, max_iter=40, tol=1e-9):
        """Chargement cyclique entre une charge résiduelle et la charge de service.

        La voie n'est jamais totalement déchargée : poids propre de la
        superstructure, essieux voisins, précontrainte du bourrage. `lam_min`
        est cette charge résiduelle, en fraction de la charge de service.

        Domaine de validité
        -------------------
        Le modèle **ne supporte pas une décharge complète**, et c'est une
        propriété de la loi, pas du solveur. La surface de charge est une
        ellipse fermée : une fois le cap écroui jusqu'à pc, redescendre
        au-dessous de p̄ = pc − X la franchit de nouveau, par le bas. Le retour
        radial y ramène l'état vers le cap, c'est-à-dire qu'il fait *résister*
        le matériau à la décharge, et le calcul diverge.

        Un modèle de cap complet (type Drucker-Prager cap) borne ce côté par une
        enveloppe de rupture en cisaillement, entre laquelle et le cap le
        domaine est élastique. Cette enveloppe n'est pas implémentée ici :
        en dessous d'environ 40 % de la charge de service, le calcul sort du
        domaine de validité et lève une erreur plutôt que de rendre un résultat.

        Ce que le modèle montre
        -----------------------
        Une **adaptation** (shakedown) après le premier cycle : le cap ne
        s'écrouit que dans le sens de la compaction et ne se dégrade jamais, si
        bien qu'il finit par envelopper l'état atteint. Les cycles suivants sont
        élastiques et le tassement se stabilise. C'est un résultat correct du
        modèle. Reproduire le tassement progressif réel de la voie demanderait
        une loi à dégradation cyclique, où pc décroît avec le nombre de cycles.
        """
        if not 0.0 <= lam_min < 1.0:
            raise ValueError("lam_min doit être dans [0, 1[")
        if lam_min < 0.35:
            raise ValueError(
                f"lam_min = {lam_min:.2f} : sous ~0,35 la décharge sort du cap par "
                "sa branche basse et le calcul diverge. Voir le domaine de "
                "validité dans la documentation de run_cyclic.")

        U = np.zeros(self.ndof)
        state = self._init_state()
        F_ref = self._load_vector(p_total, sleeper_span)
        top = np.asarray(self.bnds["top"])
        i_centre = int(np.argmin(np.abs(self.nodes[top, 0] - 0.5 * self.mesh.lx)))
        dof_centre = 2 * top[i_centre] + 1

        historique = []
        self.iterations = []

        def rampe(depart, cible, U, state):
            for k in range(1, nsteps + 1):
                lam = depart + (cible - depart) * k / nsteps
                U, state, it = self._increment(U, state, lam * F_ref, 1.0, max_iter, tol)
                self.iterations.append(it)
                if not np.all(np.isfinite(U)):
                    raise RuntimeError(
                        "le calcul cyclique a divergé : l'état est sorti du domaine "
                        "de validité du modèle de cap")
            return U, state

        # Montée initiale, charge et température ensemble.
        for k in range(1, nsteps + 1):
            lam = k / nsteps
            U, state, it = self._increment(U, state, lam * F_ref, lam, max_iter, tol)
            self.iterations.append(it)

        for cycle in range(1, n_cycles + 1):
            for phase, cible in (("decharge", lam_min), ("charge", 1.0)):
                depart = 1.0 if phase == "decharge" else lam_min
                U, state = rampe(depart, cible, U, state)
                self.state = state
                historique.append({
                    "cycle": cycle,
                    "phase": phase,
                    "lam": cible,
                    "uy_centre_um": float(U[dof_centre]) * 1e6,
                    "evp_max_pct": self.state_max("evp") * 100,
                    "pc_max_kpa": self.state_max("pc", self.layers[0].law.pc0
                                                 if hasattr(self.layers[0].law, "pc0") else 0.0) / 1e3,
                })

        self.state = state
        self.historique_cyclique = historique
        return U, historique

    # ── post-traitement ──────────────────────────────────────────────────────
    def top_profile(self, U):
        """(x trié en m, u_y en µm) sur la surface supérieure."""
        top = np.asarray(self.bnds["top"])
        x = self.nodes[top, 0]
        order = np.argsort(x)
        uy = np.array([U[2 * n + 1] for n in top[order]]) * 1e6
        return x[order], uy

    def gauss_field(self, key):
        """Champ scalaire aux points de Gauss, moyenné par élément.

        Clés utiles : 'pc', 'B', 'evp'. Les couches élastiques n'ont pas d'état
        interne : elles rendent zéro, pour que le champ reste défini sur tout le
        maillage et directement affichable.
        """
        champ = np.zeros(self.nelem)
        trouve = False
        for couche, etat in zip(self.layers, self.state):
            if isinstance(etat, dict) and key in etat:
                champ[couche.elements] = etat[key].reshape(couche.elements.size, 4).mean(axis=1)
                trouve = True
        return champ if trouve else None

    def state_max(self, key, defaut=0.0):
        """Maximum d'une variable interne sur l'ensemble des couches.

        L'état est une liste — une entrée par couche — depuis l'introduction du
        profil multicouche ; les couches élastiques n'en ont pas.
        """
        valeurs = [np.max(e[key]) for e in self.state
                   if isinstance(e, dict) and key in e]
        return float(max(valeurs)) if valeurs else defaut

    def layer_of_element(self):
        """Indice de couche de chaque élément, pour la cartographie du profil."""
        idx = np.zeros(self.nelem, dtype=int)
        for i, couche in enumerate(self.layers):
            idx[couche.elements] = i
        return idx

    def stress_field(self, U, lam=1.0):
        """Contrainte équivalente de von Mises moyennée par élément, en kPa."""
        eps = np.einsum("egji,ei->egj", self.B, U[self.edofs])
        vm = np.zeros(self.ngp)

        for couche, etat in zip(self.layers, self.state):
            idx = couche.elements
            if isinstance(etat, dict) and "sig" in etat:
                # l'état convergé porte déjà les 4 composantes, σ_zz inclus
                sig4 = etat["sig"]
            else:
                eps_th = lam * couche.alpha * self.delta_T * M_VOL
                sig, _, _ = couche.law.stress_update(
                    eps[idx].reshape(-1, 3), eps_th[None, :], None)
                # ε_zz = [σ_zz − ν(σ_xx+σ_yy)]/E + αΔT = 0
                #   ⇒ σ_zz = ν(σ_xx+σ_yy) − E αΔT
                szz = (couche.law.nu * (sig[:, 0] + sig[:, 1])
                       - couche.law.E * eps_th[0])
                sig4 = np.column_stack([sig, szz])
            _, vm[couche.gp] = CrushableCap._invariants(sig4)

        return vm.reshape(self.nelem, 4).mean(axis=1) / 1e3


# ─────────────────────────────────────────────────────────────────────────────
# Cas de référence : ballast vs béton
# ─────────────────────────────────────────────────────────────────────────────
# H_cap est le module d'écrouissage du cap, dpc/d(−ε_v^p). La valeur de 1e5 Pa
# initialement utilisée n'est pas exploitable : porter la pression de service
# (σ_v ≈ 1,2 MPa sous la traverse) demanderait alors pc ≈ 640 kPa, soit près de
# 600 % de déformation volumique plastique. Le cap ne rattrape jamais l'état
# d'essai et le retour radial échoue.
#
# 3e7 Pa correspond à une densification de l'ordre de 2 % sur une plage de
# contrainte de quelques centaines de kPa, ce qui est l'ordre de grandeur admis
# pour un ballast fraîchement bourré. On obtient un tassement d'environ 2,4 mm
# sous la traverse pour 0,5 % de compaction plastique — cohérent avec les
# tassements de voie mesurés aux premiers cycles de chargement.
BALLAST_PARAMS = dict(E0=150e6, nu0=0.25, M=1.2, pc0=50e3, X=80e3,
                      H_cap=3e7, m_soft=8.0, kB=2.5)
CONCRETE_PARAMS = dict(E=30e9, nu=0.2)


def compare(p_service=300e3, delta_T=40.0, sleeper_span=0.25,
            E_ballast=150e6, alpha_ballast=1.2e-5,
            E_concrete=30e9, alpha_concrete=1.0e-5,
            lx=0.8, ly=0.35, nx=32, ny=14, nsteps=6):
    """Lance les deux simulations et renvoie un dictionnaire de résultats."""
    mesh = Mesh(lx, ly, nx, ny)

    params = dict(BALLAST_PARAMS, E0=E_ballast)
    solver_b = FEMSolver(mesh, CrushableCap(**params), is_plastic=True,
                         delta_T=delta_T, alpha=alpha_ballast)
    U_b, _, _, _ = solver_b.run(nsteps, p_service, sleeper_span)
    x, uy_b = solver_b.top_profile(U_b)

    solver_c = FEMSolver(mesh, LinearElastic(E_concrete, CONCRETE_PARAMS["nu"]),
                         is_plastic=False, delta_T=delta_T, alpha=alpha_concrete)
    U_c, _, _, _ = solver_c.run(nsteps, p_service, sleeper_span)
    _, uy_c = solver_c.top_profile(U_c)

    # La "cuvette" est l'enfoncement sous la traverse mesuré depuis le bord :
    # elle isole la réponse mécanique du soulèvement thermique d'ensemble.
    bowl_b = float(uy_b[0] - uy_b.min())
    bowl_c = float(uy_c[0] - uy_c.min())

    return {
        "x_m": x,
        "uy_ballast": uy_b,
        "uy_beton": uy_c,
        "bowl_ballast": bowl_b,
        "bowl_beton": bowl_c,
        "ratio_bowl": bowl_b / bowl_c if bowl_c else float("inf"),
        "solver_ballast": solver_b,
        "solver_beton": solver_c,
        "U_ballast": U_b,
        "U_beton": U_c,
    }


if __name__ == "__main__":
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    Lx, Ly = 0.8, 0.35
    p_service, sleeper_span, delta_T = 300e3, 0.25, 40.0
    r = compare(p_service=p_service, delta_T=delta_T, sleeper_span=sleeper_span,
                lx=Lx, ly=Ly)
    x, uy_b, uy_c = r["x_m"], r["uy_ballast"], r["uy_beton"]

    header = "x_m,uy_microm"
    np.savetxt("U_TOP_BALLAST.csv", np.column_stack([x, uy_b]),
               delimiter=",", header=header, comments="")
    np.savetxt("U_TOP_BETON.csv", np.column_stack([x, uy_c]),
               delimiter=",", header=header, comments="")

    fig, (ax, ax2) = plt.subplots(2, 1, figsize=(9, 7.5),
                                  gridspec_kw={"height_ratios": [1.1, 1]})

    ax.plot(x * 100, uy_b, "o-", color="#e67e22", ms=4, lw=1.5,
            label="Ballast (élasto-plastique, crushable-cap)")
    ax.plot(x * 100, uy_c, "x--", color="#2980b9", ms=5, lw=1.5,
            label="Béton (linéaire élastique)")
    xmid = 0.5 * Lx
    ax.axvspan((xmid - 0.5 * sleeper_span) * 100, (xmid + 0.5 * sleeper_span) * 100,
               alpha=0.12, color="gray", label="Zone semelle (traverse)")
    ax.set_ylabel("Déplacement vertical $u_y$ (µm)")
    ax.set_title("Comparaison Ballast vs Béton — surface supérieure\n"
                 r"($p_{service}$=300 kN/m, $\Delta T$=+40 °C, EF Q4 déformation plane)")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.4)

    # Second panneau : la réponse purement mécanique, décalage thermique retiré.
    ax2.plot(x * 100, uy_b - uy_b[0], "o-", color="#e67e22", ms=4, lw=1.5,
             label=f"Ballast — cuvette {r['bowl_ballast']:.1f} µm")
    ax2.plot(x * 100, uy_c - uy_c[0], "x--", color="#2980b9", ms=5, lw=1.5,
             label=f"Béton — cuvette {r['bowl_beton']:.2f} µm")
    ax2.axvspan((xmid - 0.5 * sleeper_span) * 100, (xmid + 0.5 * sleeper_span) * 100,
                alpha=0.12, color="gray")
    ax2.set_xlabel("Position x (cm)")
    ax2.set_ylabel("Enfoncement relatif au bord (µm)")
    ax2.set_title("Réponse mécanique seule (soulèvement thermique retiré)", fontsize=10)
    ax2.legend(fontsize=9)
    ax2.grid(True, alpha=0.4)

    plt.tight_layout()
    plt.savefig("w_compare.png", dpi=150)
    plt.savefig("w_compare.pdf")

    print("=" * 62)
    print(f"  Ballast : u_y ∈ [{uy_b.min():+8.2f}, {uy_b.max():+8.2f}] µm"
          f"   cuvette = {r['bowl_ballast']:8.2f} µm")
    print(f"  Béton   : u_y ∈ [{uy_c.min():+8.2f}, {uy_c.max():+8.2f}] µm"
          f"   cuvette = {r['bowl_beton']:8.2f} µm")
    print(f"  Rapport des cuvettes ballast/béton = {r['ratio_bowl']:.0f}")
    print(f"  Rapport des modules E_béton/E_ballast = {30e9 / 150e6:.0f}")
    print("=" * 62)
