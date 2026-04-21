"""
fem_ballast_beton.py — Solveur FEM 2D avec effets thermiques
Compare le comportement mécanique du ballast (élasto-plastique)
et du béton (élastique linéaire) sous charge ferroviaire.

Lancé directement : py fem_ballast_beton.py
→ génère U_TOP_BALLAST.csv, U_TOP_BETON.csv, w_compare.png
"""

from dataclasses import dataclass, field
import numpy as np
from scipy.sparse import lil_matrix, csr_matrix
from scipy.sparse.linalg import spsolve
import matplotlib
matplotlib.use('Agg')  # pas d'affichage graphique (mode serveur)
import matplotlib.pyplot as plt


# ══════════════════════════════════════════════════════════════════════════════
# MAILLAGE
# ══════════════════════════════════════════════════════════════════════════════

@dataclass
class Mesh:
    """Maillage rectangulaire structuré en éléments Q4."""
    lx: float; ly: float; nx: int; ny: int
    nodes: np.ndarray = field(default=None, init=False)
    conn:  np.ndarray = field(default=None, init=False)
    bnds:  dict       = field(default=None, init=False)

    def __post_init__(self):
        self.rect_mesh(self.lx, self.ly, self.nx, self.ny)

    def rect_mesh(self, lx, ly, nx, ny):
        # Coordonnées des nœuds sur grille régulière
        x = np.linspace(0.0, lx, nx+1)
        y = np.linspace(0.0, ly, ny+1)
        xv, yv = np.meshgrid(x, y, indexing='xy')
        self.nodes = np.column_stack([xv.ravel(), yv.ravel()])

        # Table de connectivité : 4 nœuds par élément Q4
        conn = []
        for j in range(ny):
            for i in range(nx):
                n0 = j*(nx+1)+i
                conn.append([n0, n0+1, n0+(nx+1)+1, n0+(nx+1)])
        self.conn = np.array(conn, dtype=int)

        # Identification des bords (tolérance numérique)
        tol = 1e-12
        self.bnds = {
            "left":   np.where(np.isclose(self.nodes[:,0], 0.0, atol=tol))[0],
            "right":  np.where(np.isclose(self.nodes[:,0], lx,  atol=tol))[0],
            "bottom": np.where(np.isclose(self.nodes[:,1], 0.0, atol=tol))[0],
            "top":    np.where(np.isclose(self.nodes[:,1], ly,  atol=tol))[0],
        }

    def generate(self):
        return self.nodes, self.conn, self.bnds


# ══════════════════════════════════════════════════════════════════════════════
# FONCTIONS FEM DE BASE
# ══════════════════════════════════════════════════════════════════════════════

def gauss_points():
    """Points et poids d'intégration de Gauss 2×2."""
    g = 1 / np.sqrt(3.0)
    return [(-g,-g,1.0), (g,-g,1.0), (g,g,1.0), (-g,g,1.0)]


def shape_Q4(xi, eta):
    """Fonctions de forme bilinéaires Q4 et leurs dérivées en (xi, eta)."""
    N = 0.25 * np.array([
        (1-xi)*(1-eta), (1+xi)*(1-eta),
        (1+xi)*(1+eta), (1-xi)*(1+eta)
    ])
    dN = 0.25 * np.array([
        [-(1-eta), -(1-xi)],
        [ (1-eta), -(1+xi)],
        [ (1+eta),  (1+xi)],
        [-(1+eta),  (1-xi)]
    ])
    return N, dN


def B_matrix(dN_dx):
    """Matrice cinématique B (3×8) reliant déformations et déplacements nodaux."""
    B = np.zeros((3, 8))
    for a in range(4):
        B[0, 2*a]   = dN_dx[a, 0]  # ε_xx
        B[1, 2*a+1] = dN_dx[a, 1]  # ε_yy
        B[2, 2*a]   = dN_dx[a, 1]  # γ_xy (termes croisés)
        B[2, 2*a+1] = dN_dx[a, 0]
    return B


# ══════════════════════════════════════════════════════════════════════════════
# LOIS DE COMPORTEMENT
# ══════════════════════════════════════════════════════════════════════════════

class LinearElastic:
    """Élasticité linéaire de Hooke en déformation plane (béton)."""

    def __init__(self, E, nu, plane_strain=True):
        self.E = E
        self.nu = nu

    def D(self):
        """Tenseur d'élasticité (matrice 3×3)."""
        E, nu = self.E, self.nu
        c = E / ((1+nu) * (1-2*nu))
        return c * np.array([
            [1-nu, nu,         0],
            [nu,   1-nu,       0],
            [0,    0,   (1-2*nu)/2]
        ])

    def stress(self, deps, delta_T=0.0, alpha=0.0):
        """Contraintes = D × (déformations mécanique - dilatation thermique)."""
        return self.D() @ (deps - np.array([1,1,0]) * alpha * delta_T)


class CrushableCap:
    """
    Loi élasto-plastique crushable-cap pour le ballast granulaire.
    Surface de charge elliptique en (p, q) avec écrouissage isotrope.
    Paramètres : E0, nu0 (élastiques), M (pente déviatorique), pc0 (pression init),
                 X (taille du cap), H_cap (module écrouissage), m_soft, kB (dégradation).
    """

    def __init__(self, E0, nu0, M, pc0, X, H_cap, m_soft, kB):
        self.E0=E0; self.nu0=nu0; self.M=M; self.pc0=pc0
        self.X=X; self.H_cap=H_cap; self.m_soft=m_soft; self.kB=kB

    def init_gp_state(self):
        """État initial d'un point de Gauss (tout à zéro)."""
        return dict(eps_p=np.zeros(3), evp=0.0, pc=self.pc0, B=0.0)

    def _elastic_matrices(self, Bidx):
        """Modules élastiques dégradés selon l'endommagement B."""
        K0 = self.E0 / (3*(1-2*self.nu0))
        G0 = self.E0 / (2*(1+self.nu0))
        K  = K0 / (1.0 + self.m_soft * Bidx)
        E  = 9*K*G0 / (3*K + G0 + 1e-20)
        nu = (3*K - 2*G0) / (2*(3*K + G0) + 1e-20)
        c  = E / ((1+nu)*(1-2*nu) + 1e-20)
        D  = c * np.array([[1-nu,nu,0],[nu,1-nu,0],[0,0,(1-2*nu)/2]])
        return D, K, G0, nu

    @staticmethod
    def _inv3_from_plane_strain(sig, nu):
        """Pression sphérique p et contrainte déviatorique q depuis tenseur 2D."""
        sxx, syy, sxy = sig
        szz = nu * (sxx + syy)  # déformation plane → contrainte hors-plan
        p = (sxx + syy + szz) / 3.0
        sdev = np.array([sxx, syy, szz]) - p
        J2 = 0.5 * ((sdev[0]-sdev[1])**2 + (sdev[1]-sdev[2])**2 + (sdev[2]-sdev[0])**2)
        return p, np.sqrt(max(0.0, 3*J2))

    def stress_update(self, deps, state, delta_T=0.0, alpha=0.0):
        """
        Retour plastique au cap : si l'état trial dépasse la surface de charge,
        on projette itérativement sur l'ellipse par Newton.
        """
        D_el, K, G, nu = self._elastic_matrices(state['B'])

        # Prédicteur élastique (essai)
        sig_tr = D_el @ (deps - np.array([1,1,0]) * alpha * delta_T)
        p_tr, q_tr = self._inv3_from_plane_strain(sig_tr, nu)
        pc = state['pc']
        X  = max(self.X, 1e-12)

        # Critère de plasticité : f ≤ 0 → élastique, pas de correction
        f_tr = (p_tr-pc)**2/X**2 + (q_tr/(self.M*max(pc,1e-12)))**2 - 1.0
        if f_tr <= 0.0:
            return sig_tr, D_el, state

        # Retour plastique par Newton (max 40 itérations)
        lam = 0.0
        dp_dlam  = -K;          dq_dlam  = -0.6*G
        dpc_dlam =  self.H_cap
        for _ in range(40):
            p      = p_tr  + dp_dlam  * lam
            q      = q_tr  + dq_dlam  * lam
            pc_loc = pc    + dpc_dlam * lam
            f  = (p-pc_loc)**2/X**2 + (q/(self.M*max(pc_loc,1e-12)))**2 - 1.0
            if abs(f) < 1e-9:
                break
            df = (2*(p-pc_loc)/X**2)*(dp_dlam-dpc_dlam) \
               + 2*q/(self.M**2*max(pc_loc,1e-12)**2)*dq_dlam
            dlam = -f / (df + 1e-16)
            dlam = np.clip(dlam, -1.0, 1.0)  # limite le pas
            lam += dlam

        # Mise à jour des contraintes
        p = p_tr + dp_dlam*lam
        q = q_tr + dq_dlam*lam
        sxx_tr, syy_tr, sxy_tr = sig_tr
        dev_tr   = np.array([sxx_tr-p_tr, syy_tr-p_tr, sxy_tr])
        scale_q  = (q/q_tr) if q_tr > 1e-12 else 0.0
        sig      = np.array([p+dev_tr[0]*scale_q, p+dev_tr[1]*scale_q, dev_tr[2]*scale_q])

        # Mise à jour de l'état interne
        devp = max(0.0, lam)
        state_upd = dict(
            eps_p = state['eps_p'],
            evp   = state['evp'] + devp,
            pc    = max(self.pc0, pc + self.H_cap*devp),
            B     = min(0.999, state['B'] + self.kB*devp)
        )
        D_upd, _, _, _ = self._elastic_matrices(state_upd['B'])
        return sig, D_upd, state_upd


# ══════════════════════════════════════════════════════════════════════════════
# SOLVEUR FEM
# ══════════════════════════════════════════════════════════════════════════════

class FEMSolver:
    """Assemblage, application des CL et résolution du système K·U = F."""

    def __init__(self, mesh, constitutive, is_plastic=False, delta_T=0.0, alpha=0.0):
        self.mesh         = mesh
        self.constitutive = constitutive
        self.is_plastic   = is_plastic
        self.delta_T      = delta_T
        self.alpha        = alpha
        self.nodes, self.conn, self.bnds = mesh.generate()

        # Initialisation des états plastiques (un état par point de Gauss)
        if is_plastic:
            self.states = [constitutive.init_gp_state() for _ in range(len(self.conn)*4)]

    def assemble(self, U):
        """Assemblage global de K et du vecteur résidu R."""
        nnode = self.nodes.shape[0]
        Kglob = lil_matrix((2*nnode, 2*nnode))
        R     = np.zeros(2*nnode)
        gp_states_new = []
        ip = 0

        for e, enodes in enumerate(self.conn):
            Xe = self.nodes[enodes, :]
            Ke = np.zeros((8,8))
            Re = np.zeros(8)

            for (xi, eta, w) in gauss_points():
                _, dN_dxi = shape_Q4(xi, eta)
                J    = dN_dxi.T @ Xe
                detJ = np.linalg.det(J)
                if detJ <= 0:
                    raise RuntimeError(f"Jacobien négatif (élément {e})")

                dN_dx = dN_dxi @ np.linalg.inv(J)
                Bm    = B_matrix(dN_dx)

                # Déplacements de l'élément
                ue = np.zeros(8)
                for a in range(4):
                    ue[2*a]   = U[2*enodes[a]]
                    ue[2*a+1] = U[2*enodes[a]+1]
                deps = Bm @ ue

                # Loi de comportement
                if self.is_plastic:
                    sig, D_tan, sn = self.constitutive.stress_update(
                        deps, self.states[ip], self.delta_T, self.alpha)
                    gp_states_new.append(sn)
                else:
                    D_tan = self.constitutive.D()
                    sig   = self.constitutive.stress(deps, self.delta_T, self.alpha)

                Ke += (Bm.T @ D_tan @ Bm) * detJ * w
                Re += (Bm.T @ sig)         * detJ * w
                ip += 1

            # Assemblage dans le système global
            edofs = []
            for a in range(4):
                edofs += [2*enodes[a], 2*enodes[a]+1]
            for i in range(8):
                R[edofs[i]] += Re[i]
                for j in range(8):
                    Kglob[edofs[i], edofs[j]] += Ke[i, j]

        return csr_matrix(Kglob), R, gp_states_new

    def apply_bcs_and_loads(self, K, R, traction_nodes, F_inc):
        """Conditions aux limites et chargement nodal."""
        # Chargement vertical (vers le bas) sur les nœuds de la traverse
        for nid in traction_nodes:
            R[2*nid+1] -= F_inc

        # Conditions d'encastrement
        fixed = []
        for nid in self.bnds['bottom']:  # base bloquée en x et y
            fixed += [(2*nid, 0.0), (2*nid+1, 0.0)]
        for nid in self.bnds['left']:    # bords bloqués en x seulement
            fixed += [(2*nid, 0.0)]
        for nid in self.bnds['right']:
            fixed += [(2*nid, 0.0)]

        # Application par élimination (ligne/colonne mise à l'identité)
        K = K.tolil()
        for dof, val in fixed:
            K[dof, :] = 0; K[:, dof] = 0; K[dof, dof] = 1; R[dof] = val
        return K.tocsr(), R

    def run(self, nsteps, p_total, sleeper_span):
        """
        Résolution incrémentale en nsteps pas.
        Retourne U (déplacements), y, x (coords surface haute), top_nodes.
        """
        nnode     = self.nodes.shape[0]
        U         = np.zeros(2*nnode)
        top_nodes = self.bnds['top']
        x_top     = self.nodes[top_nodes, 0]
        xmid      = 0.5 * self.mesh.lx

        # Nœuds sous la traverse (zone centrale)
        mask           = np.abs(x_top - xmid) <= 0.5*sleeper_span
        traction_nodes = list(np.array(top_nodes)[mask])
        if not traction_nodes:
            traction_nodes = [top_nodes[np.argmin(np.abs(x_top - xmid))]]

        # Force nodale = pression × espacement inter-nœuds
        dx            = np.unique(np.diff(np.sort(x_top))).mean()
        F_node_total  = p_total * dx

        for step in range(1, nsteps+1):
            K, R, gp_new = self.assemble(np.zeros(2*nnode))  # petites déformations
            K, R = self.apply_bcs_and_loads(K, R, traction_nodes, F_node_total/nsteps)
            U    = spsolve(K, R)
            if self.is_plastic:
                self.states = gp_new

        return U, self.nodes[top_nodes, 1], self.nodes[top_nodes, 0], top_nodes


# ══════════════════════════════════════════════════════════════════════════════
# SCRIPT PRINCIPAL
# ══════════════════════════════════════════════════════════════════════════════

if __name__ == '__main__':
    # Paramètres de référence
    Lx, Ly   = 0.8, 0.35
    nx, ny   = 32, 14
    mesh     = Mesh(Lx, Ly, nx, ny)
    p_service    = 300e3   # pression de service (N/m)
    sleeper_span = 0.25    # largeur de la traverse (m)
    delta_T      = 40.0    # variation thermique (°C)
    alpha_b      = 1.2e-5  # dilatation thermique ballast
    alpha_c      = 1.0e-5  # dilatation thermique béton

    # ── Calcul ballast ────────────────────────────────────────────────────────
    ballast  = CrushableCap(150e6, 0.25, 1.2, 50e3, 80e3, 1e5, 8.0, 2.5)
    solver_b = FEMSolver(mesh, ballast, is_plastic=True, delta_T=delta_T, alpha=alpha_b)
    U_b, y_b, x_b, top_nodes = solver_b.run(6, p_service, sleeper_span)

    # ── Calcul béton ──────────────────────────────────────────────────────────
    concrete = LinearElastic(30e9, 0.2)
    solver_c = FEMSolver(mesh, concrete, is_plastic=False, delta_T=delta_T, alpha=alpha_c)
    U_c, y_c, x_c, _ = solver_c.run(6, p_service, sleeper_span)

    # ── Post-traitement : déplacements en surface ─────────────────────────────
    order  = np.argsort(x_b)
    x_plot = x_b[order]
    uy_b   = np.array([U_b[2*n+1] for n in np.array(top_nodes)[order]]) * 1e6  # µm
    uy_c   = np.array([U_c[2*n+1] for n in np.array(top_nodes)[order]]) * 1e6

    # ── Export CSV ────────────────────────────────────────────────────────────
    header = "x_m,uy_microm"
    np.savetxt("U_TOP_BALLAST.csv", np.column_stack([x_plot, uy_b]),
               delimiter=",", header=header, comments='')
    np.savetxt("U_TOP_BETON.csv",   np.column_stack([x_plot, uy_c]),
               delimiter=",", header=header, comments='')

    # ── Graphique de comparaison ──────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(9, 4.5))
    ax.plot(x_plot*100, uy_b, 'o-',  color='#e67e22', ms=4, lw=1.5,
            label='Ballast (élasto-plastique, crushable-cap)')
    ax.plot(x_plot*100, uy_c, 'x--', color='#2980b9', ms=5, lw=1.5,
            label='Béton (linéaire élastique)')
    ax.axvspan((0.5*Lx - 0.5*sleeper_span)*100, (0.5*Lx + 0.5*sleeper_span)*100,
               alpha=0.12, color='gray', label='Zone semelle (traverse)')
    ax.set_xlabel('Position x (cm)', fontsize=11)
    ax.set_ylabel('Déplacement vertical $u_y$ (µm)', fontsize=11)
    ax.set_title('Comparaison Ballast vs Béton — Surface supérieure\n'
                 r'($p_{service}$=300 kN/m, $\Delta T$=40 °C, EF Q4 plane strain)', fontsize=11)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.4)
    plt.tight_layout()
    plt.savefig("w_compare.png", dpi=150)
    plt.savefig("w_compare.pdf")

    print("=== DONE ===")
    print(f"Ballast max |uy| = {np.max(np.abs(uy_b)):.4f} µm")
    print(f"Béton   max |uy| = {np.max(np.abs(uy_c)):.4f} µm")
    print(f"Ratio ballast/béton = {np.max(np.abs(uy_b))/np.max(np.abs(uy_c)):.2f}")