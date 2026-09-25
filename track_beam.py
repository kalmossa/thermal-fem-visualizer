"""Flexion d'une traverse sur appui élastique (modèle de Winkler).

Le problème
-----------
Une traverse repose sur le ballast et reçoit les efforts des deux files de rail,
appliqués aux deux sièges distants de l'écartement (1,435 m). Le ballast est
modélisé par un appui élastique réparti de module k (N/m de flèche, par mètre de
traverse), et la semelle sous rail (EVA) par un ressort ponctuel entre le rail
et la traverse.

Deux résolutions cohabitent ici, et c'est voulu :

1. `deflection` — solution analytique de Hetényi pour une **poutre infinie** :

       w(x) = (P·β / 2k)·e^(−β|x|)·(cos β|x| + sin β|x|),   β = (k / 4EI)^(1/4)

   Elle sert de référence exacte et de cas de contrôle.

2. `SleeperBeam` — éléments finis de poutre d'Euler-Bernoulli (fonctions de
   Hermite, 2 degrés de liberté par nœud) sur fondation de Winkler. C'est ce que
   l'application utilise, parce qu'une traverse de 2,5 m pour une longueur
   caractéristique β⁻¹ de l'ordre du mètre n'est **pas** une poutre infinie :
   βL ≈ 2,5, les extrémités libres comptent. Le test
   `test_poutre_longue_rejoint_hetenyi` vérifie que le modèle numérique
   retrouve la solution analytique quand on allonge la poutre.

Ce que ce modèle apporte par rapport à une somme de cloches gaussiennes
----------------------------------------------------------------------
- **L'équilibre est vérifié** : la réaction du ballast intégrée sur la traverse
  égale exactement la charge appliquée. C'est le contrôle que ne passe aucune
  forme choisie a priori.
- **EI intervient**, via β. La rigidité de flexion fixe la largeur du bassin,
  donc la façon dont la charge se répartit sur le ballast.
- **Le moment fléchissant est disponible.** C'est lui qui dimensionne la
  traverse — positif sous les sièges, négatif au centre — et il ne se déduit
  pas d'une forme d'allure plausible.
- **Le soulèvement apparaît quand il doit apparaître, et pas avant.** La
  solution infinie présente des lobes de soulèvement en βx ≈ π ; une traverse
  réelle (βL ≈ 2,9) est trop courte pour les développer et s'enfonce partout,
  le centre moins que les sièges. Le modèle reproduit les deux régimes et la
  transition entre les deux vers βL ≈ 4,5 — là où une gaussienne ne peut, par
  construction, jamais changer de signe.

Conventions : flèche négative vers le bas, charges descendantes positives.
"""
import numpy as np

# Traverse monobloc béton : section 0,25 × 0,20 m, E ≈ 35 GPa.
EI_TRAVERSE = 35e9 * 0.25 * 0.20 ** 3 / 12.0       # ≈ 5,83 MN·m²
LONGUEUR_TRAVERSE = 2.50                            # m
ECARTEMENT_STANDARD = 1.435                         # m
K_BALLAST = 40e6                                    # N/m² (module de réaction × largeur)


# ─────────────────────────────────────────────────────────────────────────────
# Solution analytique : poutre infinie (référence)
# ─────────────────────────────────────────────────────────────────────────────
def beta(k, ei):
    """Inverse de la longueur caractéristique, en 1/m."""
    if k <= 0 or ei <= 0:
        raise ValueError("k et EI doivent être strictement positifs")
    return (k / (4.0 * ei)) ** 0.25


def deflection(x, p, k, ei, x0=0.0):
    """Flèche d'une poutre infinie sous charge ponctuelle, en m (négative = bas)."""
    b = beta(k, ei)
    bx = b * np.abs(np.asarray(x, dtype=float) - x0)
    return -(p * b / (2.0 * k)) * np.exp(-bx) * (np.cos(bx) + np.sin(bx))


def moment_infini(x, p, k, ei, x0=0.0):
    """Moment fléchissant de la poutre infinie, en N·m."""
    b = beta(k, ei)
    bx = b * np.abs(np.asarray(x, dtype=float) - x0)
    return (p / (4.0 * b)) * np.exp(-bx) * (np.cos(bx) - np.sin(bx))


# ─────────────────────────────────────────────────────────────────────────────
# Éléments finis de poutre sur fondation élastique
# ─────────────────────────────────────────────────────────────────────────────
class SleeperBeam:
    """Poutre d'Euler-Bernoulli à extrémités libres sur fondation de Winkler.

    Le système n'est pas singulier bien qu'aucun appui ne soit imposé : c'est la
    fondation elle-même qui équilibre la poutre. Un corps rigide en translation
    ou en rotation y travaille, donc la matrice reste inversible.
    """

    def __init__(self, longueur=LONGUEUR_TRAVERSE, ei=EI_TRAVERSE, k=K_BALLAST,
                 n_elem=120):
        if longueur <= 0 or n_elem < 2:
            raise ValueError("longueur et n_elem doivent être positifs")
        self.L = float(longueur)
        self.ei = float(ei)
        self.k = float(k)
        self.n_elem = int(n_elem)
        self.x = np.linspace(-self.L / 2, self.L / 2, self.n_elem + 1)
        self.he = self.L / self.n_elem
        self._assemble()

    # ── matrices élémentaires ────────────────────────────────────────────
    @staticmethod
    def _k_flexion(ei, h):
        return (ei / h ** 3) * np.array([
            [12.0,   6 * h,  -12.0,   6 * h],
            [6 * h,  4 * h * h, -6 * h, 2 * h * h],
            [-12.0, -6 * h,   12.0,  -6 * h],
            [6 * h,  2 * h * h, -6 * h, 4 * h * h],
        ])

    @staticmethod
    def _k_fondation(k, h):
        """Matrice de fondation cohérente, obtenue en intégrant k·N·Nᵀ."""
        return (k * h / 420.0) * np.array([
            [156.0,  22 * h,   54.0,  -13 * h],
            [22 * h, 4 * h * h, 13 * h, -3 * h * h],
            [54.0,   13 * h,  156.0,  -22 * h],
            [-13 * h, -3 * h * h, -22 * h, 4 * h * h],
        ])

    def _assemble(self):
        ndof = 2 * (self.n_elem + 1)
        ke = self._k_flexion(self.ei, self.he) + self._k_fondation(self.k, self.he)
        K = np.zeros((ndof, ndof))
        for e in range(self.n_elem):
            d = np.array([2 * e, 2 * e + 1, 2 * e + 2, 2 * e + 3])
            K[np.ix_(d, d)] += ke
        self.K = K

    # ── résolution ───────────────────────────────────────────────────────
    def solve(self, charges):
        """`charges` : liste de (position en m, force descendante en N).

        Chaque charge est reportée sur les deux nœuds de son élément par les
        fonctions de forme de Hermite, ce qui évite de contraindre la
        discrétisation à faire tomber un nœud sur chaque siège de rail.
        """
        F = np.zeros(2 * (self.n_elem + 1))
        for x0, p in charges:
            x0 = float(np.clip(x0, self.x[0], self.x[-1]))
            e = min(int((x0 - self.x[0]) / self.he), self.n_elem - 1)
            s = (x0 - self.x[e]) / self.he            # coordonnée locale ∈ [0, 1]
            h = self.he
            n = np.array([
                1 - 3 * s ** 2 + 2 * s ** 3,
                h * (s - 2 * s ** 2 + s ** 3),
                3 * s ** 2 - 2 * s ** 3,
                h * (-s ** 2 + s ** 3),
            ])
            F[2 * e:2 * e + 4] += -p * n              # descendant ⇒ flèche négative

        u = np.linalg.solve(self.K, F)
        w = u[0::2]
        theta = u[1::2]
        return w, theta, u

    def moments(self, u):
        """Moment fléchissant aux nœuds, en N·m (moyenne des deux éléments voisins)."""
        h = self.he
        somme = np.zeros(self.n_elem + 1)
        compte = np.zeros(self.n_elem + 1)
        for e in range(self.n_elem):
            d = u[2 * e:2 * e + 4]
            for s, noeud in ((0.0, e), (1.0, e + 1)):
                d2n = np.array([
                    (-6 + 12 * s) / h ** 2,
                    (-4 + 6 * s) / h,
                    (6 - 12 * s) / h ** 2,
                    (-2 + 6 * s) / h,
                ])
                somme[noeud] += self.ei * float(d2n @ d)
                compte[noeud] += 1
        return somme / np.maximum(compte, 1)

    def deflection_at(self, w, theta, x0):
        """Flèche interpolée en un point quelconque, par les fonctions de Hermite."""
        x0 = float(np.clip(x0, self.x[0], self.x[-1]))
        e = min(int((x0 - self.x[0]) / self.he), self.n_elem - 1)
        s = (x0 - self.x[e]) / self.he
        h = self.he
        n = np.array([
            1 - 3 * s ** 2 + 2 * s ** 3,
            h * (s - 2 * s ** 2 + s ** 3),
            3 * s ** 2 - 2 * s ** 3,
            h * (-s ** 2 + s ** 3),
        ])
        d = np.array([w[e], theta[e], w[e + 1], theta[e + 1]])
        return float(n @ d)


# ─────────────────────────────────────────────────────────────────────────────
# Cas d'application : traverse sous deux files de rail
# ─────────────────────────────────────────────────────────────────────────────
def twin_rail(f_left_kn=64.0, f_right_kn=64.0, k_dff_kn_mm=27.9,
              gauge=ECARTEMENT_STANDARD, ei=EI_TRAVERSE, k_ballast=K_BALLAST,
              longueur=LONGUEUR_TRAVERSE, n=161):
    """Réponse d'une traverse chargée par les deux files de rail.

    La semelle sous rail est un ressort ponctuel de raideur `k_dff_kn_mm` : le
    rail descend de la flèche de la traverse **plus** l'écrasement de la semelle
    sous sa propre charge. Les deux raideurs agissent donc sur des grandeurs
    distinctes, là où un modèle à un seul étage ne peut produire qu'une fraction
    fixe de la flèche.
    """
    if gauge >= longueur:
        raise ValueError("l'écartement doit être inférieur à la longueur de la traverse")

    poutre = SleeperBeam(longueur=longueur, ei=ei, k=k_ballast, n_elem=max(n - 1, 20))
    xl, xr = -gauge / 2.0, gauge / 2.0
    pl, pr = f_left_kn * 1e3, f_right_kn * 1e3

    w, theta, u = poutre.solve([(xl, pl), (xr, pr)])
    wl, tl, _ = poutre.solve([(xl, pl)])
    wr, tr, _ = poutre.solve([(xr, pr)])
    m = poutre.moments(u)

    k_pad = k_dff_kn_mm * 1e6                       # kN/mm → N/m
    ecr_l = -pl / k_pad
    ecr_r = -pr / k_pad
    w_traverse_l = poutre.deflection_at(w, theta, xl)
    w_traverse_r = poutre.deflection_at(w, theta, xr)

    q = k_ballast * w                                # réaction répartie, N/m
    b = beta(k_ballast, ei)
    # Une traverse réelle (βL ≈ 2,9) est trop courte pour que la flèche change
    # de signe : tout descend, le centre moins que les sièges. Les lobes de
    # soulèvement de la solution infinie n'apparaissent qu'au-delà de βL ≈ 2π.
    # On rapporte donc le point le moins enfoncé, et on dit explicitement s'il
    # y a soulèvement ou non plutôt que d'appeler « soulèvement » un minimum
    # d'enfoncement.
    i_haut = int(np.argmax(w))
    souleve = bool(w.max() > 0)

    return {
        "x_cm": (poutre.x * 100).round(3).tolist(),
        "w_traverse": (w * 1e6).round(3).tolist(),
        "w_left": (wl * 1e6).round(3).tolist(),
        "w_right": (wr * 1e6).round(3).tolist(),
        "reaction_kN_m": (q / 1e3).round(4).tolist(),
        "moment_kNm": (m / 1e3).round(4).tolist(),
        "summary": {
            "w_rail_gauche_um": (w_traverse_l + ecr_l) * 1e6,
            "w_rail_droit_um": (w_traverse_r + ecr_r) * 1e6,
            "w_traverse_gauche_um": w_traverse_l * 1e6,
            "w_traverse_droit_um": w_traverse_r * 1e6,
            "ecrasement_gauche_um": ecr_l * 1e6,
            "ecrasement_droit_um": ecr_r * 1e6,
            "w_haut_um": float(w[i_haut]) * 1e6,
            "w_haut_x_cm": float(poutre.x[i_haut]) * 100,
            "souleve": souleve,
            "differentiel_siege_centre_um": (w_traverse_l - float(np.interp(0.0, poutre.x, w))) * 1e6,
            "longueur_caracteristique_m": 1.0 / b,
            "beta_L": b * longueur,
            "k_ballast_MPa": k_ballast / 1e6,
            "moment_max_kNm": float(np.max(np.abs(m))) / 1e3,
            "moment_siege_kNm": float(np.interp(xl, poutre.x, m)) / 1e3,
            "moment_centre_kNm": float(np.interp(0.0, poutre.x, m)) / 1e3,
            "reaction_totale_kN": float(np.trapezoid(-q, poutre.x)) / 1e3,
        },
        "params": {
            "f_left_kn": f_left_kn, "f_right_kn": f_right_kn,
            "k_dff_kn_mm": k_dff_kn_mm, "gauge": gauge,
            "ei_MNm2": ei / 1e6, "k_ballast_MPa": k_ballast / 1e6,
            "longueur": longueur,
        },
    }


if __name__ == "__main__":
    r = twin_rail()
    s = r["summary"]
    print("=" * 66)
    print(f"  longueur caractéristique β⁻¹ = {s['longueur_caracteristique_m']:.3f} m"
          f"   (βL = {s['beta_L']:.2f})")
    print(f"  flèche de la traverse au siège = {s['w_traverse_gauche_um']:+9.1f} µm")
    print(f"  écrasement de la semelle       = {s['ecrasement_gauche_um']:+9.1f} µm")
    print(f"  flèche du rail                 = {s['w_rail_gauche_um']:+9.1f} µm")
    print(f"  point le moins enfoncé         = {s['w_haut_um']:+9.1f} µm "
          f"à x = {s['w_haut_x_cm']:+.0f} cm"
          f"  ({'soulèvement' if s['souleve'] else 'pas de soulèvement'})")
    print(f"  différentiel siège − centre    = {s['differentiel_siege_centre_um']:+9.1f} µm")
    print(f"  moment au siège / au centre    = {s['moment_siege_kNm']:+6.2f} / "
          f"{s['moment_centre_kNm']:+6.2f} kN·m")
    print(f"  réaction totale du ballast     = {s['reaction_totale_kN']:.2f} kN "
          f"(charge appliquée {r['params']['f_left_kn'] + r['params']['f_right_kn']:.0f} kN)")
    print("=" * 66)
