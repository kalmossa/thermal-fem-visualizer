"""Génère les figures du README à partir du solveur.

Des captures d'écran se périment sans prévenir ; ces figures se refont en une
commande après chaque modification du solveur :

    py tools/make_figures.py
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

RACINE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, RACINE)

from fem_ballast_beton import (  # noqa: E402
    BALLAST_PARAMS, CONCRETE_PARAMS, CrushableCap, FEMSolver, LinearElastic, Mesh,
)

SORTIE = os.path.join(RACINE, "docs")


def carte(ax, mesh, valeurs, titre, unite, cmap):
    """Trace un champ constant par élément sur le maillage."""
    nx, ny = mesh.nx, mesh.ny
    grille = np.asarray(valeurs).reshape(ny, nx)
    im = ax.imshow(grille, origin="lower", cmap=cmap, aspect="auto",
                   extent=[0, mesh.lx * 100, 0, mesh.ly * 100])
    ax.set_title(titre, fontsize=10)
    ax.set_xlabel("x (cm)", fontsize=8)
    ax.set_ylabel("y (cm)", fontsize=8)
    ax.tick_params(labelsize=7)
    cb = plt.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
    cb.set_label(unite, fontsize=8)
    cb.ax.tick_params(labelsize=7)

    # emprise de la traverse
    xm, demi = mesh.lx * 50, 12.5
    ax.axvspan(xm - demi, xm + demi, fc="none", ec="w", ls="--", lw=1.1, alpha=.7)


def main():
    os.makedirs(SORTIE, exist_ok=True)
    mesh = Mesh(0.8, 0.35, 32, 14)
    p_service, span, dT = 300e3, 0.25, 40.0

    sb = FEMSolver(mesh, CrushableCap(**BALLAST_PARAMS), is_plastic=True,
                   delta_T=dT, alpha=1.2e-5)
    U_b, *_ = sb.run(6, p_service, span)

    sc = FEMSolver(mesh, LinearElastic(30e9, CONCRETE_PARAMS["nu"]),
                   delta_T=dT, alpha=1.0e-5)
    U_c, *_ = sc.run(6, p_service, span)

    # uy moyen par élément, pour rester homogène avec les autres champs
    uy = U_b[1::2] * 1e6
    nx, ny = mesh.nx, mesh.ny
    uy_elem = np.array([
        (uy[j * (nx + 1) + i] + uy[j * (nx + 1) + i + 1]
         + uy[(j + 1) * (nx + 1) + i] + uy[(j + 1) * (nx + 1) + i + 1]) / 4
        for j in range(ny) for i in range(nx)
    ])

    fig, axes = plt.subplots(1, 3, figsize=(15, 3.6))
    carte(axes[0], mesh, uy_elem, "Déplacement vertical $u_y$ — Ballast", "µm", "RdYlBu_r")
    carte(axes[1], mesh, sb.stress_field(U_b), "Contrainte de von Mises — Ballast", "kPa", "inferno")
    carte(axes[2], mesh, sb.gauss_field("evp") * 100,
          "Déformation volumique plastique cumulée", "%", "viridis")
    fig.suptitle(f"Ballast sous traverse — p = 300 kN/m, ΔT = +{dT:.0f} °C, "
                 f"maillage Q4 {nx}×{ny}", fontsize=11)
    fig.tight_layout()
    chemin = os.path.join(SORTIE, "champs.png")
    fig.savefig(chemin, dpi=130)
    print(f"écrit : {os.path.relpath(chemin, RACINE)}")

    _, uy_b = sb.top_profile(U_b)
    _, uy_c = sc.top_profile(U_c)
    print(f"  ballast : u_y ∈ [{uy_b.min():+.1f}, {uy_b.max():+.1f}] µm, "
          f"cuvette {uy_b[0] - uy_b.min():.0f} µm")
    print(f"  béton   : u_y ∈ [{uy_c.min():+.1f}, {uy_c.max():+.1f}] µm, "
          f"cuvette {uy_c[0] - uy_c.min():.2f} µm")


if __name__ == "__main__":
    main()
