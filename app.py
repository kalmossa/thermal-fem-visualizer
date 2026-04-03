"""
API Flask — FEM Ballast vs Béton
Lance le vrai calcul FEM avec les paramètres envoyés par le frontend.

Usage :
    pip install flask flask-cors numpy scipy matplotlib
    python app.py
    → http://localhost:5000
"""

from flask import Flask, request, jsonify, send_from_directory
# CORS géré manuellement ci-dessous
import numpy as np
import os

# ── Import du solveur FEM ─────────────────────────────────────────────────────
# On importe les classes du fichier Python existant sans lancer son __main__
import importlib.util, types

def load_fem_module():
    path = os.path.join(os.path.dirname(__file__), "fem_ballast_beton.py")
    spec = importlib.util.spec_from_file_location("fem_solver", path)
    mod  = importlib.util.module_from_spec(spec)
    # On court-circuite le __main__ pour ne pas lancer le script complet
    mod.__name__ = "fem_solver"
    spec.loader.exec_module(mod)
    return mod

fem = load_fem_module()
Mesh         = fem.Mesh
LinearElastic= fem.LinearElastic
CrushableCap = fem.CrushableCap
FEMSolver    = fem.FEMSolver

# ── App Flask ─────────────────────────────────────────────────────────────────
app = Flask(__name__, static_folder=".")

# CORS manuel — autorise le HTML local à appeler l'API
@app.after_request
def add_cors(response):
    response.headers["Access-Control-Allow-Origin"]  = "*"
    response.headers["Access-Control-Allow-Methods"] = "POST, GET, OPTIONS"
    response.headers["Access-Control-Allow-Headers"] = "Content-Type"
    return response

@app.route("/api/simulate", methods=["OPTIONS"])
def options():
    return "", 204

# Maillage fixe (32×14, 0.8m × 0.35m) — on ne le recrée pas à chaque requête
MESH = Mesh(0.8, 0.35, 32, 14)

@app.route("/")
def index():
    """Sert directement le frontend HTML."""
    return send_from_directory(".", "fem_visualizer.html")

@app.route("/api/simulate", methods=["POST"])
def simulate():
    """
    Lance le calcul FEM avec les paramètres reçus en JSON.

    Corps attendu :
    {
        "p"      : 300,     // pression en kN/m
        "delta_T": 40,      // variation thermique en °C
        "E_b"    : 150,     // module Young ballast en MPa
        "alpha_b": 1.2e-5,  // coeff thermique ballast
        "E_c"    : 30000,   // module Young béton en MPa
        "alpha_c": 1.0e-5   // coeff thermique béton
    }

    Réponse :
    {
        "x_m"      : [...],  // positions x en mètres (33 valeurs)
        "uy_ballast": [...],  // déplacements verticaux ballast en µm
        "uy_beton"  : [...],  // déplacements verticaux béton en µm
        "max_ballast": float,
        "max_beton"  : float,
        "ratio"      : float
    }
    """
    data = request.get_json(force=True)

    # ── Paramètres avec valeurs par défaut ────────────────────────────────────
    p_kn_m   = float(data.get("p",       300))       # kN/m
    delta_T  = float(data.get("delta_T", 40.0))      # °C
    E_b_mpa  = float(data.get("E_b",     150))       # MPa
    alpha_b  = float(data.get("alpha_b", 1.2e-5))
    E_c_mpa  = float(data.get("E_c",     30000))     # MPa
    alpha_c  = float(data.get("alpha_c", 1.0e-5))

    # Conversion unités → SI
    p_service   = p_kn_m * 1e3          # Pa·m = N/m
    E_b         = E_b_mpa * 1e6         # Pa
    E_c         = E_c_mpa * 1e6         # Pa
    sleeper_span = 0.25                  # m (fixe)

    try:
        # ── Calcul Ballast (élasto-plastique crushable-cap) ───────────────────
        ballast  = CrushableCap(E_b, 0.25, 1.2, 50e3, 80e3, 1e5, 8.0, 2.5)
        solver_b = FEMSolver(MESH, ballast, is_plastic=True,
                             delta_T=delta_T, alpha=alpha_b)
        U_b, y_b, x_b, top_nodes = solver_b.run(6, p_service, sleeper_span)

        # ── Calcul Béton (linéaire élastique) ─────────────────────────────────
        concrete = LinearElastic(E_c, 0.2)
        solver_c = FEMSolver(MESH, concrete, is_plastic=False,
                             delta_T=delta_T, alpha=alpha_c)
        U_c, y_c, x_c, _ = solver_c.run(6, p_service, sleeper_span)

        # ── Post-traitement ───────────────────────────────────────────────────
        order  = np.argsort(x_b)
        x_plot = x_b[order]
        uy_b   = np.array([U_b[2*n+1] for n in np.array(top_nodes)[order]]) * 1e6
        uy_c   = np.array([U_c[2*n+1] for n in np.array(top_nodes)[order]]) * 1e6

        max_b = float(np.max(np.abs(uy_b)))
        max_c = float(np.max(np.abs(uy_c)))

        return jsonify({
            "x_m":        x_plot.tolist(),
            "uy_ballast": uy_b.tolist(),
            "uy_beton":   uy_c.tolist(),
            "max_ballast": max_b,
            "max_beton":   max_c,
            "ratio":       round(max_b / max_c, 3) if max_c else 0,
            "params": {
                "p": p_kn_m, "delta_T": delta_T,
                "E_b": E_b_mpa, "alpha_b": alpha_b,
                "E_c": E_c_mpa, "alpha_c": alpha_c,
            }
        })

    except Exception as e:
        return jsonify({"error": str(e)}), 500


if __name__ == "__main__":
    print("=" * 55)
    print("  FEM Visualizer — API Flask")
    print("  http://localhost:5000")
    print("  POST /api/simulate  → calcul FEM")
    print("=" * 55)
    app.run(debug=True, port=5000)
