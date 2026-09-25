"""Régénère le jeu de données de repli embarqué dans `fem_visualizer.html`.

Quand la page est ouverte sans serveur (double-clic sur le fichier), elle
affiche ce jeu de données plutôt qu'une page vide. Ce sont de vrais résultats
du solveur aux paramètres par défaut — la page l'annonce explicitement et
désactive toute prétention à recalculer.

À relancer après toute modification du solveur, sinon la démo hors-ligne
affiche les résultats de l'ancienne version :

    py tools/make_offline_data.py
"""
import io
import json
import os
import sys

RACINE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, RACINE)

from app import app  # noqa: E402  (import après ajustement du sys.path)

PAGE = os.path.join(RACINE, "fem_visualizer.html")
DEBUT = "const OFFLINE = "
FIN = ";\n"

CHAMPS = ("vm_ballast", "vm_beton", "evp_ballast", "uy_nodal_ballast", "uy_nodal_beton")
RESUME = ("uy_min_ballast", "uy_max_ballast", "uy_min_beton", "uy_max_beton",
          "bowl_ballast", "bowl_beton", "ratio_bowl", "evp_max_pct", "B_max", "pc_max_kpa")


def main():
    client = app.test_client()
    reponse = client.post("/api/run-simulation", json={})
    if reponse.status_code != 200:
        raise SystemExit(f"le calcul a échoué : {reponse.get_json()}")
    d = reponse.get_json()

    charge = {
        "params": d["params"],
        "mesh": d["mesh"],
        "x_cm": d["x_cm"],
        "uy_ballast": d["uy_ballast"],
        "uy_beton": d["uy_beton"],
        "summary": {k: d["summary"][k] for k in RESUME},
        "fields": {k: d["fields"][k] for k in CHAMPS},
    }
    charge_json = json.dumps(charge, separators=(",", ":"))

    page = io.open(PAGE, encoding="utf-8").read()
    i = page.index(DEBUT) + len(DEBUT)
    j = page.index(FIN, i)
    io.open(PAGE, "w", encoding="utf-8").write(page[:i] + charge_json + page[j:])

    s = charge["summary"]
    print(f"{len(charge_json)} octets injectés dans fem_visualizer.html")
    print(f"  cuvette ballast {s['bowl_ballast']:.1f} µm | béton {s['bowl_beton']:.2f} µm"
          f" | rapport {s['ratio_bowl']:.0f}")


if __name__ == "__main__":
    main()
