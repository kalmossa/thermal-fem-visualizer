# thermal-fem-visualizer

Application web de simulation thermo-mécanique par éléments finis pour matériaux ferroviaires. Le solveur Python calcule en temps réel les déplacements nodaux sous charge combinée mécanique et thermique — les résultats sont visualisés via une interface HTML5 connectée à une API REST Flask.

Le projet inclut aussi une analyse de matériaux à changement de phase (PCM) avec classement automatique de 180 fichiers de simulation et génération de 36 graphiques comparatifs.

---

## Stack

- **Backend** — Python, Flask, NumPy, SciPy, Matplotlib
- **Frontend** — HTML5, CSS3, JavaScript vanilla, Chart.js, Canvas 2D
- **API** — REST, JSON, Fetch API async
- **FEM** — Éléments Q4 bilinéaires, intégration de Gauss 2×2, plasticité crushable-cap, effets thermiques

---

## Lancer le projet

```powershell
py -m pip install flask numpy scipy matplotlib
py app.py
```

Ouvrir **http://localhost:5000** dans Chrome.

Sans serveur : double-clic sur `fem_visualizer.html` (mode statique avec fallback).

---

## Classement MCP + graphiques

```powershell
cd "Projet decouverte"
py organiser.py
```

Génère automatiquement les 4 arborescences de classement, les 36 graphiques PNG et extrait les 3 cas pour la page web.

---

## Cloner et faire tourner

```powershell
git clone https://github.com/kalmossa/thermal-fem-visualizer.git
cd thermal-fem-visualizer
py -m pip install flask numpy scipy matplotlib
py app.py
```

Pour le MCP :
```powershell
cd "Projet decouverte"
py organiser.py
```

---

*Elias Lallouet — 2026*
