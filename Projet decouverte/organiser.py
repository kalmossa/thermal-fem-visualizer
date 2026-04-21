"""
organiser.py — Classement des 180 CSV thermiques MCP
Lance ce script dans le dossier qui contient les fichiers THERM_*.csv

    py organiser.py

Actions :
  1. Classe les 180 CSV dans 4 arborescences (par MCP, épaisseur, convection, dT0)
  2. Génère les 36 graphiques T_surface comparatifs dans graphs/
  3. Copie les 3 cas sélectionnés dans data/processed/
  4. Crée note_classement.txt
"""

import os, re, shutil, csv
import matplotlib
matplotlib.use('Agg')  # pas d'affichage GUI
import matplotlib.pyplot as plt


# ── Détection des fichiers CSV ────────────────────────────────────────────────
script_dir = os.path.dirname(os.path.abspath(__file__))
all_csv    = [f for f in os.listdir(script_dir)
              if f.startswith('THERM_') and f.endswith('.csv')]
print(f"Fichiers CSV détectés : {len(all_csv)}")


# ── Parsing du nom de fichier ─────────────────────────────────────────────────
# Format : THERM_<MCP>_ep<ep>mm_hx<hx>_dT0<dt>.csv
PATTERN = re.compile(
    r'THERM_(?P<mcp>[^_]+(?:_[^_]+)*)_ep(?P<ep>\d+)mm_hx(?P<hx>\d+)_dT0(?P<dt>[^.]+)\.csv'
)

# Table de correspondance hx → valeur réelle
HX_MAP = {'5': '0.5', '10': '1', '20': '2'}

def parse_name(fname):
    """Extrait MCP, épaisseur, hx et dT0 depuis le nom de fichier."""
    m = PATTERN.match(fname)
    if not m:
        return None
    return {
        'mcp':    m.group('mcp'),
        'ep':     int(m.group('ep')),
        'hx':     HX_MAP.get(m.group('hx'), m.group('hx')),
        'hx_raw': m.group('hx'),
        'dt':     m.group('dt'),
        'fname':  fname
    }

parsed = [p for p in (parse_name(f) for f in all_csv) if p]
print(f"Fichiers parsés : {len(parsed)}")


# ── Utilitaires ───────────────────────────────────────────────────────────────
def mkd(*parts):
    """Crée un dossier (et ses parents) si nécessaire."""
    path = os.path.join(script_dir, *parts)
    os.makedirs(path, exist_ok=True)
    return path

def copyto(fname, *dest_parts):
    """Copie un fichier dans le dossier destination."""
    src = os.path.join(script_dir, fname)
    dst = os.path.join(mkd(*dest_parts), fname)
    shutil.copy2(src, dst)

# Correspondance dT0 brut → nom de dossier
DT_FOLDER = {'0': 'dT0_0', '-2': 'dT0-2', '2': 'dT0+2', '4': 'dT0+4'}
HX_FOLDER = {'0.5': 'hx05', '1': 'hx1', '2': 'hx2'}


# ── 1. Classement des fichiers ────────────────────────────────────────────────
for p in parsed:
    mcp_dir = p['mcp'].replace('_0C', '0C')  # ex: Eutectic_0C → Eutectic0C
    copyto(p['fname'], 'mcp',        mcp_dir)
    copyto(p['fname'], 'epaisseur',  f"{p['ep']}mm")
    copyto(p['fname'], 'convection', HX_FOLDER.get(p['hx'], f"hx{p['hx']}"))
    copyto(p['fname'], 'dt0',        DT_FOLDER.get(p['dt'], f"dT0_{p['dt']}"))

print("✅ Classement terminé")


# ── 2. Génération des 36 graphiques ──────────────────────────────────────────
# Couleurs fixes par MCP pour cohérence visuelle
COLORS = {
    'RT2':         '#3b82f6',
    'RT2HC':       '#f59e0b',
    'RT2_EG':      '#a78bfa',
    'Eutectic_0C': '#ef4444',
    'Hydrate_0C':  '#22c55e',
}
MCP_LIST = list(COLORS.keys())

# Toutes les combinaisons (3×3×4 = 36 graphiques)
EP_LIST = [10, 15, 30]
HX_LIST = [('5', '0.5'), ('10', '1'), ('20', '2')]
DT_LIST = ['4', '2', '0', '-2']

def read_ts(fpath):
    """Lit un CSV et retourne (temps en heures, T_surface en °C)."""
    t, ts = [], []
    with open(fpath, encoding='utf-8', errors='replace') as f:
        reader = csv.reader(f)
        next(reader)  # saute le header
        for row in reader:
            if len(row) >= 2:
                t.append(float(row[0]) / 3600)
                ts.append(float(row[1]))
    return t, ts

graphs_generated = 0
graphs_missing   = 0

for ep in EP_LIST:
    for hx_raw, hx_label in HX_LIST:
        for dt in DT_LIST:
            # Dossier de sortie du graphique
            dt_folder = DT_FOLDER.get(dt, f'dT0_{dt}')
            out_dir   = mkd('graphs', f'ep{ep}mm', f'hx{hx_raw}', dt_folder)
            out_file  = os.path.join(out_dir, f'graph_ep{ep}mm_hx{hx_raw}_dT0{dt}.png')

            fig, ax = plt.subplots(figsize=(10, 4.5))
            found   = False

            # Trace une courbe par MCP si le fichier existe
            for mcp in MCP_LIST:
                fname = f'THERM_{mcp}_ep{ep}mm_hx{hx_raw}_dT0{dt}.csv'
                fpath = os.path.join(script_dir, fname)
                if os.path.exists(fpath):
                    t_h, ts = read_ts(fpath)
                    ax.plot(t_h, ts, label=mcp, color=COLORS[mcp], linewidth=1.5)
                    found = True
                else:
                    graphs_missing += 1

            if found:
                dt_display = f'+{dt}' if dt not in ('0', '-2') else dt
                ax.set_title(f'Comparaison T_surface selon MCP — ep={ep}mm, hx={hx_label}, dT0={dt_display}°C')
                ax.set_xlabel('Temps (heures)')
                ax.set_ylabel('T_surface (°C)')
                ax.legend(fontsize=9)
                ax.grid(True, alpha=0.3)
                fig.tight_layout()
                fig.savefig(out_file, dpi=80)  # dpi réduit pour limiter la taille
                graphs_generated += 1
            plt.close(fig)

print(f"✅ {graphs_generated} graphiques générés ({graphs_missing} fichiers manquants ignorés)")


# ── 3. Extraction des 3 cas pour la page web ──────────────────────────────────
# A) Cas standard : ep15mm, hx=1, dT0=+2
# B) Cas extrême  : ep10mm, hx=0.5, dT0=-2
# C) Cas optimal  : ep30mm, hx=2, dT0=+4
CASES = [
    ('A_standard', 'ep15mm', 'hx10', '2'),
    ('B_extreme',  'ep10mm', 'hx5',  '-2'),
    ('C_optimal',  'ep30mm', 'hx20', '4'),
]

proc_dir = mkd('data', 'processed')

for label, ep_str, hx_str, dt in CASES:
    for mcp in MCP_LIST:
        fname = f'THERM_{mcp}_{ep_str}_{hx_str}_dT0{dt}.csv'
        src   = os.path.join(script_dir, fname)
        if os.path.exists(src):
            shutil.copy2(src, os.path.join(proc_dir, f'{label}_{fname}'))

print("✅ data/processed/ créé")


# ── 4. Note de classement ────────────────────────────────────────────────────
NOTE = """NOTE DE CLASSEMENT — THERM_*.csv
=================================
Auteur : Elias Lallouet
Date   : Mars 2026

STRUCTURE DES NOMS DE FICHIERS
-------------------------------
Format : THERM_<MCP>_ep<épaisseur>mm_hx<convection>_dT0<condition>.csv

Exemples :
  THERM_RT2_ep15mm_hx10_dT02.csv
    → MCP = RT2, épaisseur = 15 mm, convection hx = 1 W/m²K, dT0 = +2°C
  THERM_Hydrate_0C_ep30mm_hx20_dT04.csv
    → MCP = Hydrate_0C, épaisseur = 30 mm, convection hx = 2 W/m²K, dT0 = +4°C

CORRESPONDANCES DES CODES
--------------------------
Convection (hx) :
  hx5  → h = 0.5 W/m²K (faible)
  hx10 → h = 1.0 W/m²K (standard)
  hx20 → h = 2.0 W/m²K (fort)

Condition initiale (dT0) :
  dT04  → T_initiale = +4°C
  dT02  → T_initiale = +2°C
  dT00  → T_initiale =  0°C
  dT0-2 → T_initiale = -2°C

ARBORESCENCES CRÉÉES
--------------------
mcp/            → 5 dossiers (36 fichiers chacun)
epaisseur/      → 3 dossiers (60 fichiers chacun)
convection/     → 3 dossiers (60 fichiers chacun)
dt0/            → 4 dossiers (45 fichiers chacun)
graphs/         → 36 graphiques PNG organisés par ep/hx/dT0
data/processed/ → 3 cas sélectionnés pour la page web

SÉLECTION POUR LA PAGE WEB
---------------------------
A) Cas standard : ep=15mm, hx=1, dT0=+2
B) Cas extrême  : ep=10mm, hx=0.5, dT0=-2
C) Cas optimal  : ep=30mm, hx=2, dT0=+4

COMPATIBILITÉ FUTURES SIMULATIONS
----------------------------------
Seules les colonnes garanties sont utilisées : t(s), Tsurface(°C), Flux(W), E(J), f_liq
Toute colonne supplémentaire dans de futurs CSV est ignorée automatiquement.
"""

with open(os.path.join(script_dir, 'note_classement.txt'), 'w', encoding='utf-8') as f:
    f.write(NOTE)

print("✅ note_classement.txt créé")
print("\n=== TERMINÉ ===")
print("Ouvre Projet decouverte/web/index.html dans Chrome.")