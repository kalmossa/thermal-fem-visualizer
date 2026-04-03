"""
organiser.py
------------
Lance ce script dans le dossier qui contient les 180 fichiers THERM_*.csv
    py organiser.py

Il va :
1. Créer toute l'arborescence de classement (mcp/, epaisseur/, convection/, dt0/)
2. Copier chaque fichier dans les bons dossiers
3. Générer les 36 graphiques T_surface dans graphs/
4. Extraire les 3 fichiers utiles dans data/processed/
5. Créer une note de classement (note_classement.txt)
"""

import os, re, shutil, csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ── 1. Détection des fichiers CSV ─────────────────────────────────────────────
script_dir = os.path.dirname(os.path.abspath(__file__))
all_csv = [f for f in os.listdir(script_dir) if f.startswith('THERM_') and f.endswith('.csv')]
print(f"Fichiers CSV détectés : {len(all_csv)}")

# ── 2. Parsing du nom de fichier ──────────────────────────────────────────────
# THERM_<MCP>_ep<ep>mm_hx<hx>_dT0<dt>.csv
PATTERN = re.compile(
    r'THERM_(?P<mcp>[^_]+(?:_[^_]+)*)_ep(?P<ep>\d+)mm_hx(?P<hx>\d+)_dT0(?P<dt>[^.]+)\.csv'
)

def parse_name(fname):
    m = PATTERN.match(fname)
    if not m:
        return None
    mcp = m.group('mcp')
    ep  = int(m.group('ep'))
    hx_raw = m.group('hx')
    # hx5 → 0.5, hx10 → 1, hx20 → 2
    hx_map = {'5': '0.5', '10': '1', '20': '2'}
    hx  = hx_map.get(hx_raw, hx_raw)
    dt  = m.group('dt')   # ex: '4', '2', '0', '-2'
    return {'mcp': mcp, 'ep': ep, 'hx': hx, 'hx_raw': hx_raw, 'dt': dt, 'fname': fname}

parsed = [p for p in (parse_name(f) for f in all_csv) if p]
print(f"Fichiers parsés : {len(parsed)}")

# ── 3. Arborescences de classement ────────────────────────────────────────────
def mkd(*parts):
    path = os.path.join(script_dir, *parts)
    os.makedirs(path, exist_ok=True)
    return path

def copyto(fname, *dest_parts):
    src = os.path.join(script_dir, fname)
    dst = os.path.join(mkd(*dest_parts), fname)
    shutil.copy2(src, dst)

for p in parsed:
    # Par MCP
    mcp_dir = p['mcp'].replace('_0C', '0C').replace('_EG', '_EG')
    copyto(p['fname'], 'mcp', mcp_dir)

    # Par épaisseur
    copyto(p['fname'], 'epaisseur', f"{p['ep']}mm")

    # Par convection
    hx_folder = {'0.5': 'hx05', '1': 'hx1', '2': 'hx2'}.get(p['hx'], f"hx{p['hx']}")
    copyto(p['fname'], 'convection', hx_folder)

    # Par dT0
    dt_raw = p['dt']
    if dt_raw == '0':       dt_folder = 'dT0_0'
    elif dt_raw == '-2':    dt_folder = 'dT0-2'
    elif dt_raw == '2':     dt_folder = 'dT0+2'
    elif dt_raw == '4':     dt_folder = 'dT0+4'
    else:                   dt_folder = f'dT0_{dt_raw}'
    copyto(p['fname'], 'dt0', dt_folder)

print("✅ Classement terminé")

# ── 4. Génération des 36 graphiques ──────────────────────────────────────────
COLORS = {
    'RT2':         '#3b82f6',
    'RT2HC':       '#f59e0b',
    'RT2_EG':      '#a78bfa',
    'Eutectic_0C': '#ef4444',
    'Hydrate_0C':  '#22c55e',
}
EP_LIST = [10, 15, 30]
HX_LIST = [('5', '0.5'), ('10', '1'), ('20', '2')]
DT_LIST = ['4', '2', '0', '-2']
MCP_LIST = list(COLORS.keys())

def read_ts(fpath):
    """Retourne (t_heures[], Tsurface[]) depuis un CSV."""
    t, ts = [], []
    with open(fpath, encoding='utf-8', errors='replace') as f:
        reader = csv.reader(f)
        next(reader)  # header
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
            # Dossier de sortie
            dt_folder = {'4':'dT0+4','2':'dT0+2','0':'dT0_0','-2':'dT0-2'}.get(dt, f'dT0_{dt}')
            out_dir = mkd('graphs', f'ep{ep}mm', f'hx{hx_raw}', dt_folder)
            out_file = os.path.join(out_dir, f'graph_ep{ep}mm_hx{hx_raw}_dT0{dt}.png')

            fig, ax = plt.subplots(figsize=(10, 4.5))
            found_any = False

            for mcp in MCP_LIST:
                fname = f'THERM_{mcp}_ep{ep}mm_hx{hx_raw}_dT0{dt}.csv'
                fpath = os.path.join(script_dir, fname)
                if os.path.exists(fpath):
                    t, ts = read_ts(fpath)
                    ax.plot(t, ts, label=mcp, color=COLORS[mcp], linewidth=1.5)
                    found_any = True
                else:
                    graphs_missing += 1

            if found_any:
                dt_display = f'+{dt}' if dt not in ('0', '-2') else dt
                ax.set_title(f'Comparaison T_surface selon MCP — ep={ep}mm, hx={hx_label}, dT0={dt_display}°C')
                ax.set_xlabel('Temps (heures)')
                ax.set_ylabel('T_surface (°C)')
                ax.legend(fontsize=9)
                ax.grid(True, alpha=0.3)
                fig.tight_layout()
                fig.savefig(out_file, dpi=80)
                graphs_generated += 1
            plt.close(fig)

print(f"✅ {graphs_generated} graphiques générés ({graphs_missing} fichiers manquants ignorés)")

# ── 5. Extraction data/processed/ ────────────────────────────────────────────
# A) Cas standard : ep15mm, hx1 (hx10), dT0+2
# B) Cas extrême  : ep10mm, hx0.5 (hx5), dT0-2
# C) Cas optimal  : ep30mm, hx2 (hx20), dT0+4  (libre)
cases = [
    ('A_standard',  'ep15mm', 'hx10', '2'),
    ('B_extreme',   'ep10mm', 'hx5',  '-2'),
    ('C_optimal',   'ep30mm', 'hx20', '4'),
]

proc_dir = mkd('data', 'processed')
selected = {}

for label, ep_str, hx_str, dt in cases:
    ep_num = int(ep_str.replace('mm','').replace('ep',''))
    for mcp in MCP_LIST:
        fname = f'THERM_{mcp}_{ep_str}_{hx_str}_dT0{dt}.csv'
        src = os.path.join(script_dir, fname)
        if os.path.exists(src):
            dst = os.path.join(proc_dir, f'{label}_{fname}')
            shutil.copy2(src, dst)
    # retenir un exemple pour le tableau web (Eutectic_0C ou premier trouvé)
    for mcp in ['Eutectic_0C'] + MCP_LIST:
        fname = f'THERM_{mcp}_{ep_str}_{hx_str}_dT0{dt}.csv'
        src = os.path.join(script_dir, fname)
        if os.path.exists(src):
            selected[label] = (fname, src)
            break

print("✅ data/processed/ créé")

# ── 6. Note de classement ────────────────────────────────────────────────────
note = """NOTE DE CLASSEMENT — THERM_*.csv
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
mcp/          → 5 dossiers (36 fichiers chacun)
epaisseur/    → 3 dossiers (60 fichiers chacun)
convection/   → 3 dossiers (60 fichiers chacun)
dt0/          → 4 dossiers (45 fichiers chacun)
graphs/       → 36 graphiques PNG organisés par ep/hx/dT0
data/raw/     → CSV d'origine (non modifiés)
data/processed/ → 3 cas sélectionnés pour la page web

SÉLECTION POUR LA PAGE WEB
---------------------------
A) Cas standard  : ep=15mm, hx=1, dT0=+2
B) Cas extrême   : ep=10mm, hx=0.5, dT0=-2
C) Cas optimal   : ep=30mm, hx=2, dT0=+4

COMPATIBILITÉ FUTURES SIMULATIONS
----------------------------------
Les tableaux HTML sont construits sur les colonnes garanties :
  t(s), Tsurface(°C), Flux(W), E(J), f_liq
Toute colonne supplémentaire dans de futurs CSV sera ignorée
sans casser l'affichage de la page web.
"""

with open(os.path.join(script_dir, 'note_classement.txt'), 'w', encoding='utf-8') as f:
    f.write(note)

print("✅ note_classement.txt créé")
print("\n=== TERMINÉ ===")
print("Tu peux maintenant ouvrir web/index.html dans ton navigateur.")
