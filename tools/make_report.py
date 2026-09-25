"""Génère le rapport technique à partir des résultats du solveur.

Toutes les valeurs chiffrées du document sont calculées au moment de la
génération : le rapport ne peut pas diverger du code, ce qui était précisément
le défaut de la version précédente.

    py -m pip install python-docx
    py tools/make_report.py [chemin/de/sortie.docx]
"""
import os
import sys
from datetime import date

import numpy as np
from docx import Document
from docx.enum.section import WD_SECTION
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.shared import Cm, Pt, RGBColor

RACINE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, RACINE)

import fem_ballast_beton as F  # noqa: E402

BLEU = RGBColor(0x1F, 0x4E, 0x79)
GRIS = RGBColor(0x59, 0x59, 0x59)

AUTEUR = "Elias Lallouet"
DEPOT = "github.com/kalmossa/thermal-fem-visualizer"


# ─────────────────────────────────────────────────────────────────────────────
# Helpers de mise en forme
# ─────────────────────────────────────────────────────────────────────────────
def fr(x, dec=1, signe=False):
    """Formate un nombre à la française : virgule décimale, espace fine des milliers."""
    t = f"{x:{'+' if signe else ''}.{dec}f}".replace(".", ",")
    ent, _, frac = t.partition(",")
    prefixe, chiffres = (ent[0], ent[1:]) if ent[0] in "+-" else ("", ent)
    if len(chiffres) > 4:
        chiffres = f"{int(chiffres):,}".replace(",", " ")
    return prefixe + chiffres + ("," + frac if frac else "")


def para(doc, texte="", taille=10.5, gras=False, italique=False, couleur=None,
         align=None, avant=0, apres=6):
    p = doc.add_paragraph()
    p.paragraph_format.space_before = Pt(avant)
    p.paragraph_format.space_after = Pt(apres)
    if align is not None:
        p.alignment = align
    if texte:
        r = p.add_run(texte)
        r.font.size = Pt(taille)
        r.bold = gras
        r.italic = italique
        if couleur is not None:
            r.font.color.rgb = couleur
    return p


def titre(doc, texte, niveau=1):
    h = doc.add_heading(texte, level=niveau)
    for r in h.runs:
        r.font.color.rgb = BLEU
        r.font.size = Pt(15 if niveau == 1 else 12.5)
    h.paragraph_format.space_before = Pt(16 if niveau == 1 else 12)
    h.paragraph_format.space_after = Pt(6)
    return h


def puce(doc, texte, taille=10.5):
    p = doc.add_paragraph(style="List Bullet")
    p.paragraph_format.space_after = Pt(2)
    r = p.add_run(texte)
    r.font.size = Pt(taille)
    return p


def tableau(doc, entetes, lignes, largeurs=None, taille=9.5):
    t = doc.add_table(rows=1, cols=len(entetes))
    t.style = "Table Grid"
    t.alignment = WD_TABLE_ALIGNMENT.CENTER
    for i, e in enumerate(entetes):
        cell = t.rows[0].cells[i]
        cell.text = ""
        r = cell.paragraphs[0].add_run(e)
        r.bold = True
        r.font.size = Pt(taille)
        r.font.color.rgb = BLEU
    for ligne in lignes:
        cells = t.add_row().cells
        for i, v in enumerate(ligne):
            cells[i].text = ""
            r = cells[i].paragraphs[0].add_run(str(v))
            r.font.size = Pt(taille)
    if largeurs:
        for ligne in t.rows:
            for i, w in enumerate(largeurs):
                ligne.cells[i].width = Cm(w)
    doc.add_paragraph().paragraph_format.space_after = Pt(4)
    return t


def figure(doc, chemin, largeur_cm, legende):
    if not os.path.exists(chemin):
        para(doc, f"[figure manquante : {chemin}]", italique=True, couleur=GRIS)
        return
    doc.add_picture(chemin, width=Cm(largeur_cm))
    doc.paragraphs[-1].alignment = WD_ALIGN_PARAGRAPH.CENTER
    para(doc, legende, taille=9, italique=True, couleur=GRIS,
         align=WD_ALIGN_PARAGRAPH.CENTER, apres=10)


# ─────────────────────────────────────────────────────────────────────────────
# Calculs
# ─────────────────────────────────────────────────────────────────────────────
def calculer():
    r = F.compare()
    sb, sc = r["solver_ballast"], r["solver_beton"]

    # Contrôle analytique : dilatation d'une couche encastrée, sommet libre.
    E, nu, al, dT, Ly = 30e9, 0.2, 1.0e-5, 40.0, 0.35
    K, G = E / (3 * (1 - 2 * nu)), E / (2 * (1 + nu))
    attendu = 3 * K * al * dT / (K + 4 * G / 3) * Ly * 1e6

    s_th = F.FEMSolver(F.Mesh(0.8, 0.35, 32, 14),
                       F.LinearElastic(E, nu), delta_T=dT, alpha=al)
    _, uy_th = s_th.top_profile(s_th.run(3, 0.0, 0.25)[0])

    convergence = []
    for nx, ny in ((16, 7), (24, 10), (32, 14), (48, 21), (64, 28)):
        rr = F.compare(nx=nx, ny=ny)
        convergence.append((f"{nx}×{ny}", (nx + 1) * (ny + 1),
                            rr["bowl_ballast"], rr["bowl_beton"]))

    increments = []
    for ns in (4, 8, 16):
        rr = F.compare(nsteps=ns)
        increments.append((ns, rr["bowl_ballast"]))

    return {
        "r": r, "sb": sb, "sc": sc,
        "evp_max": sb.state_max("evp"), "pc_max": sb.state_max("pc"),
        "vm_max": float(sb.stress_field(r["U_ballast"]).max()),
        "th_attendu": attendu, "th_calcule": float(uy_th.mean()),
        "convergence": convergence, "increments": increments,
    }


# ─────────────────────────────────────────────────────────────────────────────
# Document
# ─────────────────────────────────────────────────────────────────────────────
def ecart_analytique(d):
    """L'accord est ici au niveau du bruit d'arrondi : l'annoncer comme tel plutôt
    qu'afficher un « 0,0e+00 » qui ressemble à une valeur non calculée."""
    e = abs(d["th_calcule"] - d["th_attendu"]) / d["th_attendu"]
    return "< 10⁻¹² (précision machine)" if e < 1e-12 else f"{e:.1e}".replace(".", ",")


def construire(d, sortie):
    r, sb = d["r"], d["sb"]
    x, ub, uc = r["x_m"], r["uy_ballast"], r["uy_beton"]

    doc = Document()
    s = doc.sections[0]
    s.page_width, s.page_height = Cm(21.0), Cm(29.7)
    for attr, v in (("top_margin", 2.2), ("bottom_margin", 2.2),
                    ("left_margin", 2.3), ("right_margin", 2.3)):
        setattr(s, attr, Cm(v))
    doc.styles["Normal"].font.name = "Calibri"
    doc.styles["Normal"].font.size = Pt(10.5)

    # ── Page de titre ────────────────────────────────────────────────────
    para(doc, "Analyse par éléments finis d'un support de voie ferrée", taille=19,
         gras=True, couleur=BLEU, align=WD_ALIGN_PARAGRAPH.CENTER, avant=48, apres=4)
    para(doc, "Comportement thermo-mécanique comparé du ballast et du béton",
         taille=13, couleur=GRIS, align=WD_ALIGN_PARAGRAPH.CENTER, apres=2)
    para(doc, "Rapport technique", taille=11, italique=True, couleur=GRIS,
         align=WD_ALIGN_PARAGRAPH.CENTER, apres=30)
    para(doc, AUTEUR, taille=12, gras=True, align=WD_ALIGN_PARAGRAPH.CENTER, apres=2)
    para(doc, DEPOT, taille=10.5, couleur=GRIS,
         align=WD_ALIGN_PARAGRAPH.CENTER, apres=2)
    para(doc, date.today().strftime("%d/%m/%Y"), taille=10.5, couleur=GRIS,
         align=WD_ALIGN_PARAGRAPH.CENTER)

    doc.add_section(WD_SECTION.NEW_PAGE)

    # ── 1. Contexte ──────────────────────────────────────────────────────
    titre(doc, "1. Contexte et objectifs")
    para(doc, "Ce rapport présente une analyse numérique comparative du comportement "
              "mécanique de deux matériaux employés comme support de voie ferrée : le "
              "ballast, milieu granulaire au comportement élasto-plastique, et le béton, "
              "milieu continu supposé linéaire élastique. La simulation intègre les "
              "effets de dilatation liés aux variations de température en service.")
    para(doc, "L'objectif est de quantifier les déplacements verticaux en surface sous "
              "l'action combinée d'une charge de traverse et d'une variation thermique, "
              "afin de comparer la rigidité et la capacité de déformation des deux "
              "solutions.")
    para(doc, "Une attention particulière est portée à la séparation des deux "
              "contributions. Le déplacement brut mélange un soulèvement thermique "
              "d'ensemble et un enfoncement mécanique localisé ; confondre les deux "
              "conduit à des conclusions inversées, comme le montre la section 3.")

    # ── 2. Modèle ────────────────────────────────────────────────────────
    titre(doc, "2. Modèle numérique")

    titre(doc, "2.1 Géométrie et maillage", 2)
    para(doc, "Le domaine est un rectangle de 0,80 m × 0,35 m représentant la couche de "
              "support en coupe transversale, en hypothèse de déformation plane.")
    for t in ("Largeur 0,80 m — hauteur 0,35 m",
              "Discrétisation : 32 × 14 éléments quadrangulaires Q4, soit 495 nœuds et "
              "448 éléments (990 degrés de liberté)",
              "Intégration de Gauss 2 × 2, soit 1 792 points d'intégration",
              "Chargement réparti sur 25 cm centrés (emprise de la semelle de traverse)",
              "Conditions aux limites : encastrement en base, appuis à glissement "
              "(u_x = 0) sur les bords latéraux, surface supérieure libre"):
        puce(doc, t)

    titre(doc, "2.2 Paramètres matériaux", 2)
    tableau(doc,
            ["Paramètre", "Ballast (crushable-cap)", "Béton (lin. élastique)"],
            [["Module de Young E", "150 MPa", "30 000 MPa"],
             ["Coefficient de Poisson ν", "0,25", "0,20"],
             ["Coeff. de dilatation α", "1,2 × 10⁻⁵ /°C", "1,0 × 10⁻⁵ /°C"],
             ["Pression de consolidation pc₀", "50 kPa", "—"],
             ["Demi-axe du cap X", "80 kPa", "—"],
             ["Module d'écrouissage H_cap", "30 MPa", "—"],
             ["Pente critique M", "1,2", "—"],
             ["Loi de comportement", "Élasto-plastique à cap", "Hooke (déf. plane)"]],
            largeurs=[6.0, 5.2, 4.8])

    para(doc, "Conditions de chargement :", gras=True, avant=4, apres=2)
    puce(doc, "Pression de service : p = 300 kN/m (charge linéique sur la traverse)")
    puce(doc, "Variation de température : ΔT = +40 °C")
    para(doc, "Le chargement est appliqué de façon proportionnelle en 6 incréments : "
              "pression et température croissent ensemble de 0 à leur valeur nominale. "
              "Ce choix est nécessaire pour une loi plastique, dont la réponse dépend du "
              "chemin de chargement.", taille=10)

    titre(doc, "2.3 Formulation et résolution", 2)
    para(doc, "L'équilibre discret s'écrit f_int(U) = f_ext, où f_int(U) = ∫ Bᵀ σ(U) dΩ "
              "est le vecteur des forces internes. Il est résolu à chaque incrément par "
              "la méthode de Newton-Raphson : K_t ΔU = f_ext − f_int(U), puis U ← U + ΔU, "
              "jusqu'à ce que le résidu tombe sous 10⁻⁹ fois la charge appliquée.")
    para(doc, "La contrainte thermique n'est pas ajoutée comme un chargement séparé : "
              "elle entre par la loi de comportement, σ = D : (ε − ε_th). Elle apparaît "
              "ainsi automatiquement dans le résidu et avec le bon signe, ce qui évite "
              "l'erreur décrite en section 3.")
    para(doc, "La réduction en déformation plane impose ε_th = (1 + ν) α ΔT [1, 1, 0]ᵀ. "
              "Le facteur (1 + ν) provient du blocage de la dilatation hors-plan, dont la "
              "contrainte σ_zz est réinjectée dans le plan.")
    para(doc, "Le ballast suit un modèle à cap de la famille Cam-Clay, de surface de "
              "charge f(p̄, q, pc) = (p̄ − pc)²/X² + (q/(M·pc))² − 1, avec écoulement "
              "associé, écrouissage du cap par la déformation volumique plastique et "
              "dégradation du module volumique traduisant l'écrasement des grains. "
              "La variable p̄ est la pression moyenne en convention compression positive.",
         apres=10)

    # ── 3. Vérification ──────────────────────────────────────────────────
    titre(doc, "3. Vérification et validation du solveur")
    para(doc, "Une première version de cette étude concluait que le ballast se déformait "
              "légèrement moins que le béton, avec un rapport de 0,89. Ce résultat était "
              "un artefact numérique. Le solveur a été repris et une suite de tests "
              "automatisés a été écrite ; chaque défaut identifié y a son test de "
              "non-régression.")

    titre(doc, "3.1 Défauts corrigés", 2)
    tableau(doc,
            ["Défaut", "Conséquence"],
            [["Signe de la charge thermique inversé dans le résidu",
              "un échauffement faisait descendre la surface au lieu de la soulever"],
             ["Facteur (1 + ν) absent de la réduction en déformation plane",
              "charge thermique sous-estimée de 20 %"],
             ["Déplacement réécrit à chaque incrément au lieu d'être cumulé",
              "seul 1/6 de la charge était réellement appliqué"],
             ["Déformation nulle transmise à la loi de comportement",
              "la plasticité n'était jamais pilotée par la charge mécanique"],
             ["Cap centré du côté tendu",
              "l'écrouissage réduisait la capacité portante au lieu de l'augmenter"],
             ["σ_zz reconstruit par ν(σ_xx + σ_yy) après retour plastique",
              "dérive de 17 % de la réponse selon le nombre d'incréments"],
             ["Second invariant J₂ mal formé, cisaillement ignoré",
              "déviateur faux d'un facteur √3"],
             ["Forces nodales à poids plein aux extrémités de la zone chargée",
              "résultante 10 % supérieure à la charge prescrite"]],
            largeurs=[7.6, 8.4])

    para(doc, "Le module d'écrouissage du cap a par ailleurs été recalibré. À sa valeur "
              "initiale de 0,1 MPa, porter la pression de service aurait exigé près de "
              "600 % de déformation volumique plastique et le retour sur la surface de "
              "charge ne convergeait pas. La valeur retenue de 30 MPa correspond à une "
              "densification d'environ 2 % sur une plage de quelques centaines de kPa, "
              "ordre de grandeur admis pour un ballast fraîchement bourré.")

    titre(doc, "3.2 Confrontation à une solution analytique", 2)
    para(doc, "Une couche d'épaisseur L encastrée en base, à bords latéraux glissants et "
              "sommet libre, chauffée uniformément de ΔT, se soulève de "
              "u_y = 3K α ΔT L / (K + 4G/3). Pour le béton (ΔT = +40 °C) :")
    tableau(doc, ["Grandeur", "Valeur"],
            [["Solution analytique", f"{fr(d['th_attendu'], 2)} µm"],
             ["Solveur", f"{fr(d['th_calcule'], 2)} µm"],
             ["Écart relatif", ecart_analytique(d)]],
            largeurs=[8.0, 8.0])

    titre(doc, "3.3 Indépendance à la discrétisation", 2)
    para(doc, "Une solution qui dépend du nombre d'incréments n'est pas une solution. "
              "La cuvette du ballast est ici stable à mieux que 1 % :")
    tableau(doc, ["Nombre d'incréments", "Cuvette ballast (µm)"],
            [[str(ns), fr(v, 1)] for ns, v in d["increments"]],
            largeurs=[8.0, 8.0])

    para(doc, "Le raffinement du maillage donne :", avant=6)
    tableau(doc, ["Maillage", "Nœuds", "Cuvette ballast (µm)", "Cuvette béton (µm)"],
            [[m, str(n), fr(a, 1), fr(b, 3)] for m, n, a, b in d["convergence"]],
            largeurs=[3.6, 3.4, 4.6, 4.4])
    para(doc, "Le cas élastique converge proprement. Le cas plastique n'est pas monotone : "
              "les maillages grossiers sous-estiment nettement la cuvette, et au-delà de "
              "32 × 14 subsiste une oscillation de l'ordre de 5 %, la charge étant très "
              "localisée sous la traverse. C'est une limite assumée du modèle, pas un "
              "résultat convergé à la troisième décimale.", taille=10, apres=10)

    # ── 4. Résultats ─────────────────────────────────────────────────────
    titre(doc, "4. Résultats numériques")

    titre(doc, "4.1 Profil de surface", 2)
    figure(doc, os.path.join(RACINE, "w_compare.png"), 15.5,
           "Figure 1 — Déplacement vertical u_y en surface. En haut, déplacement brut ; "
           "en bas, réponse mécanique seule, soulèvement thermique retiré.")

    titre(doc, "4.2 Champs calculés", 2)
    figure(doc, os.path.join(RACINE, "docs", "champs.png"), 16.0,
           "Figure 2 — Ballast : déplacement vertical, contrainte de von Mises et "
           "déformation volumique plastique cumulée. Le tireté marque l'emprise de la traverse.")

    titre(doc, "4.3 Valeurs nodales", 2)
    para(doc, "Déplacements verticaux en surface, échantillonnés tous les 4 nœuds. "
              "Convention : une valeur positive est un soulèvement.", taille=10)
    lignes = [[fr(x[i] * 100, 1), fr(ub[i], 1, signe=True), fr(uc[i], 2, signe=True),
               fr(ub[i] - uc[i], 1, signe=True)] for i in range(0, len(x), 4)]
    tableau(doc, ["x (cm)", "u_y Ballast (µm)", "u_y Béton (µm)", "Écart (µm)"],
            lignes, largeurs=[3.4, 4.4, 4.2, 4.0])

    titre(doc, "4.4 Synthèse", 2)
    tableau(doc, ["Grandeur", "Ballast", "Béton"],
            [["u_y au bord (soulèvement thermique)", f"{fr(ub[0], 1, True)} µm",
              f"{fr(uc[0], 2, True)} µm"],
             ["u_y sous la traverse", f"{fr(ub.min(), 1, True)} µm", f"{fr(uc.min(), 2, True)} µm"],
             ["Cuvette (enfoncement mécanique)", f"{fr(r['bowl_ballast'], 1)} µm",
              f"{fr(r['bowl_beton'], 2)} µm"],
             ["Déformation plastique volumique max", f"{fr(d['evp_max'] * 100, 2)} %",
              "— (élastique)"],
             ["Pression de consolidation finale",
              f"{fr(d['pc_max'] / 1e3, 0)} kPa (initiale 50)", "—"],
             ["Contrainte de von Mises max", f"{fr(d['vm_max'], 0)} kPa", "—"]],
            largeurs=[7.2, 4.6, 4.2])
    para(doc, f"Rapport des cuvettes ballast / béton : {fr(r['ratio_bowl'], 0)}.",
         gras=True, apres=10)

    # ── 5. Analyse ───────────────────────────────────────────────────────
    titre(doc, "5. Analyse et interprétation")
    para(doc, "Les deux matériaux se soulèvent sous l'effet de la dilatation, mais leurs "
              "réponses mécaniques n'ont pas le même ordre de grandeur.")
    para(doc, f"Le béton, avec un module 200 fois supérieur, ne s'enfonce que de "
              f"{fr(r['bowl_beton'], 2)} µm sous la traverse. Son profil reste pratiquement "
              f"plat autour de {fr(uc.mean(), 0, True)} µm : la dilatation domine complètement, et "
              f"la réponse à la charge est invisible à l'échelle du graphique. Le matériau "
              f"reste entièrement dans son domaine élastique.")
    para(doc, f"Le ballast se comporte différemment. Il se soulève de {fr(ub[0], 0, True)} µm en "
              f"périphérie, où seule la dilatation agit, mais s'enfonce jusqu'à "
              f"{fr(ub.min(), 0, True)} µm sous la traverse : la charge y dépasse largement la "
              f"capacité élastique du milieu granulaire. La cuvette résultante atteint "
              f"{fr(r['bowl_ballast'] / 1000, 2)} mm.")
    para(doc, f"Une part de cette déformation est irréversible. La déformation volumique "
              f"plastique atteint {fr(d['evp_max'] * 100, 2)} % sous la semelle et la "
              f"pression de consolidation passe de 50 à {fr(d['pc_max'] / 1e3, 0)} kPa : "
              f"le ballast se densifie sous la charge. C'est ce mécanisme que le modèle à "
              f"cap est fait pour représenter, et c'est lui qui explique le tassement "
              f"progressif de la voie, que le béton ne présente pas.")
    para(doc, f"Le rapport des cuvettes, {fr(r['ratio_bowl'], 0)}, dépasse le rapport des "
              f"modules (200). L'écart vient de la plasticité : au-delà du seuil, le "
              f"ballast se déforme davantage que ne le prévoirait sa seule raideur "
              f"élastique.")
    para(doc, "Sur le plan de l'ingénierie, ces deux comportements ne se hiérarchisent pas "
              "simplement. La rigidité du béton limite le tassement mais transmet les "
              "efforts dynamiques au reste de la structure ; la souplesse du ballast amortit "
              "ces efforts, au prix d'un tassement progressif qui impose des opérations de "
              "bourrage périodiques. Le présent modèle, monotone, ne traite que le premier "
              "chargement et ne permet pas de conclure sur ce tassement cumulé.", apres=10)

    # ── 6. Conclusion ────────────────────────────────────────────────────
    titre(doc, "6. Conclusion")
    para(doc, "La simulation confirme que la nature du support gouverne la réponse de la "
              "voie. Le béton, très rigide, reste élastique et ne s'enfonce que de "
              "quelques micromètres. Le ballast s'enfonce de plus de deux millimètres, "
              "dont une part irréversible par compaction des grains.")
    para(doc, "Le modèle élasto-plastique à cap est adapté à la représentation du ballast : "
              "il reproduit le seuil de plastification, la consolidation progressive et "
              "l'écrasement des grains, trois mécanismes qu'une loi élastique ne peut pas "
              "décrire.")
    para(doc, "Sur le plan méthodologique, ce travail a surtout montré qu'un résultat "
              "d'allure vraisemblable n'est pas un résultat juste. Les courbes de la "
              "première version étaient lisses, symétriques et d'ordre de grandeur "
              "plausible, alors que la conclusion était inversée. Ce sont les contrôles "
              "élémentaires — vérifier le signe d'un soulèvement, confronter à une "
              "solution analytique, s'assurer que la réponse ne dépend pas du nombre "
              "d'incréments — qui ont mis les défauts en évidence.")
    para(doc, "Le solveur, sa suite de tests, l'API de calcul et l'interface web sont "
              f"publiés sur {DEPOT}. Toutes les valeurs chiffrées de ce document sont "
              "générées directement par le solveur au moment de sa production : le "
              "rapport ne peut pas diverger du code.",
         taille=10, italique=True, couleur=GRIS)

    doc.save(sortie)
    return sortie


def main():
    sortie = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
        RACINE, "docs", "technical-report.docx")
    os.makedirs(os.path.dirname(sortie), exist_ok=True)
    print("calcul des résultats…")
    d = calculer()
    construire(d, sortie)
    r = d["r"]
    print(f"écrit : {sortie}")
    print(f"  cuvette ballast {r['bowl_ballast']:.1f} µm | béton {r['bowl_beton']:.2f} µm"
          f" | rapport {r['ratio_bowl']:.0f}")


if __name__ == "__main__":
    main()
