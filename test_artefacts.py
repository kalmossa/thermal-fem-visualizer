"""Tests des artefacts embarqués dans le dépôt.

La page HTML contient un jeu de données de repli, affiché quand elle est ouverte
sans serveur. Rien n'empêche mécaniquement quelqu'un de modifier le solveur sans
relancer `tools/make_offline_data.py` : la démo hors-ligne afficherait alors
silencieusement les résultats d'une version antérieure, exactement le genre
d'écart que ce projet a déjà payé une fois.

Ces tests valent aussi pour les fichiers CSV livrés et pour la cohérence des
dictionnaires de traduction.

    py -m pytest test_artefacts.py -v
"""
import io
import json
import os
import re

import numpy as np
import pytest

RACINE = os.path.dirname(os.path.abspath(__file__))
PAGE = os.path.join(RACINE, "fem_visualizer.html")


@pytest.fixture(scope="module")
def html():
    return io.open(PAGE, encoding="utf-8").read()


@pytest.fixture(scope="module")
def offline(html):
    debut = html.index("const OFFLINE = ") + len("const OFFLINE = ")
    fin = html.index(";\n", debut)
    return json.loads(html[debut:fin])


# ─────────────────────────────────────────────────────────────────────────────
# Jeu de repli hors-ligne
# ─────────────────────────────────────────────────────────────────────────────
def test_jeu_hors_ligne_a_jour(offline):
    """Le jeu embarqué doit correspondre à ce que rend le solveur aujourd'hui.

    En cas d'échec : relancer `py tools/make_offline_data.py` et commiter.
    """
    from fem_ballast_beton import compare

    p = offline["params"]
    r = compare(p_service=p["p"] * 1e3, delta_T=p["delta_T"], sleeper_span=p["span"],
                E_ballast=p["E_b"] * 1e6, alpha_ballast=p["alpha_b"],
                E_concrete=p["E_c"] * 1e6, alpha_concrete=p["alpha_c"],
                lx=p["lx"], ly=p["ly"], nx=p["nx"], ny=p["ny"], nsteps=p["nsteps"])

    assert np.allclose(offline["uy_ballast"], r["uy_ballast"], atol=1e-3), (
        "le jeu hors-ligne a divergé du solveur — relance tools/make_offline_data.py")
    assert np.allclose(offline["uy_beton"], r["uy_beton"], atol=1e-3)
    assert offline["summary"]["bowl_ballast"] == pytest.approx(r["bowl_ballast"], rel=1e-6)


def test_jeu_hors_ligne_complet(offline):
    n_noeuds = offline["mesh"]["n_nodes"]
    n_elem = offline["mesh"]["n_elements"]
    assert len(offline["x_cm"]) == offline["mesh"]["nx"] + 1
    assert len(offline["uy_ballast"]) == len(offline["x_cm"])
    for cle in ("vm_ballast", "vm_beton", "evp_ballast"):
        assert len(offline["fields"][cle]) == n_elem
    for cle in ("uy_nodal_ballast", "uy_nodal_beton"):
        assert len(offline["fields"][cle]) == n_noeuds


def test_jeu_hors_ligne_reste_raisonnable(offline):
    """Garde-fou de taille : la page doit rester ouvrable par double-clic."""
    assert os.path.getsize(PAGE) < 1_000_000, "la page dépasse 1 Mo"


# ─────────────────────────────────────────────────────────────────────────────
# CSV livrés
# ─────────────────────────────────────────────────────────────────────────────
@pytest.mark.parametrize("nom,colonne", [
    ("U_TOP_BALLAST.csv", "uy_ballast"),
    ("U_TOP_BETON.csv", "uy_beton"),
])
def test_csv_livres_a_jour(nom, colonne):
    """Les CSV du dépôt sont les livrables chiffrés cités par le rapport.

    En cas d'échec : relancer `py fem_ballast_beton.py` et commiter.
    """
    from fem_ballast_beton import compare

    donnees = np.loadtxt(os.path.join(RACINE, nom), delimiter=",", skiprows=1)
    r = compare()
    assert np.allclose(donnees[:, 0], r["x_m"], atol=1e-9)
    assert np.allclose(donnees[:, 1], r[colonne], rtol=1e-6), (
        f"{nom} a divergé du solveur — relance py fem_ballast_beton.py")


# ─────────────────────────────────────────────────────────────────────────────
# Traductions
# ─────────────────────────────────────────────────────────────────────────────
def dictionnaires(html):
    debut = html.index("const I18N = {")
    fin = html.index("\n};", debut)
    bloc = html[debut:fin]
    en = set(re.findall(r'^    "([a-z0-9_]+)":', bloc[:bloc.index('  fr: {')], re.M))
    fr = set(re.findall(r'^    "([a-z0-9_]+)":', bloc[bloc.index('  fr: {'):], re.M))
    return en, fr


def test_les_deux_langues_ont_les_memes_cles(html):
    en, fr = dictionnaires(html)
    assert en, "dictionnaire anglais introuvable"
    assert en == fr, (
        f"absentes du français : {sorted(en - fr)} | "
        f"absentes de l'anglais : {sorted(fr - en)}")


def test_pas_de_cle_de_traduction_manquante(html):
    """Toute clé utilisée par le markup ou par t() doit exister au dictionnaire."""
    en, _ = dictionnaires(html)
    utilisees = set(re.findall(r'data-i18n(?:-title)?="([a-z0-9_]+)"', html))
    # `t('cle')` — on écarte les faux positifs du motif, comme getContext('2d')
    utilisees |= {k for k in re.findall(r"[^A-Za-z]t\('([a-z0-9_]+)'", html)
                  if len(k) > 2 and not k[0].isdigit()}
    assert utilisees <= en, f"clés sans traduction : {sorted(utilisees - en)}"


def test_pas_de_francais_en_dur_dans_le_markup(html):
    """L'anglais est la langue par défaut : tout texte visible passe par i18n."""
    debut, fin = html.index("<body"), html.index("<script>", html.index("<body"))
    suspects = [m.group(1).strip() for m in re.finditer(r'>([^<>]+)<', html[debut:fin])
                if re.search(r'[àâäéèêëîïôöùûüç]|\b(les|des|une|dans|pour|avec)\b',
                             m.group(1), re.I)]
    assert not suspects, f"texte français dans le markup : {suspects}"
