"""Bibliothèque de matériaux pour la comparaison multi-matériaux.

Les valeurs sont des ordres de grandeur usuels pour une plateforme ferroviaire,
tirés de la littérature courante sur les structures d'assise. Elles servent à
comparer des comportements, pas à dimensionner : un dimensionnement réel
demanderait des essais sur les matériaux effectivement mis en œuvre.

Chaque entrée porte :
    E        module de Young (Pa)
    nu       coefficient de Poisson
    alpha    coefficient de dilatation thermique (1/°C)
    model    "elastic" ou "cap"
    cap      paramètres du modèle à cap, si model == "cap"
    color    couleur d'affichage (cohérente entre les graphiques du front)
"""
from fem_ballast_beton import BALLAST_PARAMS, CrushableCap, LinearElastic

MATERIALS = {
    "ballast": {
        "name": {"en": "Ballast (crushed stone)", "fr": "Ballast (pierre concassée)"},
        "note": {"en": "Freshly tamped granular layer, elasto-plastic cap model.",
                 "fr": "Couche granulaire fraîchement bourrée, modèle élasto-plastique à cap."},
        "E": 150e6, "nu": 0.25, "alpha": 1.2e-5, "model": "cap",
        "cap": dict(BALLAST_PARAMS), "color": "#f59e0b",
    },
    "ballast_consolide": {
        "name": {"en": "Consolidated ballast", "fr": "Ballast consolidé"},
        "note": {"en": "Ballast after service compaction: stiffer, higher preconsolidation.",
                 "fr": "Ballast après compaction en service : plus raide, préconsolidation plus élevée."},
        "E": 250e6, "nu": 0.25, "alpha": 1.2e-5, "model": "cap",
        "cap": dict(BALLAST_PARAMS, E0=250e6, pc0=150e3, H_cap=6e7), "color": "#fb923c",
    },
    "beton": {
        "name": {"en": "Concrete slab", "fr": "Dalle béton"},
        "note": {"en": "Plain structural concrete, assumed linear elastic.",
                 "fr": "Béton de structure, supposé linéaire élastique."},
        "E": 30e9, "nu": 0.20, "alpha": 1.0e-5, "model": "elastic", "color": "#34d399",
    },
    "grave_bitume": {
        "name": {"en": "Bituminous base course", "fr": "Grave-bitume"},
        "note": {"en": "Asphalt sub-ballast layer; high thermal expansion.",
                 "fr": "Sous-couche bitumineuse ; forte dilatation thermique."},
        "E": 9e9, "nu": 0.35, "alpha": 2.5e-5, "model": "elastic", "color": "#a78bfa",
    },
    "beton_bitumineux": {
        "name": {"en": "Asphalt concrete", "fr": "Béton bitumineux"},
        "note": {"en": "Surface asphalt mix, markedly softer than bituminous base.",
                 "fr": "Enrobé de surface, nettement plus souple que la grave-bitume."},
        "E": 5e9, "nu": 0.35, "alpha": 2.5e-5, "model": "elastic", "color": "#c084fc",
    },
    "grave_non_traitee": {
        "name": {"en": "Unbound granular sub-base", "fr": "Grave non traitée"},
        "note": {"en": "Unbound aggregate sub-base, modelled as elastic here.",
                 "fr": "Sous-couche granulaire non traitée, modélisée ici en élastique."},
        "E": 200e6, "nu": 0.30, "alpha": 1.0e-5, "model": "elastic", "color": "#60a5fa",
    },
    "sol_support": {
        "name": {"en": "Subgrade soil (silt)", "fr": "Sol support (limon)"},
        "note": {"en": "Soft natural subgrade — the worst case for track settlement.",
                 "fr": "Sol naturel peu porteur — le cas défavorable du tassement de voie."},
        "E": 50e6, "nu": 0.35, "alpha": 1.0e-5, "model": "elastic", "color": "#f87171",
    },
    "acier": {
        "name": {"en": "Steel", "fr": "Acier"},
        "note": {"en": "Reference upper bound: a metallic sleeper is three orders "
                       "of magnitude stiffer than ballast.",
                 "fr": "Borne supérieure de référence : une traverse métallique est "
                       "trois ordres de grandeur plus raide que le ballast."},
        "E": 210e9, "nu": 0.30, "alpha": 1.2e-5, "model": "elastic", "color": "#94a3b8",
    },
}

DEFAUT = ["ballast", "beton", "grave_bitume", "sol_support"]


def constitutive(key, E_override=None):
    """Instancie la loi de comportement d'un matériau de la bibliothèque."""
    m = MATERIALS[key]
    E = E_override if E_override else m["E"]
    if m["model"] == "cap":
        return CrushableCap(**dict(m["cap"], E0=E)), True
    return LinearElastic(E, m["nu"]), False


def catalogue():
    """Description sérialisable de la bibliothèque, pour le front.

    Les libellés partent dans les deux langues : le front bascule entre elles
    sans redemander le catalogue au serveur.
    """
    return [
        {
            "id": key,
            "name": m["name"],
            "note": m["note"],
            "E_MPa": m["E"] / 1e6,
            "nu": m["nu"],
            "alpha": m["alpha"],
            "model": m["model"],
            "color": m["color"],
        }
        for key, m in MATERIALS.items()
    ]
