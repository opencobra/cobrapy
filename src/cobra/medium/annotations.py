"""Provide lists and annotations for compartment names and reactions.

Please send a PR if you want to add something here :)

"""

# A dictionary having keys as reaction types and keys as prefixes of
# reaction IDs that usually indicate that the reaction is *not* a reaction
# of the specified type.
excludes = {
    "demand": ["SN_", "SK_", "sink", "EX_", "exchange"],
    "exchange": [
        "demand",
        "DM_",
        "biosynthesis",
        "transcription",
        "replication",
        "SN_",
        "SK_",
        "sink",
    ],
    "sink": [
        "demand",
        "DM_",
        "biosynthesis",
        "transcription",
        "replication",
        "EX_",
        "exchange",
    ],
}


# A dictionary having SBO term identifiers as values and boundary types
# as keys.
sbo_terms = {
    "demand": "SBO:0000628",
    "exchange": "SBO:0000627",
    "sink": "SBO:0000632",
    "biomass": "SBO:0000629",
    "pseudoreaction": "SBO:0000631",
}


# A dictionary having keys as common compartment abbreviations and values
# as alternative names.
compartment_shortlist = {
    "ce": ["cell envelope"],
    "c": [
        "cytoplasm",
        "cytosol",
        "default",
        "in",
        "intra cellular",
        "intracellular",
        "intracellular region",
        "intracellular space",
    ],
    "er": ["endoplasmic reticulum"],
    "erm": ["endoplasmic reticulum membrane"],
    "e": [
        "extracellular",
        "extraorganism",
        "out",
        "extracellular space",
        "extra organism",
        "extra cellular",
        "extra-organism",
        "external",
        "external medium",
    ],
    "f": ["flagellum", "bacterial-type flagellum"],
    "g": ["golgi", "golgi apparatus"],
    "gm": ["golgi membrane"],
    "h": ["chloroplast"],
    "l": ["lysosome"],
    "im": ["mitochondrial intermembrane space"],
    "mm": ["mitochondrial membrane"],
    "m": ["mitochondrion", "mitochondria"],
    "n": ["nucleus"],
    "p": ["periplasm", "periplasmic space"],
    "x": ["peroxisome", "glyoxysome"],
    "u": ["thylakoid"],
    "vm": ["vacuolar membrane"],
    "v": ["vacuole"],
    "w": ["cell wall"],
    "s": ["eyespot", "eyespot apparatus", "stigma"],
}
