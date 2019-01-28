# -*- coding: utf-8 -*-

"""Lists and annotations for compartment names and reactions.

Please send a PR if you want to add something here :)
"""

excludes = {
    "demand": ["SN_", "SK_", "sink", "EX_", "exchange"],
    "exchange": ["demand", "DM_", "biosynthesis", "transcription",
                 "replication", "SN_", "SK_", "sink"],
    "sink": ["demand", "DM_", "biosynthesis", "transcription",
             "replication", "EX_", "exchange"]
}
"""A list of sub-strings in reaction IDs that usually indicate
that the reaction is *not* a reaction of the specified type."""

sbo_terms = {"demand": "SBO:0000628",
             "exchange": "SBO:0000627",
             "sink": "SBO:0000632",
             "biomass": "SBO:0000629",
             "pseudoreaction": "SBO:0000631"}
"""SBO term identifiers for various boundary types."""

compartment_shortlist = {
    "ce": ["cell envelope"],
    "c": ["cytoplasm", "cytosol", "default", "in", "intra cellular",
          "intracellular", "intracellular region", "intracellular space"],
    "er": ["endoplasmic reticulum"],
    "erm": ["endoplasmic reticulum membrane"],
    "e": ["extracellular", "extraorganism", "out", "extracellular space",
          "extra organism", "extra cellular", "extra-organism", "external",
          "external medium"],
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
    "s": ["eyespot", "eyespot apparatus", "stigma"]
}
"""A list of common compartment abbreviations and alternative names."""
