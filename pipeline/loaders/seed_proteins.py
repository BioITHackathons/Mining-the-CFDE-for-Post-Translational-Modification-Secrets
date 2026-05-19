"""
Seed data: neurodegenerative disease proteins and their associations.
Keyed by UniProt accession -> gene symbol.
Disease scope from the project plan (DOID ontology).
"""

# UniProt accession -> gene symbol
SEED_PROTEINS = {
    # Alzheimer's disease
    "P05067": "APP",
    "P49768": "PSEN1",
    "P10636": "MAPT",
    "P56817": "BACE1",
    "Q9NZC2": "TREM2",
    # Parkinson's disease
    "P37840": "SNCA",
    "Q5S007": "LRRK2",
    "Q99497": "PARK7",
    "Q9BXM7": "PINK1",
    "P04062": "GBA",
    # ALS
    "P00441": "SOD1",
    "Q13148": "TARDBP",  # TDP-43
    "P35637": "FUS",
    "Q96LT7": "C9orf72",
    # Huntington's disease
    "P42858": "HTT",
    # Frontotemporal dementia (MAPT already included)
    "P28799": "GRN",
}

# Disease -> DOID mapping
SEED_DISEASES = {
    "DOID:10652": {
        "name": "Alzheimer's disease",
        "category": "neurodegenerative",
        "proteins": ["P05067", "P49768", "P10636", "P56817", "Q9NZC2"],
    },
    "DOID:14330": {
        "name": "Parkinson's disease",
        "category": "neurodegenerative",
        "proteins": ["P37840", "Q5S007", "Q99497", "Q9BXM7", "P04062"],
    },
    "DOID:332": {
        "name": "amyotrophic lateral sclerosis",
        "category": "neurodegenerative",
        "proteins": ["P00441", "Q13148", "P35637", "Q96LT7"],
    },
    "DOID:12858": {
        "name": "Huntington's disease",
        "category": "neurodegenerative",
        "proteins": ["P42858"],
    },
    "DOID:9255": {
        "name": "frontotemporal dementia",
        "category": "neurodegenerative",
        "proteins": ["P10636", "P28799", "Q96LT7"],
    },
}
