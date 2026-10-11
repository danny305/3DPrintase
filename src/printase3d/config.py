"""Config components of the 3DPrintase pipeline."""
from __future__ import annotations


VERSION = "1.3.0"


DEFAULT_OUTPUT = "out/misc"


RCSB_SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"


RCSB_CIF = "https://files.rcsb.org/download/{pdb}.cif"


PUBCHEM_SDF = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/SDF?record_type=3d"


PUBCHEM_JSON = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/Title,SMILES,ConnectivitySMILES,CanonicalSMILES,IsomericSMILES,MolecularFormula,InChI,Charge/JSON"


PUBCHEM_NAME = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/cids/JSON"


PRESETS = {
    "amylose": {"keywords": ["maltopentaose", "maltohexaose", "maltoheptaose", "maltodextrin"], "components": ["GLC", "BGC", "MAL"]},
    "cellulose": {"keywords": ["cellopentaose", "cellohexaose", "cellooligosaccharide"], "components": ["BGC", "CBI"]},
    "glucan": {"keywords": ["glucan oligosaccharide"], "components": ["GLC", "BGC"]},
    "mannan": {"keywords": ["mannooligosaccharide", "mannan"], "components": ["MAN", "BMA"]},
    "xylan": {"keywords": ["xylooligosaccharide", "xylan"], "components": ["XYP", "XYS"]},
    "chitin": {"keywords": ["chitooligosaccharide", "chitin"], "components": ["NAG", "NDG"]},
    "glycan": {"keywords": ["oligosaccharide complex"], "components": []},
}


RENDER_PRESETS = {
    # Retain the original small-spheres proportions as an explicit option.
    "small-spheres": {"spacefill": 45.0, "wireframe": 0.4},
    # Proportions selected against the supplied MakerWorld reference photograph.
    "ball-stick": {"spacefill": 35.0, "wireframe": 0.3},
}


FRUCTAN_TYPES = {
    "inulin": {
        "linkage": "b2-1",
        "linkage_label": "beta-D-Fruf-(2→1)",
        "name": "inulin-oligosaccharide",
        "conformation_profile": "inulin-beta-2-1",
    },
    "levan": {
        "linkage": "b2-6",
        "linkage_label": "beta-D-Fruf-(2→6)",
        "name": "levan-oligosaccharide",
        "conformation_profile": "levan-beta-2-6",
    },
}


GLUCAN_TYPES = {
    "maltose": {
        "linkage": "a1-4",
        "linkage_label": "alpha-D-Glcp-(1→4)",
        "name": "maltose-oligosaccharide",
        "series_name": "maltooligosaccharide",
        "residue": "Glc",
        "conformation_profile": "amylose-alpha-1-4",
    },
    "cellulose": {
        "linkage": "b1-4",
        "linkage_label": "beta-D-Glcp-(1→4)",
        "name": "cellulose-oligosaccharide",
        "series_name": "cellooligosaccharide",
        "residue": "Glc",
        "conformation_profile": "cellulose-beta-1-4",
    },
    "galactan": {
        "linkage": "b1-4",
        "linkage_label": "beta-D-Galp-(1→4)",
        "name": "beta-galactan-oligosaccharide",
        "series_name": "beta-1-4-galactooligosaccharide",
        "residue": "Gal",
        "conformation_profile": "galactan-beta-1-4",
    },
    "mannan": {
        "linkage": "b1-4",
        "linkage_label": "beta-D-Manp-(1→4)",
        "name": "beta-mannan-oligosaccharide",
        "series_name": "beta-1-4-mannooligosaccharide",
        "residue": "Man",
        "conformation_profile": "mannan-beta-1-4",
    },
}


GLYCAN_CONFORMATION_PROFILES = {
    "amylose-alpha-1-4": {
        "phi": 100.0, "psi": 115.0,
        "helix_residues_per_turn": 6.0, "helix_pitch_angstrom": 7.9, "helix_radius_angstrom": 5.0,
        "shape": "left-handed V-amylose-like helix (about six residues per turn)",
        "source": "https://pmc.ncbi.nlm.nih.gov/articles/PMC16317/",
    },
    "cellulose-beta-1-4": {
        "phi": -80.0, "psi": 100.0,
        "shape": "extended cellulose-like 2/1 screw ribbon",
        "source": "https://link.springer.com/article/10.1007/s10570-025-06424-y",
    },
    "galactan-beta-1-4": {
        # Converted from the reported proton convention (about 46, 22 deg)
        # to O5-C1-O4-C4 / C1-O4-C4-C3 for beta-D-Galp stereochemistry.
        "phi": -68.0, "psi": -99.0,
        "helix_residues_per_turn": 6.0, "helix_pitch_angstrom": 14.6, "helix_radius_angstrom": 3.2,
        "shape": "right-handed beta-1,4-galactan solution-like helix",
        "source": "https://pmc.ncbi.nlm.nih.gov/articles/PMC2975222/",
    },
    "mannan-beta-1-4": {
        # Crystal values are phi=-90 and psi(C5)=-149 deg; psi(C3) is the
        # equivalent approximately +91 deg heavy-atom convention used here.
        "phi": -90.0, "psi": 91.0,
        "shape": "extended beta-mannan-like 2/1 screw ribbon",
        "source": "https://users.man.poznan.pl/wojtekr/pdf/wr-lit-19.pdf",
    },
    "inulin-beta-2-1": {
        "phi": 68.0, "psi": 159.0, "omega": -87.0,
        "helix_residues_per_turn": 6.0, "helix_pitch_angstrom": 14.7, "helix_radius_angstrom": 5.0,
        "shape": "idealized sixfold inulin helix",
        "source": "https://pmc.ncbi.nlm.nih.gov/articles/PMC4268547/",
    },
    "levan-beta-2-6": {
        "phi": -70.0, "psi": 170.0, "omega": 50.0,
        "shape": "low-energy beta-2,6-levan repeat conformation",
        "source": "https://pmc.ncbi.nlm.nih.gov/articles/PMC1302875/",
    },
}


WATER = {"HOH", "WAT", "DOD"}


IONS = {"NA", "K", "CL", "BR", "CA", "MG", "MN", "ZN", "FE", "CU", "CO", "NI", "CD", "IOD"}
