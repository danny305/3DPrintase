"""Build uncapped linear peptides from standard one-letter amino-acid codes."""
from .utils import PipelineError


AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")


def peptide_sequence(value):
    """Normalize case, rejecting ambiguous residues and non-sequence syntax."""
    sequence = value.strip().upper()
    if not sequence:
        raise ValueError("Peptide sequence must contain at least one amino acid")
    invalid = sorted(set(sequence) - AMINO_ACIDS)
    if invalid:
        raise ValueError(
            f"Invalid peptide residue(s): {', '.join(repr(c) for c in invalid)}. "
            "Use standard one-letter codes ACDEFGHIKLMNPQRSTVWY in N-to-C order "
            "without spaces, separators, or FASTA headers."
        )
    return sequence


def build_peptide(sequence):
    """Return a peptide graph with L stereochemistry and free N/C termini.

    Coordinates, explicit hydrogens, and optional pH processing are supplied by
    prepare_molecule, just as for other molecular sources.
    """
    from rdkit import Chem

    sequence = peptide_sequence(sequence)
    mol = Chem.MolFromSequence(sequence, sanitize=True, flavor=0)
    if mol is None:
        raise PipelineError(f"RDKit could not build peptide {sequence}")
    return mol, {
        "sequence": sequence,
        "sequence_direction": "N-to-C",
        "residue_count": len(sequence),
        "stereochemistry": "L amino acids (glycine is achiral)",
        "topology": "linear",
        "termini": "uncapped amino and carboxyl termini; protonation controlled by --ph",
        "builder": "RDKit MolFromSequence, flavor=0",
    }
