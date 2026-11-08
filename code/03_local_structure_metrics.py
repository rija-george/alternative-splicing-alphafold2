"""
Step 3 - Local structural comparison around splice sites (run inside PyMOL
or with `pymol -cq`). Structures (AlphaFold2 rank_001 models / PDB entries)
must already be loaded as PyMOL objects.

Code excerpt as published in Chapter 4.3 of the thesis. `get_sequence_length`
was not printed in the thesis; a minimal version is marked TODO.
"""
import numpy as np
from pymol import cmd


def calculate_local_rmsd(isoform_A, isoform_B, range_A, range_B):
    """Local RMSD between two residue ranges (e.g. "475-504" vs "1-30")."""
    cmd.super(f"{isoform_A} and resi {range_A} and name CA",
              f"{isoform_B} and resi {range_B} and name CA")
    cmd.select("region_A", f"{isoform_A} and resi {range_A} and name CA+CB+N+C+O")
    cmd.select("region_B", f"{isoform_B} and resi {range_B} and name CA+CB+N+C+O")
    model_A = cmd.get_model("region_A")
    model_B = cmd.get_model("region_B")

    rmsd_sum = 0.0
    num_atoms = len(model_A.atom)
    for atom_A, atom_B in zip(model_A.atom, model_B.atom):
        dx = atom_A.coord[0] - atom_B.coord[0]
        dy = atom_A.coord[1] - atom_B.coord[1]
        dz = atom_A.coord[2] - atom_B.coord[2]
        rmsd_sum += dx**2 + dy**2 + dz**2
    return (rmsd_sum / num_atoms) ** 0.5


def get_sequence_length(structure, chain_id="A"):
    """TODO: original implementation not printed in the thesis."""
    return cmd.count_atoms(f"model {structure} and chain {chain_id} and name CA")


def calculate_local_gdt_ts(reference_structure, model_structure,
                           splice_site_range, region_size=10, chain_id="A"):
    """Local GDT_TS (%) around a splice site, CA atoms, cut-offs 1/2/4/8 A."""
    splice_start, splice_end = splice_site_range
    ref_length = get_sequence_length(reference_structure, chain_id)
    model_length = get_sequence_length(model_structure, chain_id)

    if splice_start == 1:                                   # N-terminal event
        selection_start, selection_end = splice_end, splice_end + region_size
    elif splice_end in (ref_length, model_length):          # C-terminal event
        selection_start, selection_end = splice_start - region_size, splice_start
    else:                                                   # internal event
        selection_start, selection_end = splice_start - region_size, splice_end + region_size

    ref_sel = (f"(chain {chain_id} and resi {selection_start}-{selection_end}) "
               f"and name CA and model {reference_structure}")
    model_sel = (f"(chain {chain_id} and resi {selection_start}-{selection_end}) "
                 f"and name CA and model {model_structure}")

    cmd.align(model_sel, ref_sel)
    ref_coords = cmd.get_coords(ref_sel)
    model_coords = cmd.get_coords(model_sel)

    distances = np.linalg.norm(np.array(ref_coords) - np.array(model_coords), axis=1)
    counts = [np.sum(distances <= t) for t in (1.0, 2.0, 4.0, 8.0)]
    gdt_ts = np.mean([c / len(distances) for c in counts])
    return gdt_ts * 100
