#!/usr/bin/env python3
"""
Generate a biologically and geometrically reasonable mock dataset:

1. `mock.pdb`:
   - 15-residue alpha helix (chain A).
   - Backbone atoms N, CA, C, O for every residue.
   - Coordinates computed using standard alpha-helix parameters:
       * Rise      = 1.5 Å per residue
       * Rotation  = 100° per residue
       * Radius    = 2.3 Å (CA positions)
   - B-factors initialized to 0.00 for all atoms.
   - Strict PDB column formatting.

2. `mock.fasta`:
   - Nucleotide coding sequence (CDS) that translates exactly to the
     protein sequence: MKTIIALSYIFCLVF

3. `mock.gtf`:
   - Single CDS entry covering the full coding sequence.

The script prints:

    Generated mock.pdb with [X] atoms.
"""

from __future__ import annotations

import math
from pathlib import Path
from textwrap import wrap


ROOT_DIR = Path(__file__).resolve().parent

PDB_FILENAME = "mock.pdb"
FASTA_FILENAME = "mock.fasta"
GTF_FILENAME = "mock.gtf"

CHROM_NAME = "mock_chrom"
GENE_ID = "GENE1"
TRANSCRIPT_ID = "TX1"

# Fixed protein sequence (15 residues, as requested)
PROTEIN_SEQUENCE = "MKTIIALSYIFCLVF"

# One simple coding DNA triplet per amino acid
CODON_TABLE = {
    "M": "ATG",
    "K": "AAA",
    "T": "ACC",
    "I": "ATA",
    "A": "GCT",
    "L": "CTG",
    "S": "TCT",
    "Y": "TAC",
    "F": "TTC",
    "C": "TGC",
    "V": "GTG",
}

AA3 = {
    "A": "ALA",
    "R": "ARG",
    "N": "ASN",
    "D": "ASP",
    "C": "CYS",
    "Q": "GLN",
    "E": "GLU",
    "G": "GLY",
    "H": "HIS",
    "I": "ILE",
    "L": "LEU",
    "K": "LYS",
    "M": "MET",
    "F": "PHE",
    "P": "PRO",
    "S": "SER",
    "T": "THR",
    "W": "TRP",
    "Y": "TYR",
    "V": "VAL",
}


def protein_to_cds(protein_seq: str) -> str:
    """Return a coding DNA sequence that translates exactly to `protein_seq`."""
    codons = []
    for aa in protein_seq:
        if aa not in CODON_TABLE:
            raise ValueError(f"No codon defined for amino acid '{aa}'")
        codons.append(CODON_TABLE[aa])
    return "".join(codons)


def write_fasta(path: Path, dna_seq: str) -> None:
    """Write the coding DNA sequence as a single FASTA record."""
    lines = [">mock_cds"]
    lines.extend(wrap(dna_seq, 60))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_gtf(path: Path, dna_seq: str) -> None:
    """
    Write a minimal GTF file with a single CDS entry that spans the
    full length of the nucleotide sequence.
    """
    length = len(dna_seq)
    attributes = f'gene_id "{GENE_ID}"; transcript_id "{TRANSCRIPT_ID}";'
    line = (
        f"{CHROM_NAME}\tmock\tCDS\t1\t{length}\t.\t+\t0\t"
        f"{attributes}"
    )
    path.write_text(line + "\n", encoding="utf-8")


def _format_atom_line(
    serial: int,
    atom_name: str,
    resname: str,
    chain_id: str,
    resseq: int,
    x: float,
    y: float,
    z: float,
    element: str,
) -> str:
    """
    Format an ATOM record according to PDB fixed-width columns.

    Columns (1-based):
      1-6   Record name   "ATOM  "
      7-11  Serial        (integer, right-justified)
      13-16 Atom name     (right/center-justified)
      17    AltLoc        (blank)
      18-20 Residue name
      22    Chain ID
      23-26 Residue sequence number
      31-38 X coordinate  (8.3)
      39-46 Y coordinate  (8.3)
      47-54 Z coordinate  (8.3)
      55-60 Occupancy     (6.2)
      61-66 B-factor      (6.2)
      77-78 Element
    """
    occupancy = 1.00
    bfactor = 0.00
    return (
        f"ATOM  "
        f"{serial:5d} "
        f"{atom_name:^4}"
        f" "
        f"{resname:>3} "
        f"{chain_id}"
        f"{resseq:4d}"
        f"    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}"
        f"{occupancy:6.2f}{bfactor:6.2f}"
        f"          "
        f"{element:>2}\n"
    )


def _build_alpha_helix_backbone(protein_seq: str):
    """
    Build backbone atom coordinates (N, CA, C, O) for a simple alpha helix.

    - CA atoms lie on an ideal helix with:
        radius  = 2.3 Å
        rise    = 1.5 Å / residue
        rotate  = 100° / residue
    - N and C are placed approximately along the helix tangent, at standard
      N-CA and CA-C distances.
    - O is placed at a reasonable offset from C (1.24 Å bond length).

    Geometry is approximate but chemically sensible and sufficient for
    structural tooling like Biopython to parse and visualize.
    """
    atoms = []
    serial = 1
    chain_id = "A"

    radius = 2.3
    rise_per_res = 1.5
    angle_per_res = math.radians(100.0)

    n_ca_dist = 1.46
    ca_c_dist = 1.53
    c_o_dist = 1.24

    for i, aa in enumerate(protein_seq):
        res_index = i + 1
        resname = AA3[aa]

        theta = i * angle_per_res
        z = i * rise_per_res

        # CA lies on the helix.
        ca_x = radius * math.cos(theta)
        ca_y = radius * math.sin(theta)
        ca_z = z

        # Tangent direction along the helix (approximate derivative).
        # For x = r cos(theta), y = r sin(theta):
        #   dx/dtheta = -r sin(theta), dy/dtheta = r cos(theta)
        t_x = -math.sin(theta)
        t_y = math.cos(theta)
        t_z = 0.0
        t_norm = math.sqrt(t_x * t_x + t_y * t_y + t_z * t_z) or 1.0
        t_x /= t_norm
        t_y /= t_norm
        t_z /= t_norm

        # Axis direction (z-axis).
        axis_x, axis_y, axis_z = 0.0, 0.0, 1.0

        # N slightly "before" CA along the helix, with a small shift along the axis.
        n_x = ca_x - n_ca_dist * t_x
        n_y = ca_y - n_ca_dist * t_y
        n_z = ca_z - 0.2

        # C slightly "after" CA along the helix, with a small shift along the axis.
        c_x = ca_x + ca_c_dist * t_x
        c_y = ca_y + ca_c_dist * t_y
        c_z = ca_z + 0.2

        # O roughly positioned off C, pointing toward the helix axis and upward.
        v_rad_x = -c_x  # from C toward axis (0,0,*)
        v_rad_y = -c_y
        v_rad_z = 0.0
        o_dir_x = v_rad_x + axis_x
        o_dir_y = v_rad_y + axis_y
        o_dir_z = v_rad_z + axis_z
        o_norm = math.sqrt(o_dir_x * o_dir_x + o_dir_y * o_dir_y + o_dir_z * o_dir_z) or 1.0
        o_dir_x /= o_norm
        o_dir_y /= o_norm
        o_dir_z /= o_norm

        o_x = c_x + c_o_dist * o_dir_x
        o_y = c_y + c_o_dist * o_dir_y
        o_z = c_z + c_o_dist * o_dir_z

        # Append backbone atoms in standard order: N, CA, C, O.
        atoms.append(
            {
                "serial": serial,
                "name": "N",
                "resname": resname,
                "chain": chain_id,
                "resseq": res_index,
                "x": n_x,
                "y": n_y,
                "z": n_z,
                "element": "N",
            }
        )
        serial += 1

        atoms.append(
            {
                "serial": serial,
                "name": "CA",
                "resname": resname,
                "chain": chain_id,
                "resseq": res_index,
                "x": ca_x,
                "y": ca_y,
                "z": ca_z,
                "element": "C",
            }
        )
        serial += 1

        atoms.append(
            {
                "serial": serial,
                "name": "C",
                "resname": resname,
                "chain": chain_id,
                "resseq": res_index,
                "x": c_x,
                "y": c_y,
                "z": c_z,
                "element": "C",
            }
        )
        serial += 1

        atoms.append(
            {
                "serial": serial,
                "name": "O",
                "resname": resname,
                "chain": chain_id,
                "resseq": res_index,
                "x": o_x,
                "y": o_y,
                "z": o_z,
                "element": "O",
            }
        )
        serial += 1

    return atoms


def write_pdb(path: Path, protein_seq: str) -> int:
    """Write the alpha-helix PDB and return the number of atoms written."""
    atoms = _build_alpha_helix_backbone(protein_seq)

    lines = []
    lines.append("HEADER    MOCK ALPHA HELIX FOR TESTING\n")
    lines.append(
        "TITLE     MOCK 15-RESIDUE ALPHA HELIX GENERATED BY generate_mock_data.py\n"
    )

    if protein_seq:
        first_resname = AA3[protein_seq[0]]
        last_resname = AA3[protein_seq[-1]]
        length = len(protein_seq)
        # Simple HELIX record describing the full chain as an alpha helix.
        lines.append(
            f"HELIX    1  H1 {first_resname:>3} A   1 "
            f"{last_resname:>3} A{length:4d}  1"
            f"                                  {length:2d}\n"
        )

    for atom in atoms:
        lines.append(
            _format_atom_line(
                serial=atom["serial"],
                atom_name=atom["name"],
                resname=atom["resname"],
                chain_id=atom["chain"],
                resseq=atom["resseq"],
                x=atom["x"],
                y=atom["y"],
                z=atom["z"],
                element=atom["element"],
            )
        )

    if atoms:
        last = atoms[-1]
        lines.append(
            f"TER   {last['serial'] + 1:5d}      "
            f"{last['resname']:>3} {last['chain']}{last['resseq']:4d}\n"
        )

    lines.append("END\n")

    path.write_text("".join(lines), encoding="utf-8")
    return len(atoms)


def main() -> None:
    protein_seq = PROTEIN_SEQUENCE
    dna_seq = protein_to_cds(protein_seq)

    fasta_path = ROOT_DIR / FASTA_FILENAME
    gtf_path = ROOT_DIR / GTF_FILENAME
    pdb_path = ROOT_DIR / PDB_FILENAME

    write_fasta(fasta_path, dna_seq)
    write_gtf(gtf_path, dna_seq)
    atom_count = write_pdb(pdb_path, protein_seq)

    print(f"Generated mock.pdb with {atom_count} atoms.")


if __name__ == "__main__":
    main()

