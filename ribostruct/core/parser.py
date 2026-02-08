import os
from typing import List, Tuple, Optional

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.PDB import PDBParser
import numpy as np


def _parse_gtf_attributes(attr_str: str) -> dict:
    attrs = {}
    for part in attr_str.strip().split(";"):
        part = part.strip()
        if not part:
            continue
        if " " in part:
            key, val = part.split(" ", 1)
            val = val.strip().strip('"')
            attrs[key] = val
        else:
            attrs[part] = ""
    return attrs


def _find_seq_record(key: str, seq_records: dict):
    if key in seq_records:
        return seq_records[key]
    if len(seq_records) == 1:
        # fallback: return the only sequence
        return next(iter(seq_records.values()))
    # try matching against description or id startswith
    for rec in seq_records.values():
        if key == rec.id or key == rec.name:
            return rec
        if key in rec.description:
            return rec
        if rec.id.split()[0] == key:
            return rec
    return None


def parse_genomic_data(fasta_path: str, gtf_path: str, gene_id: Optional[str] = None, transcript_id: Optional[str] = None) -> Tuple[str, List[Tuple[str, int]]]:
    """
    Parse FASTA and GTF to extract the concatenated CDS nucleotide sequence and a coordinate map.

    Args:
        fasta_path: Path to FASTA file containing genomic sequences (chromosomes/contigs).
        gtf_path: Path to GTF/GFF file containing CDS features.
        gene_id: Optional gene_id to select. If provided, only transcripts belonging to this gene are considered.
        transcript_id: Optional transcript_id to select. If provided, that transcript's CDS is used.

    Returns:
        tuple: (nucleotide_sequence: str, coordinate_map: list of (seqname, genomic_pos) for each nucleotide in transcript order)

    Raises:
        FileNotFoundError: if either input file is missing.
        ValueError: if no CDS features are found or fasta sequence cannot be matched.
    """
    if not os.path.exists(fasta_path):
        raise FileNotFoundError(f"FASTA file not found: {fasta_path}")
    if not os.path.exists(gtf_path):
        raise FileNotFoundError(f"GTF file not found: {gtf_path}")

    seq_records = {rec.id: rec for rec in SeqIO.parse(fasta_path, "fasta")}
    if not seq_records:
        raise ValueError(f"No sequences found in FASTA: {fasta_path}")

    # Parse GTF and collect CDS features
    cds_features = []
    with open(gtf_path, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 9:
                continue
            seqname, source, feature, start, end, score, strand, frame, attrs = parts[:9]
            if feature != "CDS":
                continue
            attrd = _parse_gtf_attributes(attrs)
            cds_features.append({
                "seqname": seqname,
                "start": int(start),
                "end": int(end),
                "strand": strand,
                "attributes": attrd,
            })

    if not cds_features:
        raise ValueError(f"No CDS features found in GTF: {gtf_path}")

    # Group by transcript_id
    transcripts = {}
    for f in cds_features:
        tid = f["attributes"].get("transcript_id") or f["attributes"].get("transcript") or ""
        gid = f["attributes"].get("gene_id") or f["attributes"].get("gene") or ""
        if gene_id and gid and gid != gene_id:
            continue
        transcripts.setdefault(tid, {"records": [], "gene_id": gid, "seqname": f["seqname"], "strand": f["strand"]})
        transcripts[tid]["records"].append((f["start"], f["end"]))

    if transcript_id:
        if transcript_id not in transcripts:
            raise ValueError(f"transcript_id '{transcript_id}' not found in GTF (after filtering by gene_id if provided).")
        chosen_tid = transcript_id
    else:
        # Choose the transcript with the longest total CDS length
        best = None
        best_len = -1
        for tid, v in transcripts.items():
            total = sum(e - s + 1 for s, e in v["records"])
            if total > best_len:
                best = tid
                best_len = total
        if best is None:
            raise ValueError("No transcript found after filtering GTF (maybe gene_id filter removed all transcripts).")
        chosen_tid = best

    chosen = transcripts[chosen_tid]
    seqname = chosen["seqname"]
    strand = chosen["strand"]
    seq_rec = _find_seq_record(seqname, seq_records)
    if seq_rec is None:
        raise ValueError(f"Sequence '{seqname}' referenced in GTF not found in FASTA.")

    # Normalize and sort intervals in transcript order
    intervals = chosen["records"]
    if strand == "+":
        intervals = sorted(intervals, key=lambda x: x[0])
    else:
        intervals = sorted(intervals, key=lambda x: x[0], reverse=True)

    seq_parts: List[str] = []
    coordinate_map: List[Tuple[str, int]] = []

    for s, e in intervals:
        if s < 1 or e > len(seq_rec.seq):
            # allow but warn — here we raise to keep behavior strict
            raise ValueError(f"CDS coordinates {s}-{e} are out of bounds for sequence '{seqname}' (length {len(seq_rec.seq)}).")
        seg = seq_rec.seq[s - 1 : e]
        if strand == "+":
            seq_parts.append(str(seg))
            coordinate_map.extend([(seqname, pos) for pos in range(s, e + 1)])
        else:
            # for negative strand, take reverse-complement of the exon and map coordinates from end->start
            seq_parts.append(str(seg.reverse_complement()))
            coordinate_map.extend([(seqname, pos) for pos in range(e, s - 1, -1)])

    nucleotide_seq = "".join(seq_parts)

    if len(nucleotide_seq) != len(coordinate_map):
        raise RuntimeError("Internal error: nucleotide sequence length disagrees with coordinate map length.")

    return nucleotide_seq, coordinate_map


def parse_pdb_structure(pdb_file: str, chain_id: Optional[str] = None) -> Tuple[str, List]:
    """
    Parse a PDB file to extract the amino acid sequence from ATOM records.
    
    This function extracts the actual amino acid sequence from the 3D structure
    coordinates (ATOM records), not from the SEQRES header. This is critical
    because SEQRES may contain residues that are not physically present in the
    structure (e.g., disordered/missing regions).
    
    Args:
        pdb_file: Path to the PDB file
        chain_id: Optional chain identifier (e.g., 'A', 'B'). If None, uses the first chain.
    
    Returns:
        tuple: (sequence_string, residue_list)
            - sequence_string: One-letter amino acid sequence (e.g., "MKTIIAL...")
            - residue_list: List of Bio.PDB.Residue objects in sequence order
    
    Raises:
        FileNotFoundError: If PDB file doesn't exist
        ValueError: If no residues are found or specified chain doesn't exist
    
    Example:
        >>> seq, residues = parse_pdb_structure("1abc.pdb", chain_id="A")
        >>> print(seq)
        'MKTIIALSY'
        >>> print(len(residues))
        9
    """
    if not os.path.exists(pdb_file):
        raise FileNotFoundError(f"PDB file not found: {pdb_file}")
    
    # Standard three-letter to one-letter amino acid code
    aa_dict = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
        'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
        'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
    }
    
    # Parse the PDB structure
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('protein', pdb_file)
    
    # Collect all available chains
    available_chains = set()
    for model in structure:
        for chain in model:
            available_chains.add(chain.id)
    
    if not available_chains:
        raise ValueError(f"No chains found in PDB file: {pdb_file}")
    
    # Select the chain to use
    if chain_id is None:
        # Use the first chain found
        selected_chain_id = sorted(available_chains)[0]
    else:
        if chain_id not in available_chains:
            raise ValueError(
                f"Chain '{chain_id}' not found in PDB. Available chains: {sorted(available_chains)}"
            )
        selected_chain_id = chain_id
    
    # Extract residues from the selected chain
    sequence = []
    residues = []
    
    for model in structure:
        for chain in model:
            if chain.id == selected_chain_id:
                for residue in chain:
                    # Only consider standard amino acids (hetero flag is blank space)
                    if residue.id[0] == ' ':
                        res_name = residue.get_resname()
                        if res_name in aa_dict:
                            sequence.append(aa_dict[res_name])
                            residues.append(residue)
                # After finding the chain, we can break
                break
        # Typically we only process the first model
        break
    
    if not sequence:
        raise ValueError(f"No standard amino acid residues found in chain '{selected_chain_id}' of PDB file: {pdb_file}")
    
    sequence_string = ''.join(sequence)
    return sequence_string, residues



def parse_ribo_density(density_file: str, chrom: str, start: int, end: int) -> np.ndarray:
    """
    Parse a bedGraph file and extract ribosome density scores for a specific genomic region.
    
    BedGraph format is tab-separated: chrom start end score
    This function reads the file, filters for the specified chromosome and coordinate range,
    and returns a density vector with one score per nucleotide position.
    
    Args:
        density_file: Path to bedGraph format file
        chrom: Chromosome/contig name (e.g., 'chr1', 'mock_chrom')
        start: Start coordinate (1-based, inclusive)
        end: End coordinate (1-based, inclusive)
    
    Returns:
        numpy array of density scores, one per nucleotide position in [start, end]
        Missing positions are filled with 0.0
    
    Raises:
        FileNotFoundError: If density file doesn't exist
        ValueError: If start > end or coordinates are invalid
    
    Example:
        >>> density = parse_ribo_density("sample.bedgraph", "chr1", 100, 150)
        >>> len(density)
        51  # 150 - 100 + 1
    
    Note:
        BedGraph uses 0-based half-open coordinates [start, end), but this function
        expects and returns 1-based inclusive coordinates for consistency with GTF.
    """
    if not os.path.exists(density_file):
        raise FileNotFoundError(f"Density file not found: {density_file}")
    
    if start > end:
        raise ValueError(f"Invalid coordinates: start ({start}) > end ({end})")
    
    if start < 1:
        raise ValueError(f"Start coordinate must be >= 1, got {start}")
    
    # Initialize density vector with zeros
    length = end - start + 1
    density_vector = np.zeros(length, dtype=float)
    
    # Read bedGraph file and populate density vector
    with open(density_file, 'r') as fh:
        for line in fh:
            line = line.strip()
            
            # Skip empty lines and comments
            if not line or line.startswith('#') or line.startswith('track') or line.startswith('browser'):
                continue
            
            parts = line.split('\t')
            if len(parts) < 4:
                continue
            
            bed_chrom = parts[0]
            bed_start = int(parts[1])  # 0-based
            bed_end = int(parts[2])    # 0-based, exclusive
            score = float(parts[3])
            
            # Skip if different chromosome
            if bed_chrom != chrom:
                continue
            
            # Convert bedGraph 0-based to 1-based coordinates
            bed_start_1based = bed_start + 1
            bed_end_1based = bed_end  # bedGraph end is already exclusive, so it becomes inclusive in 1-based
            
            # Calculate overlap with our target region
            overlap_start = max(start, bed_start_1based)
            overlap_end = min(end, bed_end_1based)
            
            # If there's overlap, fill those positions
            if overlap_start <= overlap_end:
                # Convert to indices in our density vector (0-based)
                vector_start_idx = overlap_start - start
                vector_end_idx = overlap_end - start + 1
                
                density_vector[vector_start_idx:vector_end_idx] = score
    
    return density_vector


__all__ = ["parse_genomic_data", "parse_pdb_structure", "parse_ribo_density"]
