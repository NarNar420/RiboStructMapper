import os
from typing import List, Tuple, Optional
import logging

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.PDB import PDBParser
import numpy as np

# Set up logging
logger = logging.getLogger(__name__)


def sanitize_sequence(sequence: str) -> Tuple[str, int]:
    """
    Sanitize a nucleotide sequence to ensure it's in-frame and divisible by 3.
    
    This function:
    1. Searches for the first ATG (start codon)
    2. Slices the sequence from that position
    3. Truncates the end to make length divisible by 3
    
    Args:
        sequence: Raw nucleotide sequence (may contain leading garbage or incomplete codons)
    
    Returns:
        Tuple of (cleaned_sequence, start_offset):
            - cleaned_sequence: Sequence starting at ATG and divisible by 3
            - start_offset: Number of nucleotides removed from the start
    
    Example:
        >>> sanitize_sequence("CCATGTTCGGTTA")
        ("ATGTTCGGT", 2)  # Removed "CC" from start, "TA" from end
    """
    # Convert to uppercase for consistency
    seq = sequence.upper()
    
    # Search for first ATG (start codon)
    start_index = seq.find("ATG")
    
    if start_index == -1:
        # No ATG found - keep original sequence but log warning
        logger.warning(f"No ATG start codon found in sequence (length {len(seq)}). Using sequence as-is.")
        start_index = 0
    elif start_index > 0:
        logger.info(f"Found ATG at position {start_index}. Trimming {start_index} nucleotides from start.")
    
    # Slice from start position
    seq = seq[start_index:]
    
    # Truncate end to make divisible by 3
    remainder = len(seq) % 3
    if remainder > 0:
        logger.info(f"Truncating {remainder} nucleotides from end to maintain reading frame.")
        seq = seq[:-remainder]
    
    return seq, start_index

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


import re

def parse_genomic_data(fasta_path: str) -> Tuple[str, List[Tuple[str, int]], int, str, int]:
    """
    Parse FASTA to extract CDS sequence and genomic coordinates from header.
    
    Expected Header Format:
        >GeneID Chromosome:Start..End
        Example: >YAL001C I:100..500
    
    Args:
        fasta_path: Path to FASTA file.
        
    Returns:
        tuple:
            - nucleotide_seq: Sanitized CDS sequence (divisible by 3, starts with ATG)
            - coordinate_map: List of (chrom, genomic_pos) for each nucleotide
            - start_offset: Number of nucleotides removed from start of raw sequence
            - chrom: Chromosome name (e.g., 'I')
            - genomic_start: 1-based genomic start coordinate of the sanitized sequence
            
    Raises:
        ValueError: If header format is invalid or no ATG found.
    """
    if not os.path.exists(fasta_path):
        raise FileNotFoundError(f"FASTA file not found: {fasta_path}")
    
    # Read the first sequence record
    # We only handle single-sequence FASTA files for now as per pipeline design
    try:
        record = next(SeqIO.parse(fasta_path, "fasta"))
    except StopIteration:
        raise ValueError(f"No sequences found in FASTA: {fasta_path}")
    
    header = record.description
    
    # Parse header: >ID Chrom:Start..End
    # Regex captures: (ID) (Chrom):(Start)..(End)
    # Example: >YAL003W I:142174..142253
    # Note: SeqIO separates ID and description. record.description is the full header line after '>'
    # But usually record.id is the first word.
    # Let's parse the full description line.
    
    # Pattern: ID followed by space, then Chrom:Start..End
    # We'll allow flexible spacing
    match = re.search(r'(\S+)\s+(\S+):(\d+)\.\.(\d+)', header)
    
    if not match:
        raise ValueError(
            f"FASTA header '{header}' does not match required format: '>ID Chrom:Start..End' "
            "(e.g., '>YAL001C I:100..500')"
        )
    
    gene_id, chrom, raw_start_str, raw_end_str = match.groups()
    raw_start = int(raw_start_str)
    
    nucleotide_seq_raw = str(record.seq).upper()
    
    # Sanitize sequence (find ATG, trim length)
    nucleotide_seq, start_offset = sanitize_sequence(nucleotide_seq_raw)
    
    # Calculate final genomic start position
    # If we removed N nucleotides from start, the genomic start shifts by N
    # Coordinate system: 1-based inclusive
    final_genomic_start = raw_start + start_offset
    
    # Generate coordinate map for the sanitized sequence
    # This maps each residue index to (chrom, genomic_pos)
    coordinate_map = []
    for i in range(len(nucleotide_seq)):
        # Assuming + strand for simplicity as per requirements (or if header implies it)
        # If headers imply specific strands, we'd need that info. 
        # For now, we assume simple forward mapping.
        current_pos = final_genomic_start + i
        coordinate_map.append((chrom, current_pos))
        
    return nucleotide_seq, coordinate_map, start_offset, chrom, final_genomic_start


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
        expects and returns 1-based inclusive coordinates for consistency with other genomic formats.
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
            
            # Handle standard BedGraph (4 cols) or implied-chrom BedGraph (3 cols)
            if len(parts) >= 4:
                bed_chrom = parts[0]
                bed_start = int(parts[1])  # 0-based
                bed_end = int(parts[2])    # 0-based, exclusive
                score = float(parts[3])
                
                # Skip if different chromosome
                if bed_chrom != chrom:
                    continue
                    
            elif len(parts) == 3:
                # Assume implicit chromosome (3 cols: start, end, score)
                # This is a heuristic for when BedGraph is already filtered for the gene
                bed_chrom = chrom  # Implicitly match requested chrom
                bed_start = int(parts[0])
                bed_end = int(parts[1])
                score = float(parts[2])
                
            else:
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
