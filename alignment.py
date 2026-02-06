"""
Alignment Engine for RiboStructMapper.

This module provides functions to translate nucleotide sequences to amino acids
and perform pairwise alignment between genomic and PDB sequences, generating
a mapping dictionary for downstream B-factor injection.
"""

from typing import Dict
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner


def translate_sequence(nucleotide_seq: str) -> str:
    """
    Translate a DNA/RNA nucleotide sequence to an amino acid sequence.
    
    Args:
        nucleotide_seq: DNA or RNA sequence string
        
    Returns:
        Amino acid sequence as a string
        
    Example:
        >>> translate_sequence("ATGAAACCC")
        'MKP'
    """
    seq_obj = Seq(nucleotide_seq)
    aa_seq = seq_obj.translate()
    return str(aa_seq)


def align_sequences(genomic_aa: str, pdb_aa: str) -> tuple[float, Dict[int, int]]:
    """
    Perform global pairwise alignment between genomic (theoretical) and PDB amino acid sequences.
    
    This function uses Biopython's PairwiseAligner with custom scoring to handle gaps
    and mismatches that may occur due to missing residues in the PDB structure or 
    sequence variations.
    
    Args:
        genomic_aa: Amino acid sequence translated from the genome
        pdb_aa: Amino acid sequence extracted from the PDB structure
        
    Returns:
        tuple containing:
            - alignment_score (float): The best alignment score
            - mapping_dict (Dict[int, int]): Dictionary mapping genomic indices (0-based) 
              to PDB residue indices (0-based). Gaps in the PDB sequence are NOT included
              in the mapping.
              
    Example:
        >>> genomic = "MKTIIALSY"
        >>> pdb = "MKT-IALSY"  # Missing one residue
        >>> score, mapping = align_sequences(genomic, pdb)
        >>> # mapping[3] will NOT exist because it's a gap in PDB
    """
    # Configure PairwiseAligner with specified parameters
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    
    # Perform alignment
    alignments = aligner.align(genomic_aa, pdb_aa)
    best_alignment = alignments[0]  # Get the best alignment
    alignment_score = best_alignment.score
    
    # Extract aligned sequences
    aligned_genomic = best_alignment[0]  # First sequence with gaps
    aligned_pdb = best_alignment[1]      # Second sequence with gaps
    
    # Build mapping dictionary
    # Key: genomic index (0-based in original genomic_aa)
    # Value: PDB residue index (0-based in original pdb_aa)
    mapping_dict: Dict[int, int] = {}
    
    genomic_idx = 0
    pdb_idx = 0
    
    for i in range(len(aligned_genomic)):
        genomic_char = aligned_genomic[i]
        pdb_char = aligned_pdb[i]
        
        # If both positions have residues (not gaps), record mapping
        if genomic_char != '-' and pdb_char != '-':
            mapping_dict[genomic_idx] = pdb_idx
            genomic_idx += 1
            pdb_idx += 1
        # If gap in genomic but residue in PDB (insertion in PDB)
        elif genomic_char == '-' and pdb_char != '-':
            pdb_idx += 1
        # If residue in genomic but gap in PDB (deletion in PDB / missing residue)
        elif genomic_char != '-' and pdb_char == '-':
            # Do NOT add to mapping - this genomic residue has no PDB counterpart
            genomic_idx += 1
        # Both gaps (shouldn't happen in global alignment, but handle it)
        # No index increment needed
    
    return alignment_score, mapping_dict


__all__ = ["translate_sequence", "align_sequences"]
