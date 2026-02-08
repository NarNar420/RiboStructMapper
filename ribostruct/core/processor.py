"""
RiboStructMapper - Data Processor Module

This module handles ribosome density data processing, including:
- Nucleotide-level offset application (ribosome position correction)
- Nucleotide-to-amino acid aggregation
- Multiple aggregation methods (mean, max, sum, median)

CRITICAL: Offsets are applied at the NUCLEOTIDE level BEFORE aggregation
to amino acids. This ensures proper spatial mapping of ribosome positions.
"""

import numpy as np
from typing import Dict, List


def shift_nucleotide_density(density: np.ndarray, offset: int) -> np.ndarray:
    """
    Shift nucleotide density values by the specified offset.
    
    This function applies ribosome position correction at the nucleotide level.
    Negative offsets shift density upstream (towards 5' end).
    
    Args:
        density: Nucleotide-level density array
        offset: Number of nucleotides to shift (negative = upstream)
    
    Returns:
        Shifted density array with zero-filled edges (no wrapping)
    
    Example:
        Input:  [10, 20, 30, 40, 50]
        Offset: -2 (shift upstream by 2 nt)
        Output: [30, 40, 50, 0, 0]
        
        The value at position i moves to position i+offset.
        For offset=-2: position 2 moves to position 0, etc.
    """
    if offset == 0:
        return density.copy()
    
    # Create output array filled with zeros
    shifted = np.zeros_like(density, dtype=float)
    
    if offset < 0:
        # Negative offset: shift left (upstream)
        # Values from positions [-offset:] move to positions [0:len+offset]
        shifted[:len(density) + offset] = density[-offset:]
    else:
        # Positive offset: shift right (downstream)
        # Values from positions [:-offset] move to positions [offset:]
        shifted[offset:] = density[:-offset]
    
    return shifted


def aggregate_density(
    nucleotide_density: np.ndarray,
    method: str = 'mean'
) -> np.ndarray:
    """
    Aggregate nucleotide-level density to amino acid level.
    
    Triplets of nucleotides (codons) are aggregated using the specified method.
    
    Args:
        nucleotide_density: Array of nucleotide density values (length must be divisible by 3)
        method: Aggregation method ('mean', 'max', 'sum', 'median')
    
    Returns:
        Amino acid-level density array (length = len(nucleotide_density) / 3)
    
    Raises:
        ValueError: If nucleotide density length is not divisible by 3
        ValueError: If method is not supported
    
    Example:
        Input: [10, 20, 30, 40, 50, 60]  (2 codons)
        Method: 'mean'
        Output: [20.0, 50.0]  (mean of each triplet)
    """
    if len(nucleotide_density) % 3 != 0:
        raise ValueError(
            f"Nucleotide density length ({len(nucleotide_density)}) "
            f"must be divisible by 3"
        )
    
    # Reshape into codons (N/3 rows, 3 columns)
    n_codons = len(nucleotide_density) // 3
    codons = nucleotide_density.reshape(n_codons, 3)
    
    # Apply aggregation method
    if method == 'mean':
        aa_density = np.mean(codons, axis=1)
    elif method == 'max':
        aa_density = np.max(codons, axis=1)
    elif method == 'sum':
        aa_density = np.sum(codons, axis=1)
    elif method == 'median':
        aa_density = np.median(codons, axis=1)
    else:
        raise ValueError(
            f"Unsupported aggregation method: {method}. "
            f"Use 'mean', 'max', 'sum', or 'median'."
        )
    
    return aa_density


def process_offsets(
    raw_density: np.ndarray,
    offsets: List[int],
    method: str = 'mean'
) -> Dict[int, np.ndarray]:
    """
    Process multiple offset values and return amino acid densities for each.
    
    This is the main processing function that:
    1. Applies each offset to the raw nucleotide density
    2. Aggregates shifted nucleotide density to amino acid level
    3. Returns a dictionary mapping offset → AA density
    
    Args:
        raw_density: Raw nucleotide-level density array
        offsets: List of offset values to apply
        method: Aggregation method for nucleotide→AA conversion
    
    Returns:
        Dictionary mapping offset value to amino acid density array
    
    Example:
        raw_density = np.array([10, 10, 10, 20, 20, 20])
        offsets = [0, -3]
        method = 'mean'
        
        Result:
        {
            0: [10.0, 20.0],      # No shift: codon1=(10,10,10), codon2=(20,20,20)
            -3: [16.67, 6.67]     # Shift -3: codon1=(20,20,20), codon2=(0,0,0)
        }
    """
    results = {}
    
    for offset in offsets:
        # Step 1: Shift nucleotide density
        shifted_nt_density = shift_nucleotide_density(raw_density, offset)
        
        # Step 2: Aggregate to amino acid level
        aa_density = aggregate_density(shifted_nt_density, method=method)
        
        # Store result
        results[offset] = aa_density
    
    return results


# Legacy function name for backward compatibility
def apply_offsets(
    amino_acid_scores: np.ndarray,
    offsets: List[int]
) -> Dict[int, np.ndarray]:
    """
    DEPRECATED: This function is kept for backward compatibility only.
    
    New code should use process_offsets() which correctly applies offsets
    at the nucleotide level before aggregation.
    
    This function incorrectly applies offsets at the amino acid level.
    """
    import warnings
    warnings.warn(
        "apply_offsets() is deprecated and applies offsets incorrectly at the AA level. "
        "Use process_offsets() instead, which applies offsets at the NT level.",
        DeprecationWarning,
        stacklevel=2
    )
    
    results = {}
    for offset in offsets:
        # Old (incorrect) logic: shift AA scores
        shifted = np.zeros_like(amino_acid_scores, dtype=float)
        aa_offset = offset // 3
        
        if aa_offset == 0:
            shifted = amino_acid_scores.copy()
        elif aa_offset < 0:
            shifted[:len(amino_acid_scores) + aa_offset] = amino_acid_scores[-aa_offset:]
        else:
            shifted[aa_offset:] = amino_acid_scores[:-aa_offset]
        
        results[offset] = shifted
    
    return results
