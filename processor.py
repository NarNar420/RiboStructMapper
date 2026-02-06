"""
Data Processor for RiboStructMapper.

This module provides functions to aggregate ribosome density data from nucleotide-level
to amino acid-level, and apply offset corrections to the resulting vectors.
"""

from typing import List, Dict, Literal
import numpy as np


def aggregate_density(
    raw_density_vector: List[float], 
    method: Literal['mean', 'max', 'sum', 'median'] = 'mean'
) -> List[float]:
    """
    Aggregate nucleotide-level density scores to amino acid-level scores.
    
    Each amino acid is encoded by a codon (3 nucleotides). This function groups
    the raw density scores into triplets and applies an aggregation method to
    produce one score per amino acid.
    
    Args:
        raw_density_vector: List of density scores at nucleotide resolution
        method: Aggregation method - 'mean', 'max', 'sum', or 'median'
        
    Returns:
        List of aggregated scores, one per amino acid (length = len(raw_density_vector) // 3)
        
    Raises:
        ValueError: If the length of raw_density_vector is not divisible by 3
        ValueError: If method is not one of the supported methods
        
    Example:
        >>> aggregate_density([1, 2, 3, 4, 5, 6], method='mean')
        [2.0, 5.0]
        >>> aggregate_density([1, 2, 3, 4, 5, 6], method='max')
        [3.0, 6.0]
    """
    if len(raw_density_vector) % 3 != 0:
        raise ValueError(
            f"Density vector length ({len(raw_density_vector)}) must be divisible by 3 "
            f"to aggregate codons into amino acids."
        )
    
    valid_methods = ['mean', 'max', 'sum', 'median']
    if method not in valid_methods:
        raise ValueError(f"Method must be one of {valid_methods}, got '{method}'")
    
    aa_scores = []
    
    # Iterate through the density vector in steps of 3 (codons)
    for i in range(0, len(raw_density_vector), 3):
        codon_scores = raw_density_vector[i:i+3]
        
        if method == 'mean':
            aggregated_score = np.mean(codon_scores)
        elif method == 'max':
            aggregated_score = np.max(codon_scores)
        elif method == 'sum':
            aggregated_score = np.sum(codon_scores)
        elif method == 'median':
            aggregated_score = np.median(codon_scores)
        
        aa_scores.append(float(aggregated_score))
    
    return aa_scores


def apply_offsets(
    aa_score_vector: List[float], 
    offset_values: List[int]
) -> Dict[int, List[float]]:
    """
    Apply multiple offset shifts to an amino acid score vector.
    
    Ribosome profiling data may have a phase offset due to the ribosome's position
    relative to the codon being translated. This function applies various offset
    shifts to test different hypotheses about the ribosome's position.
    
    Offset logic:
    - Positive offset (+N): Shift scores to the RIGHT (later positions)
      Example: offset=+1 means score at index 0 moves to index 1
    - Negative offset (-N): Shift scores to the LEFT (earlier positions)
      Example: offset=-1 means score at index 1 moves to index 0
    - Positions that fall outside bounds are filled with 0.0
    
    Args:
        aa_score_vector: List of scores at amino acid resolution
        offset_values: List of integer offsets to apply (in nucleotide positions,
                       will be converted to amino acid positions by dividing by 3)
        
    Returns:
        Dictionary mapping each offset value to its corresponding shifted vector
        
    Example:
        >>> apply_offsets([1.0, 2.0, 3.0], [0, -3])
        {0: [1.0, 2.0, 3.0], -3: [2.0, 3.0, 0.0]}
    """
    result = {}
    
    for offset_nt in offset_values:
        # Convert nucleotide offset to amino acid offset
        # offset_nt is in nucleotide positions (e.g., -15 nt = -5 aa)
        offset_aa = offset_nt // 3
        
        # Create a new vector filled with zeros
        shifted_vector = [0.0] * len(aa_score_vector)
        
        # Apply the shift
        for i in range(len(aa_score_vector)):
            # Calculate the source index
            source_index = i - offset_aa
            
            # If source index is within bounds, copy the score
            if 0 <= source_index < len(aa_score_vector):
                shifted_vector[i] = aa_score_vector[source_index]
            # Otherwise, leave as 0.0 (already initialized)
        
        result[offset_nt] = shifted_vector
    
    return result


__all__ = ["aggregate_density", "apply_offsets"]
