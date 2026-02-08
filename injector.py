"""
PDB Injector for RiboStructMapper.

This module provides functionality to inject ribosome density scores into
PDB B-factor fields for visualization in molecular viewers like PyMOL and Chimera.
"""

import os
from typing import Dict, List, Optional
from Bio.PDB import Structure, PDBIO


def inject_bfactors(
    structure: Structure.Structure,
    mapping_index: Dict[int, int],
    score_vector: List[float],
    output_path: str,
    unmapped_value: float = 0.0
) -> None:
    """
    Inject ribosome density scores into PDB B-factor fields.
    
    This function takes a PDB structure, a mapping from genomic indices to PDB
    residue indices, and a vector of scores. It sets the B-factor of every atom
    in each mapped residue to the corresponding score. Unmapped residues are
    set to a default value (typically 0.0).
    
    Args:
        structure: Bio.PDB.Structure object (typically from parse_pdb_structure)
        mapping_index: Dictionary mapping genomic AA indices (0-based) to PDB 
                      residue indices (0-based). This comes from align_sequences().
        score_vector: List of density scores at amino acid resolution.
                     Index corresponds to genomic AA position.
        output_path: Path where the modified PDB file will be saved
        unmapped_value: B-factor value to assign to residues not in the mapping
                       (default: 0.0)
    
    Returns:
        None (writes to output_path)
    
    Raises:
        ValueError: If output directory doesn't exist or scores are invalid
        IndexError: If mapping references invalid indices
    
    Example:
        >>> from parser import parse_pdb_structure
        >>> from alignment import align_sequences
        >>> from processor import aggregate_density
        >>> 
        >>> # Load structure
        >>> seq, residues = parse_pdb_structure("protein.pdb")
        >>> structure = parser.get_structure('protein', 'protein.pdb')
        >>> 
        >>> # Get alignment mapping
        >>> score, mapping = align_sequences(genomic_aa, seq)
        >>> 
        >>> # Get scores
        >>> scores = aggregate_density(raw_density, method='mean')
        >>> 
        >>> # Inject B-factors
        >>> inject_bfactors(structure, mapping, scores, "output.pdb")
    """
    # Validate output directory exists
    output_dir = os.path.dirname(output_path)
    if output_dir and not os.path.exists(output_dir):
        raise ValueError(f"Output directory does not exist: {output_dir}")
    
    # Build a list of all residues in the structure (in order)
    # We need to match them to the PDB indices from parse_pdb_structure
    all_residues = []
    for model in structure:
        for chain in model:
            for residue in chain:
                # Only consider standard amino acids (hetero flag is blank space)
                if residue.id[0] == ' ':
                    all_residues.append(residue)
        # Only process the first model
        break
    
    # Create reverse mapping: PDB residue index -> genomic index
    pdb_to_genomic = {pdb_idx: genomic_idx 
                      for genomic_idx, pdb_idx in mapping_index.items()}
    
    # Iterate through all residues and set B-factors
    for pdb_idx, residue in enumerate(all_residues):
        # Determine the score for this residue
        if pdb_idx in pdb_to_genomic:
            genomic_idx = pdb_to_genomic[pdb_idx]
            
            # Validate genomic index
            if genomic_idx < 0 or genomic_idx >= len(score_vector):
                raise IndexError(
                    f"Genomic index {genomic_idx} from mapping is out of bounds "
                    f"for score_vector of length {len(score_vector)}"
                )
            
            score = score_vector[genomic_idx]
        else:
            # Residue not in mapping - use unmapped value
            score = unmapped_value
        
        # Set B-factor for all atoms in this residue
        for atom in residue:
            atom.set_bfactor(score)
    
    # Save the modified structure to a new PDB file
    io = PDBIO()
    io.set_structure(structure)
    io.save(output_path)


__all__ = ["inject_bfactors"]
