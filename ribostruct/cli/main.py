#!/usr/bin/env python3
"""
RiboStructMapper CLI - End-to-End Pipeline Integration

This script demonstrates the complete RiboStructMapper pipeline by processing
genomic data, ribosome density, and PDB structure files to generate modified
PDB files with density scores in the B-factor field.

Usage:
    python main_cli.py

This will process the mock data files and generate output PDB files for
different offset values.
"""

import os
import sys
from typing import List, Dict
from Bio.PDB import PDBParser

# Import all modules
from ribostruct.core.parser import parse_genomic_data, parse_pdb_structure, parse_ribo_density
from ribostruct.core.alignment import translate_sequence, align_sequences
from ribostruct.core.processor import process_offsets
from ribostruct.core.injector import inject_bfactors



def run_pipeline(
    pdb_path: str,
    fasta_path: str,
    bedgraph_path: str = None,
    offsets: List[int] = None,
    chain_id: str = None,
    aggregation_method: str = 'mean',
    output_dir: str = '.'
) -> Dict[int, str]:
    """
    Run the complete RiboStructMapper pipeline.
    
    Args:
        pdb_path: Path to input PDB file
        fasta_path: Path to FASTA file (header format: >ID Chrom:Start..End)
        bedgraph_path: Path to bedGraph density file
        offsets: List of offset values (in nucleotides) to apply
        chain_id: Optional PDB chain ID (default: first chain)
        aggregation_method: Density aggregation method ('mean', 'max', 'sum', 'median')
        output_dir: Directory to save output files (default: current directory)
    
    Returns:
        Dictionary mapping offset values to output file paths
    """
    print("=" * 70)
    print("RIBOSTRUCTMAPPER - END-TO-END PIPELINE")
    print("=" * 70)
    
    # ========================================================================
    # STEP 1: Parse Genomic Data (FASTA Header for Coordinates)
    # ========================================================================
    print("\n[STEP 1] Parsing genomic data...")
    print(f"  FASTA: {fasta_path}")
    
    nucleotide_seq, coord_map, start_offset, chrom, genomic_start = parse_genomic_data(
        fasta_path
    )
    
    print(f"  ✓ Extracted CDS: {len(nucleotide_seq)} nucleotides")
    print(f"    Chromosome: {chrom}")
    print(f"    Genomic Start: {genomic_start} (1-based, inclusive)")
    if start_offset > 0:
        print(f"    ⚠ Sanitized: Removed {start_offset} nucleotides from start")
    print(f"    Sequence: {nucleotide_seq[:60]}...")
    
    # ========================================================================
    # STEP 2: Parse Ribosome Density (bedGraph)
    # ========================================================================
    print("\n[STEP 2] Parsing ribosome density...")
    print(f"  bedGraph: {bedgraph_path}")
    
    # Calculate end coordinate based on sequence length
    # start is 1-based inclusive, so end = start + length - 1
    # Example: Start 100, Len 3 -> 100, 101, 102 (End 102) -> 100 + 3 - 1 = 102
    end = genomic_start + len(nucleotide_seq) - 1
    
    raw_density = parse_ribo_density(bedgraph_path, chrom, genomic_start, end)
    
    print(f"  ✓ Loaded raw density: {len(raw_density)} values")
    print(f"    Chromosome: {chrom}, Region: {genomic_start}-{end}")
    print(f"    Density range: [{raw_density.min():.1f}, {raw_density.max():.1f}]")
    
    # Density vector is already aligned to the sanitized sequence because
    # we fetched it using the *final* genomic start and length.
    aligned_density = raw_density
    
    # Safety check
    assert len(aligned_density) == len(nucleotide_seq), \
        f"Density length ({len(aligned_density)}) must match sequence length ({len(nucleotide_seq)})"
    
    print(f"  ✓ Final density: {len(aligned_density)} values (matches sequence)")
    
    
    # ========================================================================
    # STEP 3: Parse PDB Structure
    # ========================================================================
    print("\n[STEP 3] Parsing PDB structure...")
    print(f"  PDB: {pdb_path}")
    
    pdb_aa_sequence, pdb_residues = parse_pdb_structure(pdb_path, chain_id=chain_id)
    
    # Also load the full structure for later injection
    parser = PDBParser(QUIET=True)
    pdb_structure = parser.get_structure('protein', pdb_path)
    
    print(f"  ✓ Extracted sequence: {len(pdb_aa_sequence)} amino acids")
    print(f"    Sequence: {pdb_aa_sequence}")
    if chain_id:
        print(f"    Chain: {chain_id}")
    
    # ========================================================================
    # STEP 4: Translate Genomic Sequence
    # ========================================================================
    print("\n[STEP 4] Translating genomic sequence...")
    
    genomic_aa = translate_sequence(nucleotide_seq)
    
    print(f"  ✓ Translated: {len(genomic_aa)} amino acids")
    print(f"    Sequence: {genomic_aa}")
    
    # ========================================================================
    # STEP 5: Align Sequences (Genomic vs PDB)
    # ========================================================================
    print("\n[STEP 5] Aligning sequences...")
    
    alignment_score, mapping_index = align_sequences(genomic_aa, pdb_aa_sequence)
    
    print(f"  ✓ Alignment complete")
    print(f"    Score: {alignment_score:.1f}")
    print(f"    Mapped residues: {len(mapping_index)}/{len(genomic_aa)}")
    print(f"    Coverage: {len(mapping_index)/len(genomic_aa)*100:.1f}%")
    
    # ========================================================================
    # STEP 6: Process Offsets (NT-level shift + AA aggregation)
    # ========================================================================
    print(f"\n[STEP 6] Processing offsets at nucleotide level...")
    print(f"  Offsets: {offsets}")
    print(f"  Aggregation method: {aggregation_method}")
    
    # CRITICAL: Apply offsets BEFORE aggregation (at NT level, not AA level)
    offset_scores = process_offsets(aligned_density, offsets, method=aggregation_method)
    
    print(f"  ✓ Generated {len(offset_scores)} offset versions")
    for offset, scores in offset_scores.items():
        print(f"    Offset {offset:4d} nt: {len(scores)} AA scores, range [{scores.min():.1f}, {scores.max():.1f}]")
    
    # ========================================================================
    # STEP 7: Inject B-factors and Save Output Files
    # ========================================================================
    print(f"\n[STEP 7] Injecting B-factors and saving output files...")
    
    output_files = {}
    
    for offset_value in offsets:
        # Get the score vector for this offset
        score_vector = offset_scores[offset_value]
        
        # Generate output filename
        output_filename = f"output_offset_{offset_value}.pdb"
        output_path = os.path.join(output_dir, output_filename)
        
        # Inject B-factors
        # Note: We need to reload the structure for each offset to avoid modifying the same object
        fresh_structure = parser.get_structure('protein', pdb_path)
        inject_bfactors(fresh_structure, mapping_index, score_vector, output_path)
        
        output_files[offset_value] = output_path
        print(f"  ✓ Saved: {output_filename}")
    
    # ========================================================================
    # SUMMARY
    # ========================================================================
    print("\n" + "=" * 70)
    print("PIPELINE COMPLETE")
    print("=" * 70)
    print(f"Generated {len(output_files)} output files:")
    for offset, path in sorted(output_files.items()):
        print(f"  Offset {offset:4d} nt: {path}")
    print("=" * 70 + "\n")
    
    return output_files


def verify_output(output_path: str) -> None:
    """
    Verify that an output PDB file has non-zero B-factors.
    
    Args:
        output_path: Path to PDB file to verify
    """
    print(f"\nVerifying output: {output_path}")
    
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('output', output_path)
    
    bfactors = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[0] == ' ':
                    for atom in residue:
                        bfactors.append(atom.get_bfactor())
        break
    
    non_zero_count = sum(1 for bf in bfactors if abs(bf) > 0.01)
    
    print(f"  Total atoms: {len(bfactors)}")
    print(f"  Non-zero B-factors: {non_zero_count} ({non_zero_count/len(bfactors)*100:.1f}%)")
    print(f"  B-factor range: [{min(bfactors):.2f}, {max(bfactors):.2f}]")
    
    if non_zero_count > 0:
        print(f"  ✓ Output file has non-zero B-factors")
    else:
        print(f"  ⚠ Warning: All B-factors are zero")


def main():
    """Main execution function."""
    print("\nRiboStructMapper CLI - Running with mock data\n")
    
    # Define paths to mock data
    pdb_path = "data/mock/mock.pdb"
    fasta_path = "data/mock/mock.fasta"
    # Run the pipeline
    output_files = run_pipeline(
        pdb_path=pdb_path,
        fasta_path=fasta_path,
        bedgraph_path=bedgraph_path,
        offsets=offsets,
        aggregation_method='mean',
        output_dir='.'
    )
    
    # Verify the output files
    print("\n" + "=" * 70)
    print("VERIFICATION")
    print("=" * 70)
    
    for offset, path in sorted(output_files.items()):
        verify_output(path)
    
    print("\n✓ Pipeline execution complete!\n")


if __name__ == "__main__":
    main()
