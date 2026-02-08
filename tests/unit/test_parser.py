#!/usr/bin/env python3
"""
Test script for parser.py using mock data.
Tests the parse_genomic_data function with mock FASTA and GTF files.
"""

import os
import sys
from ribostruct.core.parser import parse_genomic_data


def test_parse_genomic_data_basic():
    """Test basic parsing of genomic data with mock files."""
    print("=" * 60)
    print("Test 1: Basic parse_genomic_data with mock files")
    print("=" * 60)
    
    fasta_path = "data/mock/mock.fasta"
    gtf_path = "data/mock/mock.gtf"
    
    try:
        nucleotide_seq, coord_map = parse_genomic_data(fasta_path, gtf_path)
        
        print(f"✓ Successfully parsed genomic data")
        print(f"  Nucleotide sequence length: {len(nucleotide_seq)}")
        print(f"  Nucleotide sequence: {nucleotide_seq}")
        print(f"  Coordinate map length: {len(coord_map)}")
        print(f"  Coordinate map: {coord_map}")
        
        # Validate that sequence length matches coordinate map
        assert len(nucleotide_seq) == len(coord_map), \
            f"Length mismatch: seq={len(nucleotide_seq)}, coord_map={len(coord_map)}"
        print(f"✓ Sequence and coordinate map lengths match: {len(nucleotide_seq)}")
        
        # Expected: 45 nucleotides (positions 1-45 in GTF)
        assert len(nucleotide_seq) == 45, \
            f"Expected 45 nucleotides, got {len(nucleotide_seq)}"
        print(f"✓ Sequence length is correct (45 bp)")
        
        # Expected sequence from mock.fasta
        expected_seq = "ATGAAAACCATAATAGCTCTGTCTTACATATTCTGCCTGGTGTTC"
        assert nucleotide_seq == expected_seq, \
            f"Sequence mismatch:\nExpected: {expected_seq}\nGot:      {nucleotide_seq}"
        print(f"✓ Sequence matches expected mock data")
        
        print("\n✓✓✓ Test 1 PASSED ✓✓✓\n")
        return True
        
    except Exception as e:
        print(f"✗ Test 1 FAILED: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_parse_genomic_data_with_gene_id():
    """Test parsing with explicit gene_id filter."""
    print("=" * 60)
    print("Test 2: parse_genomic_data with gene_id filter")
    print("=" * 60)
    
    fasta_path = "data/mock/mock.fasta"
    gtf_path = "data/mock/mock.gtf"
    
    try:
        nucleotide_seq, coord_map = parse_genomic_data(
            fasta_path, gtf_path, gene_id="GENE1"
        )
        
        print(f"✓ Successfully parsed with gene_id='GENE1'")
        print(f"  Sequence length: {len(nucleotide_seq)}")
        assert len(nucleotide_seq) == 45
        print(f"✓ Sequence length is correct (45 bp)")
        
        print("\n✓✓✓ Test 2 PASSED ✓✓✓\n")
        return True
        
    except Exception as e:
        print(f"✗ Test 2 FAILED: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_parse_genomic_data_with_transcript_id():
    """Test parsing with explicit transcript_id filter."""
    print("=" * 60)
    print("Test 3: parse_genomic_data with transcript_id filter")
    print("=" * 60)
    
    fasta_path = "data/mock/mock.fasta"
    gtf_path = "data/mock/mock.gtf"
    
    try:
        nucleotide_seq, coord_map = parse_genomic_data(
            fasta_path, gtf_path, transcript_id="TX1"
        )
        
        print(f"✓ Successfully parsed with transcript_id='TX1'")
        print(f"  Sequence length: {len(nucleotide_seq)}")
        assert len(nucleotide_seq) == 45
        print(f"✓ Sequence length is correct (45 bp)")
        
        print("\n✓✓✓ Test 3 PASSED ✓✓✓\n")
        return True
        
    except Exception as e:
        print(f"✗ Test 3 FAILED: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_missing_fasta_file():
    """Test error handling for missing FASTA file."""
    print("=" * 60)
    print("Test 4: Error handling - missing FASTA file")
    print("=" * 60)
    
    try:
        nucleotide_seq, coord_map = parse_genomic_data(
            "nonexistent.fasta", 
            "data/mock/mock.gtf"
        )
        print(f"✗ Should have raised FileNotFoundError")
        return False
        
    except FileNotFoundError as e:
        print(f"✓ Correctly raised FileNotFoundError: {e}")
        print("\n✓✓✓ Test 4 PASSED ✓✓✓\n")
        return True
    except Exception as e:
        print(f"✗ Unexpected error: {e}")
        return False


def test_missing_gtf_file():
    """Test error handling for missing GTF file."""
    print("=" * 60)
    print("Test 5: Error handling - missing GTF file")
    print("=" * 60)
    
    try:
        nucleotide_seq, coord_map = parse_genomic_data(
            "data/mock/mock.fasta",
            "nonexistent.gtf"
        )
        print(f"✗ Should have raised FileNotFoundError")
        return False
        
    except FileNotFoundError as e:
        print(f"✓ Correctly raised FileNotFoundError: {e}")
        print("\n✓✓✓ Test 5 PASSED ✓✓✓\n")
        return True
    except Exception as e:
        print(f"✗ Unexpected error: {e}")
        return False


def main():
    """Run all tests."""
    print("\n" + "=" * 60)
    print("PARSER TEST SUITE")
    print("=" * 60 + "\n")
    
    # Change to the project directory
    project_dir = os.path.dirname(os.path.abspath(__file__))
    os.chdir(project_dir)
    
    results = []
    results.append(("Test 1: Basic parse", test_parse_genomic_data_basic()))
    results.append(("Test 2: With gene_id", test_parse_genomic_data_with_gene_id()))
    results.append(("Test 3: With transcript_id", test_parse_genomic_data_with_transcript_id()))
    results.append(("Test 4: Missing FASTA", test_missing_fasta_file()))
    results.append(("Test 5: Missing GTF", test_missing_gtf_file()))
    
    # Summary
    print("\n" + "=" * 60)
    print("TEST SUMMARY")
    print("=" * 60)
    passed = sum(1 for _, result in results if result)
    total = len(results)
    
    for name, result in results:
        status = "PASS" if result else "FAIL"
        print(f"{status}: {name}")
    
    print(f"\nTotal: {passed}/{total} tests passed")
    print("=" * 60 + "\n")
    
    return all(result for _, result in results)


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
