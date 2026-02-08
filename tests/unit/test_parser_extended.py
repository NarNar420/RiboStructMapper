#!/usr/bin/env python3
"""
Extended test script for parser.py - Testing PDB structure and density parsers.

This script tests the new parse_pdb_structure() and parse_ribo_density() 
functions added to parser.py.
"""

import os
import sys
import numpy as np
from ribostruct.core.parser import parse_pdb_structure, parse_ribo_density



def test_parse_pdb_structure_basic():
    """Test basic PDB structure parsing without chain specification."""
    print("=" * 60)
    print("Test 1: Basic parse_pdb_structure (no chain filter)")
    print("=" * 60)
    
    pdb_path = "data/mock/mock.pdb"
    
    try:
        sequence, residues = parse_pdb_structure(pdb_path)
        
        print(f"✓ Successfully parsed PDB structure")
        print(f"  Sequence: {sequence}")
        print(f"  Length: {len(sequence)} residues")
        print(f"  Number of residue objects: {len(residues)}")
        
        # Validate that sequence length matches residue count
        assert len(sequence) == len(residues), \
            f"Length mismatch: sequence={len(sequence)}, residues={len(residues)}"
        print(f"✓ Sequence and residue list lengths match: {len(sequence)}")
        
        # Expected: 15 residues (MKTIIALSY + IFCLVF)
        assert len(sequence) == 15, \
            f"Expected 15 residues, got {len(sequence)}"
        print(f"✓ Sequence length is correct (15 residues)")
        
        # Expected sequence from mock.pdb (extracted manually)
        expected_seq = "MKTIIALSY" "IFCLVF"  # 15 amino acids
        assert sequence == expected_seq, \
            f"Sequence mismatch:\\nExpected: {expected_seq}\\nGot:      {sequence}"
        print(f"✓ Sequence matches expected mock data")
        
        print("\n✓✓✓ Test 1 PASSED ✓✓✓\n")
        return True
        
    except Exception as e:
        print(f"✗ Test 1 FAILED: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_parse_pdb_structure_with_chain():
    """Test PDB structure parsing with explicit chain filter."""
    print("=" * 60)
    print("Test 2: parse_pdb_structure with chain_id='A'")
    print("=" * 60)
    
    pdb_path = "data/mock/mock.pdb"
    
    try:
        sequence, residues = parse_pdb_structure(pdb_path, chain_id='A')
        
        print(f"✓ Successfully parsed with chain_id='A'")
        print(f"  Sequence: {sequence}")
        print(f"  Length: {len(sequence)}")
        
        # Verify all residues belong to chain A
        for res in residues:
            chain = res.get_parent().id
            assert chain == 'A', f"Found residue in chain {chain}, expected A"
        print(f"✓ All residues belong to chain A")
        
        assert len(sequence) == 15
        print(f"✓ Sequence length is correct (15 residues)")
        
        print("\n✓✓✓ Test 2 PASSED ✓✓✓\n")
        return True
        
    except Exception as e:
        print(f"✗ Test 2 FAILED: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_parse_pdb_structure_invalid_chain():
    """Test error handling for invalid chain ID."""
    print("=" * 60)
    print("Test 3: Error handling - invalid chain ID")
    print("=" * 60)
    
    pdb_path = "data/mock/mock.pdb"
    
    try:
        sequence, residues = parse_pdb_structure(pdb_path, chain_id='Z')
        print(f"✗ Should have raised ValueError for invalid chain")
        return False
        
    except ValueError as e:
        print(f"✓ Correctly raised ValueError: {e}")
        print("\n✓✓✓ Test 3 PASSED ✓✓✓\n")
        return True
    except Exception as e:
        print(f"✗ Unexpected error: {e}")
        return False


def test_parse_pdb_structure_missing_file():
    """Test error handling for missing PDB file."""
    print("=" * 60)
    print("Test 4: Error handling - missing PDB file")
    print("=" * 60)
    
    try:
        sequence, residues = parse_pdb_structure("nonexistent.pdb")
        print(f"✗ Should have raised FileNotFoundError")
        return False
        
    except FileNotFoundError as e:
        print(f"✓ Correctly raised FileNotFoundError: {e}")
        print("\n✓✓✓ Test 4 PASSED ✓✓✓\n")
        return True
    except Exception as e:
        print(f"✗ Unexpected error: {e}")
        return False


def test_parse_ribo_density_basic():
    """Test basic bedGraph density parsing."""
    print("=" * 60)
    print("Test 5: Basic parse_ribo_density")
    print("=" * 60)
    
    bedgraph_path = "data/mock/mock.bedgraph"
    
    try:
        # Parse density for the full CDS region (1-45)
        density = parse_ribo_density(bedgraph_path, "mock_chrom", 1, 45)
        
        print(f"✓ Successfully parsed bedGraph density")
        print(f"  Density vector length: {len(density)}")
        print(f"  First 10 values: {density[:10]}")
        print(f"  Last 10 values: {density[-10:]}")
        
        # Should return 45 values (positions 1-45)
        assert len(density) == 45, \
            f"Expected 45 density values, got {len(density)}"
        print(f"✓ Density vector length is correct (45)")
        
        # Check that density values match the mock bedGraph
        # mock.bedgraph has:
        # mock_chrom 0  15  10.5  -> 1-based positions 1-15
        # mock_chrom 15 30  25.3  -> 1-based positions 16-30
        # mock_chrom 30 45  18.7  -> 1-based positions 31-45
        
        # Check first region (positions 1-15)
        assert np.all(density[0:15] == 10.5), \
            f"Expected first 15 values to be 10.5, got {density[0:15]}"
        print(f"✓ Positions 1-15 correctly filled with 10.5")
        
        # Check second region (positions 16-30)
        assert np.all(density[15:30] == 25.3), \
            f"Expected positions 16-30 to be 25.3, got {density[15:30]}"
        print(f"✓ Positions 16-30 correctly filled with 25.3")
        
        # Check third region (positions 31-45)
        assert np.all(density[30:45] == 18.7), \
            f"Expected positions 31-45 to be 18.7, got {density[30:45]}"
        print(f"✓ Positions 31-45 correctly filled with 18.7")
        
        print("\n✓✓✓ Test 5 PASSED ✓✓✓\n")
        return True
        
    except Exception as e:
        print(f"✗ Test 5 FAILED: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_parse_ribo_density_partial_overlap():
    """Test density parsing with partial coordinate overlap."""
    print("=" * 60)
    print("Test 6: parse_ribo_density with partial overlap")
    print("=" * 60)
    
    bedgraph_path = "data/mock/mock.bedgraph"
    
    try:
        # Request positions 20-35 (spans two bedGraph regions)
        density = parse_ribo_density(bedgraph_path, "mock_chrom", 20, 35)
        
        print(f"✓ Successfully parsed partial region")
        print(f"  Density vector: {density}")
        print(f"  Length: {len(density)}")
        
        # Should return 16 values (positions 20-35 inclusive)
        assert len(density) == 16, \
            f"Expected 16 density values, got {len(density)}"
        print(f"✓ Density vector length is correct (16)")
        
        # Positions 20-30 should be 25.3 (11 positions)
        # Positions 31-35 should be 18.7 (5 positions)
        assert np.all(density[0:11] == 25.3), \
            f"Expected first 11 values to be 25.3"
        assert np.all(density[11:16] == 18.7), \
            f"Expected last 5 values to be 18.7"
        print(f"✓ Partial overlap handled correctly")
        
        print("\n✓✓✓ Test 6 PASSED ✓✓✓\n")
        return True
        
    except Exception as e:
        print(f"✗ Test 6 FAILED: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_parse_ribo_density_no_overlap():
    """Test density parsing with no overlap (should return zeros)."""
    print("=" * 60)
    print("Test 7: parse_ribo_density with no overlap")
    print("=" * 60)
    
    bedgraph_path = "data/mock/mock.bedgraph"
    
    try:
        # Request positions way outside the bedGraph data (1000-1010)
        density = parse_ribo_density(bedgraph_path, "mock_chrom", 1000, 1010)
        
        print(f"✓ Successfully parsed non-overlapping region")
        print(f"  Density vector: {density}")
        
        # Should return all zeros
        assert len(density) == 11
        assert np.all(density == 0.0), \
            f"Expected all zeros, got {density}"
        print(f"✓ Non-overlapping region correctly filled with zeros")
        
        print("\n✓✓✓ Test 7 PASSED ✓✓✓\n")
        return True
        
    except Exception as e:
        print(f"✗ Test 7 FAILED: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_parse_ribo_density_missing_file():
    """Test error handling for missing bedGraph file."""
    print("=" * 60)
    print("Test 8: Error handling - missing bedGraph file")
    print("=" * 60)
    
    try:
        density = parse_ribo_density("nonexistent.bedgraph", "chr1", 1, 100)
        print(f"✗ Should have raised FileNotFoundError")
        return False
        
    except FileNotFoundError as e:
        print(f"✓ Correctly raised FileNotFoundError: {e}")
        print("\n✓✓✓ Test 8 PASSED ✓✓✓\n")
        return True
    except Exception as e:
        print(f"✗ Unexpected error: {e}")
        return False


def test_parse_ribo_density_invalid_coords():
    """Test error handling for invalid coordinates."""
    print("=" * 60)
    print("Test 9: Error handling - invalid coordinates")
    print("=" * 60)
    
    bedgraph_path = "data/mock/mock.bedgraph"
    
    try:
        # start > end
        density = parse_ribo_density(bedgraph_path, "mock_chrom", 100, 50)
        print(f"✗ Should have raised ValueError for start > end")
        return False
        
    except ValueError as e:
        print(f"✓ Correctly raised ValueError: {e}")
        print("\n✓✓✓ Test 9 PASSED ✓✓✓\n")
        return True
    except Exception as e:
        print(f"✗ Unexpected error: {e}")
        return False


def main():
    """Run all tests."""
    print("\n" + "=" * 60)
    print("EXTENDED PARSER TEST SUITE")
    print("Testing parse_pdb_structure() and parse_ribo_density()")
    print("=" * 60 + "\n")
    
    # Change to the project directory
    project_dir = os.path.dirname(os.path.abspath(__file__))
    os.chdir(project_dir)
    
    results = []
    
    # PDB structure parsing tests
    results.append(("Test 1: Basic PDB parse", test_parse_pdb_structure_basic()))
    results.append(("Test 2: PDB with chain", test_parse_pdb_structure_with_chain()))
    results.append(("Test 3: Invalid chain", test_parse_pdb_structure_invalid_chain()))
    results.append(("Test 4: Missing PDB file", test_parse_pdb_structure_missing_file()))
    
    # Ribosome density parsing tests
    results.append(("Test 5: Basic density parse", test_parse_ribo_density_basic()))
    results.append(("Test 6: Partial overlap", test_parse_ribo_density_partial_overlap()))
    results.append(("Test 7: No overlap", test_parse_ribo_density_no_overlap()))
    results.append(("Test 8: Missing bedGraph", test_parse_ribo_density_missing_file()))
    results.append(("Test 9: Invalid coords", test_parse_ribo_density_invalid_coords()))
    
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
