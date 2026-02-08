#!/usr/bin/env python3
"""
Test script for the PDB Injector module.

This script verifies that:
1. B-factors are correctly injected into PDB structures
2. Mapped residues receive the correct scores
3. Unmapped residues are set to 0.0
4. All atoms in a residue get the same B-factor
"""

import os
import sys
import tempfile
from Bio.PDB import PDBParser
from ribostruct.core.parser import parse_pdb_structure
from ribostruct.core.injector import inject_bfactors


def test_basic_injection():
    """Test basic B-factor injection with simple mapping."""
    print("=" * 60)
    print("Test 1: Basic B-factor injection")
    print("=" * 60)
    
    pdb_path = "data/mock/mock.pdb"
    
    try:
        # Parse the structure
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('test', pdb_path)
        sequence, residues = parse_pdb_structure(pdb_path)
        
        print(f"✓ Loaded PDB with {len(residues)} residues")
        print(f"  Sequence: {sequence}")
        
        # Create a simple mapping and score vector
        # Map first 3 residues to specific scores
        mapping = {
            0: 0,  # Genomic index 0 -> PDB index 0
            1: 1,  # Genomic index 1 -> PDB index 1
            2: 2,  # Genomic index 2 -> PDB index 2
        }
        scores = [99.9, 88.8, 77.7]
        
        print(f"\n  Mapping: {mapping}")
        print(f"  Scores: {scores}")
        
        # Create temporary output file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as tmp:
            output_path = tmp.name
        
        try:
            # Inject B-factors
            inject_bfactors(structure, mapping, scores, output_path)
            print(f"✓ B-factors injected to temporary file")
            
            # Re-parse the output file
            new_parser = PDBParser(QUIET=True)
            new_structure = new_parser.get_structure('output', output_path)
            
            # Extract residues from new structure
            new_residues = []
            for model in new_structure:
                for chain in model:
                    for residue in chain:
                        if residue.id[0] == ' ':
                            new_residues.append(residue)
                break
            
            print(f"✓ Reloaded modified PDB with {len(new_residues)} residues")
            
            # Verify B-factors for mapped residues
            print(f"\n  Verifying B-factors for mapped residues:")
            for genomic_idx, pdb_idx in mapping.items():
                expected_score = scores[genomic_idx]
                residue = new_residues[pdb_idx]
                
                # Check all atoms in this residue
                atom_bfactors = [atom.get_bfactor() for atom in residue]
                
                print(f"    PDB residue {pdb_idx}: atoms={len(atom_bfactors)}, "
                      f"bfactors={atom_bfactors[:3]}...")
                
                # All atoms should have the same B-factor
                assert all(abs(bf - expected_score) < 0.01 for bf in atom_bfactors), \
                    f"B-factors don't match expected {expected_score}: {atom_bfactors}"
            
            print(f"  ✓ All mapped residues have correct B-factors")
            
            # Verify unmapped residues have B-factor of 0.0
            print(f"\n  Verifying B-factors for unmapped residues:")
            unmapped_count = 0
            for i in range(len(new_residues)):
                if i not in mapping.values():
                    residue = new_residues[i]
                    atom_bfactors = [atom.get_bfactor() for atom in residue]
                    
                    # Should all be 0.0
                    assert all(abs(bf) < 0.01 for bf in atom_bfactors), \
                        f"Unmapped residue {i} has non-zero B-factors: {atom_bfactors}"
                    unmapped_count += 1
            
            print(f"    Checked {unmapped_count} unmapped residues: all set to 0.0")
            print(f"  ✓ All unmapped residues have B-factor = 0.0")
            
            print("\n✓✓✓ Test 1 PASSED ✓✓✓\n")
            return True
            
        finally:
            # Clean up temporary file
            if os.path.exists(output_path):
                os.remove(output_path)
        
    except Exception as e:
        print(f"✗ Test 1 FAILED: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_all_residues_mapped():
    """Test when all residues are mapped."""
    print("=" * 60)
    print("Test 2: All residues mapped")
    print("=" * 60)
    
    pdb_path = "data/mock/mock.pdb"
    
    try:
        # Parse the structure
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('test', pdb_path)
        sequence, residues = parse_pdb_structure(pdb_path)
        
        num_residues = len(residues)
        print(f"✓ Loaded PDB with {num_residues} residues")
        
        # Create mapping for ALL residues
        mapping = {i: i for i in range(num_residues)}
        scores = [float(i * 10) for i in range(num_residues)]
        
        print(f"  Mapping: all {num_residues} residues")
        print(f"  Scores: {scores[:5]}... (first 5)")
        
        # Create temporary output file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as tmp:
            output_path = tmp.name
        
        try:
            # Inject B-factors
            inject_bfactors(structure, mapping, scores, output_path)
            print(f"✓ B-factors injected")
            
            # Re-parse and verify
            new_parser = PDBParser(QUIET=True)
            new_structure = new_parser.get_structure('output', output_path)
            
            new_residues = []
            for model in new_structure:
                for chain in model:
                    for residue in chain:
                        if residue.id[0] == ' ':
                            new_residues.append(residue)
                break
            
            # Verify all residues
            for i in range(num_residues):
                expected = scores[i]
                residue = new_residues[i]
                actual = list(residue)[0].get_bfactor()  # Get first atom's B-factor
                
                assert abs(actual - expected) < 0.01, \
                    f"Residue {i}: expected {expected}, got {actual}"
            
            print(f"✓ All {num_residues} residues have correct B-factors")
            print("\n✓✓✓ Test 2 PASSED ✓✓✓\n")
            return True
            
        finally:
            if os.path.exists(output_path):
                os.remove(output_path)
        
    except Exception as e:
        print(f"✗ Test 2 FAILED: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_no_residues_mapped():
    """Test when no residues are mapped (all should be 0.0)."""
    print("=" * 60)
    print("Test 3: No residues mapped")
    print("=" * 60)
    
    pdb_path = "data/mock/mock.pdb"
    
    try:
        # Parse the structure
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('test', pdb_path)
        sequence, residues = parse_pdb_structure(pdb_path)
        
        num_residues = len(residues)
        print(f"✓ Loaded PDB with {num_residues} residues")
        
        # Empty mapping
        mapping = {}
        scores = []
        
        print(f"  Mapping: empty (no residues mapped)")
        
        # Create temporary output file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as tmp:
            output_path = tmp.name
        
        try:
            # Inject B-factors
            inject_bfactors(structure, mapping, scores, output_path)
            print(f"✓ B-factors injected")
            
            # Re-parse and verify all are 0.0
            new_parser = PDBParser(QUIET=True)
            new_structure = new_parser.get_structure('output', output_path)
            
            new_residues = []
            for model in new_structure:
                for chain in model:
                    for residue in chain:
                        if residue.id[0] == ' ':
                            new_residues.append(residue)
                break
            
            # Verify all residues are 0.0
            for i, residue in enumerate(new_residues):
                for atom in residue:
                    actual = atom.get_bfactor()
                    assert abs(actual) < 0.01, \
                        f"Residue {i}, atom {atom.name}: expected 0.0, got {actual}"
            
            print(f"✓ All {num_residues} residues have B-factor = 0.0")
            print("\n✓✓✓ Test 3 PASSED ✓✓✓\n")
            return True
            
        finally:
            if os.path.exists(output_path):
                os.remove(output_path)
        
    except Exception as e:
        print(f"✗ Test 3 FAILED: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_custom_unmapped_value():
    """Test using a custom unmapped value instead of 0.0."""
    print("=" * 60)
    print("Test 4: Custom unmapped value")
    print("=" * 60)
    
    pdb_path = "data/mock/mock.pdb"
    
    try:
        # Parse the structure
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('test', pdb_path)
        sequence, residues = parse_pdb_structure(pdb_path)
        
        print(f"✓ Loaded PDB with {len(residues)} residues")
        
        # Map only the first residue
        mapping = {0: 0}
        scores = [99.9]
        unmapped_value = -1.0
        
        print(f"  Mapping: only first residue")
        print(f"  Unmapped value: {unmapped_value}")
        
        # Create temporary output file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as tmp:
            output_path = tmp.name
        
        try:
            # Inject B-factors with custom unmapped value
            inject_bfactors(structure, mapping, scores, output_path, 
                          unmapped_value=unmapped_value)
            print(f"✓ B-factors injected")
            
            # Re-parse and verify
            new_parser = PDBParser(QUIET=True)
            new_structure = new_parser.get_structure('output', output_path)
            
            new_residues = []
            for model in new_structure:
                for chain in model:
                    for residue in chain:
                        if residue.id[0] == ' ':
                            new_residues.append(residue)
                break
            
            # Check first residue (mapped)
            first_atom = list(new_residues[0])[0]
            assert abs(first_atom.get_bfactor() - 99.9) < 0.01
            print(f"✓ Mapped residue has correct B-factor (99.9)")
            
            # Check second residue (unmapped)
            second_atom = list(new_residues[1])[0]
            assert abs(second_atom.get_bfactor() - unmapped_value) < 0.01
            print(f"✓ Unmapped residue has custom value ({unmapped_value})")
            
            print("\n✓✓✓ Test 4 PASSED ✓✓✓\n")
            return True
            
        finally:
            if os.path.exists(output_path):
                os.remove(output_path)
        
    except Exception as e:
        print(f"✗ Test 4 FAILED: {e}")
        import traceback
        traceback.print_exc()
        return False


def main():
    """Run all tests."""
    print("\n" + "=" * 60)
    print("PDB INJECTOR TEST SUITE")
    print("=" * 60 + "\n")
    
    # Change to the project directory
    project_dir = os.path.dirname(os.path.abspath(__file__))
    os.chdir(project_dir)
    
    results = []
    results.append(("Test 1: Basic injection", test_basic_injection()))
    results.append(("Test 2: All mapped", test_all_residues_mapped()))
    results.append(("Test 3: None mapped", test_no_residues_mapped()))
    results.append(("Test 4: Custom unmapped", test_custom_unmapped_value()))
    
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
