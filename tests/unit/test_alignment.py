"""
QA Test script for the Alignment Engine: Robustness Evaluation

This script verifies that:
1. Nucleotide sequences can be translated to amino acids
2. Pairwise semi-global alignment works correctly for various realistic PDB states:
   - Scenario A: Perfect Global Match
   - Scenario B: Missing N-Terminus in PDB (Unmodeled start)
   - Scenario C: Missing C-Terminus in PDB (Unmodeled end)
   - Scenario D: Missing Internal Loop in PDB
   - Scenario E: Single Residue Deletion / Mutation
"""
import sys
import unittest
from ribostruct.core.alignment import translate_sequence, align_sequences

class TestAlignmentEngineQA(unittest.TestCase):
    
    def test_translate_sequence(self):
        nucl = "ATGAAACCC"
        self.assertEqual(translate_sequence(nucl), "MKP")
        
    def test_scenario_a_perfect_match(self):
        genomic = "MKTIIALSY"
        pdb = "MKTIIALSY"
        score, mapping = align_sequences(genomic, pdb)
        
        # Every index should map 1:1
        self.assertEqual(len(mapping), 9)
        for i in range(9):
            self.assertEqual(mapping[i], i)
            
    def test_scenario_b_missing_n_terminus(self):
        # The PDB structure is missing the first 3 residues (MKT)
        genomic = "MKTIIALSY"
        pdb = "IIALSY"
        
        score, mapping = align_sequences(genomic, pdb)
        
        # The first 3 genomic residues shouldn't map
        self.assertNotIn(0, mapping)
        self.assertNotIn(1, mapping)
        self.assertNotIn(2, mapping)
        
        # Genomic index 3 ('I') should map to PDB index 0 ('I')
        self.assertEqual(mapping[3], 0)
        self.assertEqual(mapping[8], 5)
        self.assertEqual(len(mapping), 6)
        
    def test_scenario_c_missing_c_terminus(self):
        # The PDB structure is missing the last 3 residues (LSY)
        genomic = "MKTIIALSY"
        pdb = "MKTIIA"
        
        score, mapping = align_sequences(genomic, pdb)
        
        # Last 3 genomic residues shouldn't map
        self.assertNotIn(6, mapping)
        self.assertNotIn(7, mapping)
        self.assertNotIn(8, mapping)
        
        # First 6 should map perfectly
        self.assertEqual(mapping[0], 0)
        self.assertEqual(mapping[5], 5)
        self.assertEqual(len(mapping), 6)
        
    def test_scenario_d_missing_internal_loop(self):
        # The PDB structure is missing a 3-residue loop in the middle (IIA)
        genomic = "MKTIIALSY"
        pdb = "MKTLSY"
        
        score, mapping = align_sequences(genomic, pdb)
        
        # The internal loop indices (3, 4, 5) should be skipped
        self.assertNotIn(3, mapping)
        self.assertNotIn(4, mapping)
        self.assertNotIn(5, mapping)
        
        # Before loop
        self.assertEqual(mapping[0], 0)
        self.assertEqual(mapping[2], 2)
        
        # After loop
        # Genomic index 6 ('L') should map to PDB index 3 ('L')
        self.assertEqual(mapping[6], 3)
        self.assertEqual(mapping[8], 5)
        self.assertEqual(len(mapping), 6)
        
    def test_scenario_e_single_mutation_deletion(self):
        # Similar to the original test script: PDB is missing exactly one residue (I at index 4) due to mutation
        genomic = "MKTIIALSY"
        pdb = "MKTIALSY"
        
        score, mapping = align_sequences(genomic, pdb)
        
        # Wait, the algorithm will naturally try to align the gap correctly.
        # It's an ambiguous alignment whether it's the first or second I missing.
        # But we know total mapped should be 8
        self.assertEqual(len(mapping), 8)
        
        # Check mapping of end of chain to ensure correct shift
        self.assertEqual(mapping[8], 7) # Genomic 'Y' (8) to PDB 'Y' (7)
        

if __name__ == '__main__':
    unittest.main(verbosity=2)
