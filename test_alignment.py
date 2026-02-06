"""
Test script for the Alignment Engine.

This script verifies that:
1. Nucleotide sequences can be translated to amino acids
2. Pairwise alignment works correctly
3. The mapping dictionary properly handles gaps (missing residues in PDB)
"""

from Bio import SeqIO
from Bio.PDB import PDBParser
from alignment import translate_sequence, align_sequences


def extract_pdb_sequence(pdb_path: str, chain_id: str = 'A') -> str:
    """
    Extract the amino acid sequence from a PDB file.
    
    Args:
        pdb_path: Path to PDB file
        chain_id: Chain identifier (default 'A')
        
    Returns:
        Amino acid sequence as a string
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('protein', pdb_path)
    
    # Standard three-letter to one-letter amino acid code
    aa_dict = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
        'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
        'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
    }
    
    sequence = []
    for model in structure:
        for chain in model:
            if chain.id == chain_id:
                for residue in chain:
                    # Only consider standard amino acids
                    if residue.id[0] == ' ':  # Not a heteroatom
                        res_name = residue.get_resname()
                        if res_name in aa_dict:
                            sequence.append(aa_dict[res_name])
    
    return ''.join(sequence)


def main():
    print("="*60)
    print("RiboStructMapper - Alignment Engine Test")
    print("="*60)
    
    # 1. Load and translate the genomic sequence
    print("\n[1] Loading mock genomic sequence...")
    fasta_path = "mock_data/mock.fasta"
    
    fasta_records = list(SeqIO.parse(fasta_path, "fasta"))
    nucleotide_seq = str(fasta_records[0].seq)
    print(f"    Nucleotide sequence: {nucleotide_seq}")
    print(f"    Length: {len(nucleotide_seq)} bp")
    
    genomic_aa = translate_sequence(nucleotide_seq)
    print(f"\n[2] Translated to amino acids:")
    print(f"    Genomic AA sequence: {genomic_aa}")
    print(f"    Length: {len(genomic_aa)} residues")
    
    # 2. Extract PDB sequence
    print(f"\n[3] Extracting sequence from PDB...")
    pdb_path = "mock_data/mock.pdb"
    pdb_aa_original = extract_pdb_sequence(pdb_path)
    print(f"    Original PDB AA sequence: {pdb_aa_original}")
    print(f"    Length: {len(pdb_aa_original)} residues")
    
    # 3. Introduce a mutation (delete one residue in the middle)
    print(f"\n[4] Introducing mutation (deleting residue at position 7)...")
    mutation_pos = 7  # 0-based index, so this is the 8th residue
    pdb_aa_mutated = pdb_aa_original[:mutation_pos] + pdb_aa_original[mutation_pos + 1:]
    print(f"    Deleted residue: {pdb_aa_original[mutation_pos]} at position {mutation_pos}")
    print(f"    Mutated PDB AA sequence: {pdb_aa_mutated}")
    print(f"    New length: {len(pdb_aa_mutated)} residues")
    
    # 4. Perform alignment
    print(f"\n[5] Performing pairwise alignment...")
    alignment_score, mapping_dict = align_sequences(genomic_aa, pdb_aa_mutated)
    
    print(f"\n{'='*60}")
    print(f"ALIGNMENT RESULTS")
    print(f"{'='*60}")
    print(f"Alignment Score: {alignment_score}")
    print(f"Total mappings: {len(mapping_dict)}")
    print(f"\nFirst 5 mappings (Genomic Index -> PDB Index):")
    
    for i, (genomic_idx, pdb_idx) in enumerate(sorted(mapping_dict.items())[:5]):
        genomic_residue = genomic_aa[genomic_idx]
        pdb_residue = pdb_aa_mutated[pdb_idx]
        print(f"  {genomic_idx:2d} -> {pdb_idx:2d}  |  {genomic_residue} -> {pdb_residue}")
    
    # 5. Verify gap handling
    print(f"\n[6] Gap handling verification:")
    print(f"    Expected: Position {mutation_pos} in genomic should NOT map to PDB")
    
    if mutation_pos in mapping_dict:
        print(f"    ❌ FAILED: Position {mutation_pos} is in mapping (should be skipped)")
    else:
        print(f"    ✓ PASSED: Position {mutation_pos} correctly omitted from mapping")
    
    # Check that positions after the gap are shifted correctly
    print(f"\n    Checking shift after gap:")
    for genomic_idx in range(mutation_pos + 1, min(mutation_pos + 4, len(genomic_aa))):
        if genomic_idx in mapping_dict:
            expected_pdb_idx = genomic_idx - 1  # Should be shifted by 1 due to deletion
            actual_pdb_idx = mapping_dict[genomic_idx]
            match_status = "✓" if actual_pdb_idx == expected_pdb_idx else "❌"
            print(f"    {match_status} Genomic[{genomic_idx}] -> PDB[{actual_pdb_idx}] (expected {expected_pdb_idx})")
    
    print(f"\n{'='*60}")
    print("Test completed!")
    print(f"{'='*60}\n")


if __name__ == "__main__":
    main()
