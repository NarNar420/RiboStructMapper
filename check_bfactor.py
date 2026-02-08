
import sys
from Bio.PDB import PDBParser

def check_bfactors(pdb_path):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('test', pdb_path)
    
    bfactors = []
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    bfactors.append(atom.get_bfactor())
                    
    min_b = min(bfactors)
    max_b = max(bfactors)
    avg_b = sum(bfactors) / len(bfactors)
    
    print(f"File: {pdb_path}")
    print(f"B-factors: Min={min_b:.2f}, Max={max_b:.2f}, Avg={avg_b:.2f}")
    
    if max_b == 0:
        print("FAIL: All B-factors are zero!")
    else:
        print("PASS: Non-zero B-factors found.")

if __name__ == "__main__":
    if len(sys.argv) > 1:
        check_bfactors(sys.argv[1])
    else:
        print("Usage: python check_bfactor.py <pdb_file>")
