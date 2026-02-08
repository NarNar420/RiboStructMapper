
import os
import sys
import numpy as np
from ribostruct.cli.main import run_pipeline
from ribostruct.core.parser import parse_genomic_data, parse_ribo_density, parse_pdb_structure
from ribostruct.core.alignment import translate_sequence, align_sequences

def reproduce():
    print("Reproducing issue with real data...")
    
    pdb_path = "data/real/AF-P10664-F1-model_v6.pdb"
    fasta_path = "data/real/RPL4A_sequence.fasta"
    bedgraph_path = "data/real/GSE100626_SRX2965777.bedgraph"
    
    # Check files
    if not os.path.exists(pdb_path): print(f"Missing {pdb_path}"); return
    if not os.path.exists(fasta_path): print(f"Missing {fasta_path}"); return
    if not os.path.exists(bedgraph_path): print(f"Missing {bedgraph_path}"); return

    # 1. Parse Genomic
    print("\n[1] Parser Check")
    nucleotide_seq, coord_map, start_offset, chrom, genomic_start = parse_genomic_data(fasta_path)
    print(f"Seq Len: {len(nucleotide_seq)}")
    print(f"Start Offset: {start_offset}")
    print(f"Genomic Start: {genomic_start}")
    print(f"First 10 NT: {nucleotide_seq[:10]}")
    
    # 2. Parse Density
    end = genomic_start + len(nucleotide_seq) - 1
    print(f"\n[2] Density Check ({chrom}:{genomic_start}-{end})")
    density = parse_ribo_density(bedgraph_path, chrom, genomic_start, end)
    print(f"Density Vector Len: {len(density)}")
    print(f"Density Sum: {np.sum(density)}")
    print(f"Density Max: {np.max(density)}")
    print(f"First 10 vals: {density[:10]}")
    
    if np.sum(density) == 0:
        print("!!! DENSITY IS ALL ZERO !!!")
        
    # 3. Translate
    print("\n[3] Translation Check")
    aa_seq = translate_sequence(nucleotide_seq)
    print(f"AA Seq Len: {len(aa_seq)}")
    print(f"AA Seq: {aa_seq[:10]}...")
    
    # 4. PDB Parse
    print("\n[4] PDB Check")
    pdb_seq, residues = parse_pdb_structure(pdb_path)
    print(f"PDB Seq Len: {len(pdb_seq)}")
    print(f"PDB Seq: {pdb_seq[:10]}...")
    
    # 5. Alignment
    print("\n[5] Alignment Check")
    mapping = align_sequences(aa_seq, pdb_seq)
    print(f"Mapped residues: {len(mapping)}")
    
    # 6. Run Pipeline
    print("\n[6] Full Pipeline Run")
    try:
        run_pipeline(
            pdb_path=pdb_path,
            fasta_path=fasta_path,
            bedgraph_path=bedgraph_path,
            offsets=[0],
            output_dir="tests/output_real_debug"
        )
        print("Pipeline finished.")
    except Exception as e:
        print(f"Pipeline failed: {e}")

if __name__ == "__main__":
    reproduce()
