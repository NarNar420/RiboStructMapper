#!/usr/bin/env python3
"""
Test CLI with FASTA-only mode (no GTF).
"""
from ribostruct.cli.main import run_pipeline

# Run pipeline with CDS FASTA only (no GTF)
output_files = run_pipeline(
    pdb_path="data/mock/mock.pdb",
    fasta_path="data/mock/cds_only.fasta",
    gtf_path=None,  # No GTF!
    bedgraph_path="data/mock/mock.bedgraph",
    offsets=[0],
    aggregation_method='mean',
    output_dir='.'
)

print("\n✅ FASTA-only mode test successful!")
print(f"Generated: {list(output_files.values())}")
