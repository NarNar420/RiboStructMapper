"""
Tests for FASTA-only mode (without GTF).
"""
import pytest
import os
from ribostruct.core.parser import parse_genomic_data


def test_parse_genomic_data_fasta_only():
    """Test parsing CDS FASTA with coordinate header."""
    fasta_path = "data/mock/cds_only.fasta"
    
    # Parse (no extra args needed)
    nucleotide_seq, coord_map, start_offset, chrom, start = parse_genomic_data(
        fasta_path=fasta_path
    )
    
    # Verify sequence extracted
    assert len(nucleotide_seq) == 45
    assert nucleotide_seq.startswith("ATG")
    
    # Verify metadata extracted from header >MOCK_CDS mock_chrom:1..45
    assert chrom == "mock_chrom"
    assert start == 1
    
    # Verify coordinate map
    assert coord_map[0] == ("mock_chrom", 1)
    assert coord_map[-1] == ("mock_chrom", 45)
    assert len(coord_map) == 45
    
    # Verify no offset (clean sequence)
    assert start_offset == 0


def test_parse_genomic_data_invalid_header():
    """Test that invalid header format raises ValueError."""
    # Create temporary file with bad header
    bad_fasta = "data/mock/bad_header.fasta"
    with open(bad_fasta, "w") as f:
        f.write(">BAD_HEADER just_some_text\nATGAAA")
        
    try:
        with pytest.raises(ValueError, match="does not match required format"):
            parse_genomic_data(bad_fasta)
    finally:
        if os.path.exists(bad_fasta):
            os.remove(bad_fasta)

