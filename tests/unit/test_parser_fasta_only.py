"""
Tests for FASTA-only mode (without GTF).
"""
import pytest
from ribostruct.core.parser import parse_genomic_data


def test_parse_genomic_data_fasta_only():
    """Test parsing CDS FASTA without GTF."""
    fasta_path = "data/mock/cds_only.fasta"
    
    # Parse without GTF
    nucleotide_seq, coord_map = parse_genomic_data(
        fasta_path=fasta_path,
        gtf_path=None  # No GTF
    )
    
    # Verify sequence extracted
    assert len(nucleotide_seq) == 45
    assert nucleotide_seq.startswith("ATG")
    
    # Verify coordinate map is simple 1-based
    assert coord_map[0] == ("MOCK_CDS", 1)
    assert coord_map[1] == ("MOCK_CDS", 2)
    assert coord_map[-1] == ("MOCK_CDS", 45)
    assert len(coord_map) == 45


def test_parse_genomic_data_fasta_only_with_gene_id():
    """Test FASTA-only mode with gene_id to select specific sequence."""
    fasta_path = "data/mock/cds_only.fasta"
    
    # Parse with gene_id (which acts as sequence ID in FASTA-only mode)
    nucleotide_seq, coord_map = parse_genomic_data(
        fasta_path=fasta_path,
        gtf_path=None,
        gene_id="MOCK_CDS"
    )
    
    # Should work the same
    assert len(nucleotide_seq) == 45
    assert coord_map[0][0] == "MOCK_CDS"


def test_parse_genomic_data_fasta_only_sequence_not_found():
    """Test FASTA-only mode with non-existent gene_id."""
    fasta_path = "data/mock/cds_only.fasta"
    
    with pytest.raises(ValueError, match="Sequence ID 'NONEXISTENT' not found"):
        parse_genomic_data(
            fasta_path=fasta_path,
            gtf_path=None,
            gene_id="NONEXISTENT"
        )
