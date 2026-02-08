"""
Test sanitize_sequence function for frame-error handling.
"""
import pytest
from ribostruct.core.parser import sanitize_sequence


def test_sanitize_dirty_sequence():
    """Test sanitization with leading garbage and incomplete codons."""
    # "CC" is garbage (offset=2)
    # "ATG" is start codon
    # "TTC" is codon
    # "GGT" is codon
    # "TA" is 2 nts (remainder - should be removed)
    dirty_seq = "CCATGTTCGGTTA"
    
    clean_seq, offset = sanitize_sequence(dirty_seq)
    
    # Should remove "CC" from start and "TA" from end
    assert clean_seq == "ATGTTCGGT", f"Expected 'ATGTTCGGT', got '{clean_seq}'"
    assert len(clean_seq) == 9, f"Expected length 9, got {len(clean_seq)}"
    assert len(clean_seq) % 3 == 0, "Sequence must be divisible by 3"
    assert offset == 2, f"Expected offset 2, got {offset}"


def test_sanitize_clean_sequence():
    """Test sanitization with already-clean sequence."""
    clean_seq_input = "ATGAAACCC"
    
    clean_seq, offset = sanitize_sequence(clean_seq_input)
    
    # Should remain unchanged
    assert clean_seq == "ATGAAACCC"
    assert offset == 0
    assert len(clean_seq) % 3 == 0


def test_sanitize_no_atg():
    """Test sanitization when no ATG is found."""
    no_atg_seq = "CCCGGGTTTAGC"  # 12 nts, divisible by 3, but no ATG
    
    clean_seq, offset = sanitize_sequence(no_atg_seq)
    
    # Should keep sequence as-is (since divisible by 3)
    assert clean_seq == "CCCGGGTTTAGC"
    assert offset == 0
    assert len(clean_seq) % 3 == 0


def test_sanitize_no_atg_with_remainder():
    """Test sanitization when no ATG and has remainder."""
    no_atg_seq = "CCCGGGTTA"  # 9 nts + "GC" = 11 total
    no_atg_seq = "CCCGGGTTAG C"  # Add 2 more for remainder
    
    clean_seq, offset = sanitize_sequence(no_atg_seq.replace(" ", ""))
    
    # Should remove 2 from end
    assert clean_seq == "CCCGGGTTA"
    assert offset == 0
    assert len(clean_seq) % 3 == 0


def test_sanitize_atg_at_end():
    """Test when ATG appears late in sequence."""
    late_atg = "CCCCCCCCCCATGAAA"  # 10 garbage + ATG + AAA
    
    clean_seq, offset = sanitize_sequence(late_atg)
    
    # Should trim 10 from start, keep ATGAAA (6 nts)
    assert clean_seq == "ATGAAA"
    assert offset == 10
    assert len(clean_seq) % 3 == 0


def test_sanitize_lowercase_input():
    """Test with lowercase input."""
    lower_seq = "ccatgttcggtta"
    
    clean_seq, offset = sanitize_sequence(lower_seq)
    
    # Should uppercase and sanitize
    assert clean_seq == "ATGTTCGGT"
    assert offset == 2
    assert clean_seq.isupper()


def test_sanitize_multiple_atg():
    """Test with multiple ATG codons - should use first."""
    multi_atg = "CCATGAAATGATG"
    
    clean_seq, offset = sanitize_sequence(multi_atg)
    
    # Should start at first ATG
    assert clean_seq.startswith("ATG")
    assert offset == 2
    assert len(clean_seq) % 3 == 0
