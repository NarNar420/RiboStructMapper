#!/usr/bin/env python3
"""
Test suite for the Data Processor module.

Tests nucleotide-level offset application and aggregation to amino acids.
"""

import numpy as np
import pytest
from ribostruct.core.processor import shift_nucleotide_density, aggregate_density, process_offsets



def test_shift_nucleotide_density_no_offset():
    """Test that zero offset returns unchanged array."""
    density = np.array([10.0, 20.0, 30.0, 40.0, 50.0, 60.0])
    shifted = shift_nucleotide_density(density, 0)
    
    np.testing.assert_array_equal(shifted, density)
    print("✓ Zero offset: no change")


def test_shift_nucleotide_density_negative():
    """Test negative offset (upstream shift)."""
    # Original: [10, 20, 30, 40, 50, 60]
    # Offset -2: shift left by 2
    # Result: [0, 0, 10, 20, 30, 40]
    density = np.array([10.0, 20.0, 30.0, 40.0, 50.0, 60.0])
    shifted = shift_nucleotide_density(density, -2)
    
    expected = np.array([0.0, 0.0, 10.0, 20.0, 30.0, 40.0])
    np.testing.assert_array_equal(shifted, expected)
    print("✓ Negative offset -2: shifted upstream correctly")


def test_shift_nucleotide_density_positive():
    """Test positive offset (downstream shift)."""
    # Original: [10, 20, 30, 40, 50, 60]
    # Offset +2: shift right by 2
    # Result: [30, 40, 50, 60, 0, 0]
    density = np.array([10.0, 20.0, 30.0, 40.0, 50.0, 60.0])
    shifted = shift_nucleotide_density(density, 2)
    
    expected = np.array([30.0, 40.0, 50.0, 60.0, 0.0, 0.0])
    np.testing.assert_array_equal(shifted, expected)
    print("✓ Positive offset +2: shifted downstream correctly")


def test_shift_nucleotide_density_by_codon():
    """Test shift by one codon (3 nucleotides)."""
    # Original: [10, 10, 10, 20, 20, 20]  (2 codons)
    # Offset -3: shift left by 3 (one codon)
    # Result: [0, 0, 0, 10, 10, 10]
    density = np.array([10.0, 10.0, 10.0, 20.0, 20.0, 20.0])
    shifted = shift_nucleotide_density(density, -3)
    
    expected = np.array([0.0, 0.0, 0.0, 10.0, 10.0, 10.0])
    np.testing.assert_array_equal(shifted, expected)
    print("✓ Offset -3 (one codon): shifted correctly")


def test_aggregate_density_mean():
    """Test mean aggregation method."""
    # Two codons: [10, 20, 30] and [40, 50, 60]
    # Means: 20.0 and 50.0
    nt_density = np.array([10.0, 20.0, 30.0, 40.0, 50.0, 60.0])
    aa_density = aggregate_density(nt_density, method='mean')
    
    expected = np.array([20.0, 50.0])
    np.testing.assert_array_almost_equal(aa_density, expected)
    print("✓ Mean aggregation: correct")


def test_aggregate_density_max():
    """Test max aggregation method."""
    nt_density = np.array([10.0, 20.0, 30.0, 40.0, 50.0, 60.0])
    aa_density = aggregate_density(nt_density, method='max')
    
    expected = np.array([30.0, 60.0])
    np.testing.assert_array_almost_equal(aa_density, expected)
    print("✓ Max aggregation: correct")


def test_aggregate_density_sum():
    """Test sum aggregation method."""
    nt_density = np.array([10.0, 20.0, 30.0, 40.0, 50.0, 60.0])
    aa_density = aggregate_density(nt_density, method='sum')
    
    expected = np.array([60.0, 150.0])
    np.testing.assert_array_almost_equal(aa_density, expected)
    print("✓ Sum aggregation: correct")


def test_aggregate_density_median():
    """Test median aggregation method."""
    nt_density = np.array([10.0, 20.0, 30.0, 40.0, 50.0, 60.0])
    aa_density = aggregate_density(nt_density, method='median')
    
    expected = np.array([20.0, 50.0])
    np.testing.assert_array_almost_equal(aa_density, expected)
    print("✓ Median aggregation: correct")


def test_aggregate_density_invalid_length():
    """Test that non-divisible-by-3 length raises error."""
    nt_density = np.array([10.0, 20.0, 30.0, 40.0])  # Length 4
    
    with pytest.raises(ValueError, match="must be divisible by 3"):
        aggregate_density(nt_density)
    
    print("✓ Invalid length: error raised correctly")


def test_aggregate_density_invalid_method():
    """Test that invalid method raises error."""
    nt_density = np.array([10.0, 20.0, 30.0])
    
    with pytest.raises(ValueError, match="Unsupported aggregation method"):
        aggregate_density(nt_density, method='invalid')
    
    print("✓ Invalid method: error raised correctly")


def test_critical_offset_before_aggregation():
    """
    CRITICAL TEST: Verify offsets are applied BEFORE aggregation.
    
    This test proves we're doing nucleotide-level shifting, not AA-level.
    
    Setup:
    - 6 nucleotides: [10, 10, 10, 20, 20, 20]
    - 2 codons: Met (10,10,10) and Pro (20,20,20)
    
    Offset -1 (shift left by 1 nt):
    - Shifted NT: [10, 10, 20, 20, 20, 0]
    - New codons: (10,10,20) and (20,20,0)
    - AA means: 13.33 and 13.33
    
    If we incorrectly applied offset at AA level:
    - AA before: [10, 20]
    - Offset -1 nt = -0.33 AA (rounds to 0)
    - Would get: [10, 20] (no change!)
    
    This test MUST pass to prove correctness.
    """
    raw_density = np.array([10.0, 10.0, 10.0, 20.0, 20.0, 20.0])
    
    # Apply offset -1 (shift RIGHT by 1 nt)
    shifted_nt = shift_nucleotide_density(raw_density, -1)
    expected_shifted = np.array([0.0, 10.0, 10.0, 10.0, 20.0, 20.0])
    np.testing.assert_array_equal(shifted_nt, expected_shifted)
    
    # Aggregate to AA level
    aa_density = aggregate_density(shifted_nt, method='mean')
    expected_aa = np.array([20.0/3.0, 50.0/3.0])  # 6.67, 16.67
    np.testing.assert_array_almost_equal(aa_density, expected_aa, decimal=2)
    
    print("✓ CRITICAL: Offset -1 applied at NT level BEFORE aggregation")
    print(f"  Raw NT: {raw_density}")
    print(f"  Shifted NT: {shifted_nt}")
    print(f"  AA density: {aa_density}")


def test_process_offsets_multiple():
    """Test processing multiple offsets."""
    raw_density = np.array([10.0, 10.0, 10.0, 20.0, 20.0, 20.0])
    offsets = [0, -3]
    
    results = process_offsets(raw_density, offsets, method='mean')
    
    # Offset 0: no shift
    # Codons: (10,10,10) and (20,20,20)
    # Means: 10.0 and 20.0
    expected_0 = np.array([10.0, 20.0])
    np.testing.assert_array_almost_equal(results[0], expected_0)
    
    # Offset -3: shift RIGHT by one codon
    # Shifted: [0, 0, 0, 10, 10, 10]
    # Codons: (0,0,0) and (10,10,10)
    # Means: 0.0 and 10.0
    expected_minus3 = np.array([0.0, 10.0])
    np.testing.assert_array_almost_equal(results[-3], expected_minus3)
    
    print("✓ Multiple offsets processed correctly")


def test_process_offsets_different_methods():
    """Test that process_offsets works with different aggregation methods."""
    raw_density = np.array([10.0, 20.0, 30.0, 40.0, 50.0, 60.0])
    
    # Test with max
    results_max = process_offsets(raw_density, [0], method='max')
    expected = np.array([30.0, 60.0])
    np.testing.assert_array_almost_equal(results_max[0], expected)
    
    # Test with sum
    results_sum = process_offsets(raw_density, [0], method='sum')
    expected = np.array([60.0, 150.0])
    np.testing.assert_array_almost_equal(results_sum[0], expected)
    
    print("✓ process_offsets works with different methods")


def test_edge_case_large_offset():
    """Test offset larger than array length."""
    density = np.array([10.0, 20.0, 30.0])
    
    # Offset -10 (larger than array)
    shifted = shift_nucleotide_density(density, -10)
    expected = np.array([0.0, 0.0, 0.0])  # All zeros
    np.testing.assert_array_equal(shifted, expected)
    
    print("✓ Large offset: all zeros (correct)")


def test_real_world_example():
    """
    Test with realistic ribosome profiling data.
    
    Simulates P-site offset of -12 nucleotides.
    """
    # 15 nucleotides (5 codons)
    raw_density = np.array([
        5.0, 10.0, 15.0,   # Codon 1
        20.0, 25.0, 30.0,  # Codon 2
        35.0, 40.0, 45.0,  # Codon 3
        50.0, 55.0, 60.0,  # Codon 4
        65.0, 70.0, 75.0   # Codon 5
    ])
    
    # Apply P-site offset (-12 nt = -4 codons)
    results = process_offsets(raw_density, [0, -12], method='mean')
    
    # Offset 0: means of each codon
    expected_0 = np.array([10.0, 25.0, 40.0, 55.0, 70.0])
    np.testing.assert_array_almost_equal(results[0], expected_0)
    
    # Offset -12: shift RIGHT by 12
    # Shifted: [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 10, 15]
    # Codons: (0,0,0), (0,0,0), (0,0,0), (0,0,0), (5,10,15)
    # Means: 0.0, 0.0, 0.0, 0.0, 10.0
    expected_minus12 = np.array([0.0, 0.0, 0.0, 0.0, 10.0])
    np.testing.assert_array_almost_equal(results[-12], expected_minus12)
    
    print("✓ Real-world P-site offset (-12): correct")


def main():
    """Run all tests."""
    print("=" * 70)
    print("PROCESSOR MODULE TEST SUITE (NEW LOGIC)")
    print("=" * 70)
    print()
    
    tests = [
        ("Zero offset", test_shift_nucleotide_density_no_offset),
        ("Negative offset", test_shift_nucleotide_density_negative),
        ("Positive offset", test_shift_nucleotide_density_positive),
        ("Shift by codon", test_shift_nucleotide_density_by_codon),
        ("Mean aggregation", test_aggregate_density_mean),
        ("Max aggregation", test_aggregate_density_max),
        ("Sum aggregation", test_aggregate_density_sum),
        ("Median aggregation", test_aggregate_density_median),
        ("Invalid length error", test_aggregate_density_invalid_length),
        ("Invalid method error", test_aggregate_density_invalid_method),
        ("⭐ CRITICAL: NT-level offset", test_critical_offset_before_aggregation),
        ("Multiple offsets", test_process_offsets_multiple),
        ("Different methods", test_process_offsets_different_methods),
        ("Large offset edge case", test_edge_case_large_offset),
        ("Real-world P-site", test_real_world_example),
    ]
    
    passed = 0
    failed = 0
    
    for name, test_func in tests:
        try:
            test_func()
            passed += 1
        except Exception as e:
            print(f"✗ FAILED: {name}")
            print(f"  Error: {e}")
            import traceback
            traceback.print_exc()
            failed += 1
        print()
    
    print("=" * 70)
    print(f"RESULTS: {passed}/{len(tests)} tests passed")
    if failed > 0:
        print(f"❌ {failed} tests FAILED")
    else:
        print("✅ ALL TESTS PASSED")
    print("=" * 70)
    
    return failed == 0


if __name__ == "__main__":
    import sys
    success = main()
    sys.exit(0 if success else 1)
