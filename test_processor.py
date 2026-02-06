"""
Test script for the Data Processor module.

This script verifies that:
1. Density aggregation works correctly with different methods
2. Offset application shifts vectors properly
3. Boundary conditions are handled correctly (filling with 0.0)
"""

from processor import aggregate_density, apply_offsets


def test_aggregate_density():
    """Test the aggregate_density function with various methods."""
    print("="*60)
    print("Testing aggregate_density()")
    print("="*60)
    
    # Test vector: [1, 2, 3, 4, 5, 6] represents 2 codons
    test_vector = [1, 2, 3, 4, 5, 6]
    
    print(f"\nInput vector (6 nucleotides = 2 codons): {test_vector}")
    print(f"  Codon 1: {test_vector[0:3]} = [1, 2, 3]")
    print(f"  Codon 2: {test_vector[3:6]} = [4, 5, 6]")
    
    # Test mean aggregation
    print(f"\n[1] Testing method='mean':")
    result_mean = aggregate_density(test_vector, method='mean')
    expected_mean = [2.0, 5.0]  # (1+2+3)/3 = 2.0, (4+5+6)/3 = 5.0
    print(f"    Result: {result_mean}")
    print(f"    Expected: {expected_mean}")
    assert result_mean == expected_mean, f"Mean aggregation failed: {result_mean} != {expected_mean}"
    print(f"    ✓ PASSED")
    
    # Test max aggregation
    print(f"\n[2] Testing method='max':")
    result_max = aggregate_density(test_vector, method='max')
    expected_max = [3.0, 6.0]  # max(1,2,3) = 3.0, max(4,5,6) = 6.0
    print(f"    Result: {result_max}")
    print(f"    Expected: {expected_max}")
    assert result_max == expected_max, f"Max aggregation failed: {result_max} != {expected_max}"
    print(f"    ✓ PASSED")
    
    # Test sum aggregation
    print(f"\n[3] Testing method='sum':")
    result_sum = aggregate_density(test_vector, method='sum')
    expected_sum = [6.0, 15.0]  # 1+2+3 = 6.0, 4+5+6 = 15.0
    print(f"    Result: {result_sum}")
    print(f"    Expected: {expected_sum}")
    assert result_sum == expected_sum, f"Sum aggregation failed: {result_sum} != {expected_sum}"
    print(f"    ✓ PASSED")
    
    # Test median aggregation
    print(f"\n[4] Testing method='median':")
    result_median = aggregate_density(test_vector, method='median')
    expected_median = [2.0, 5.0]  # median(1,2,3) = 2.0, median(4,5,6) = 5.0
    print(f"    Result: {result_median}")
    print(f"    Expected: {expected_median}")
    assert result_median == expected_median, f"Median aggregation failed: {result_median} != {expected_median}"
    print(f"    ✓ PASSED")
    
    print(f"\n{'='*60}")
    print("All aggregation tests PASSED!")
    print(f"{'='*60}\n")


def test_apply_offsets():
    """Test the apply_offsets function with various offset values."""
    print("="*60)
    print("Testing apply_offsets()")
    print("="*60)
    
    # Start with the aggregated vector from previous test
    aa_vector = [2.0, 5.0]
    print(f"\nInput AA vector: {aa_vector}")
    print(f"  Index 0: {aa_vector[0]}")
    print(f"  Index 1: {aa_vector[1]}")
    
    # Test offset 0 (no shift)
    print(f"\n[1] Testing offset=0 (no shift):")
    result = apply_offsets(aa_vector, [0])
    print(f"    Result: {result[0]}")
    print(f"    Expected: {aa_vector}")
    assert result[0] == aa_vector, f"Offset 0 failed: {result[0]} != {aa_vector}"
    print(f"    ✓ PASSED")
    
    # Test offset -3 (shift left by 1 AA = -3 nucleotides)
    print(f"\n[2] Testing offset=-3 (shift left by 1 AA):")
    print(f"    Logic: Index 1 moves to index 0, index 0 is filled with 0.0")
    result = apply_offsets(aa_vector, [-3])
    expected = [5.0, 0.0]  # Shift left: position 1→0, position 0 gets 0.0
    print(f"    Result: {result[-3]}")
    print(f"    Expected: {expected}")
    assert result[-3] == expected, f"Offset -3 failed: {result[-3]} != {expected}"
    print(f"    ✓ PASSED")
    
    # Test offset +3 (shift right by 1 AA = +3 nucleotides)
    print(f"\n[3] Testing offset=+3 (shift right by 1 AA):")
    print(f"    Logic: Index 0 moves to index 1, index 1 is filled with 0.0")
    result = apply_offsets(aa_vector, [3])
    expected = [0.0, 2.0]  # Shift right: position 0→1, position 1 gets 0.0
    print(f"    Result: {result[3]}")
    print(f"    Expected: {expected}")
    assert result[3] == expected, f"Offset +3 failed: {result[3]} != {expected}"
    print(f"    ✓ PASSED")
    
    # Test multiple offsets at once
    print(f"\n[4] Testing multiple offsets: [0, -3, 3]:")
    result = apply_offsets(aa_vector, [0, -3, 3])
    print(f"    Offset 0: {result[0]}")
    print(f"    Offset -3: {result[-3]}")
    print(f"    Offset +3: {result[3]}")
    assert len(result) == 3, f"Expected 3 results, got {len(result)}"
    assert result[0] == [2.0, 5.0], "Offset 0 in multiple test failed"
    assert result[-3] == [5.0, 0.0], "Offset -3 in multiple test failed"
    assert result[3] == [0.0, 2.0], "Offset +3 in multiple test failed"
    print(f"    ✓ PASSED")
    
    # Test with a longer vector
    print(f"\n[5] Testing with longer vector:")
    long_vector = [1.0, 2.0, 3.0, 4.0, 5.0]
    print(f"    Input: {long_vector}")
    result = apply_offsets(long_vector, [-6])  # -6 nt = -2 AA
    expected = [3.0, 4.0, 5.0, 0.0, 0.0]  # Shift left by 2
    print(f"    Offset -6 (-2 AA): {result[-6]}")
    print(f"    Expected: {expected}")
    assert result[-6] == expected, f"Long vector offset failed: {result[-6]} != {expected}"
    print(f"    ✓ PASSED")
    
    print(f"\n{'='*60}")
    print("All offset tests PASSED!")
    print(f"{'='*60}\n")


def test_integration():
    """Test the full pipeline: raw density → aggregation → offset."""
    print("="*60)
    print("Integration Test: Full Pipeline")
    print("="*60)
    
    # Simulate raw density data for 5 amino acids (15 nucleotides)
    raw_density = [
        10, 12, 11,  # Codon 1
        20, 22, 21,  # Codon 2
        30, 32, 31,  # Codon 3
        40, 42, 41,  # Codon 4
        50, 52, 51   # Codon 5
    ]
    
    print(f"\n[1] Raw density (15 nt = 5 codons):")
    for i in range(0, 15, 3):
        codon_num = i // 3 + 1
        print(f"    Codon {codon_num}: {raw_density[i:i+3]}")
    
    # Aggregate to AA level
    print(f"\n[2] Aggregate with method='mean':")
    aa_scores = aggregate_density(raw_density, method='mean')
    print(f"    Result: {aa_scores}")
    print(f"    Length: {len(aa_scores)} amino acids")
    
    # Apply offsets
    print(f"\n[3] Apply offsets: [0, -15, -30]:")
    offsets = [0, -15, -30]  # 0, -5 AA, -10 AA
    offset_results = apply_offsets(aa_scores, offsets)
    
    for offset_nt in offsets:
        offset_aa = offset_nt // 3
        print(f"\n    Offset {offset_nt} nt ({offset_aa} AA):")
        print(f"      {offset_results[offset_nt]}")
    
    print(f"\n{'='*60}")
    print("Integration test PASSED!")
    print(f"{'='*60}\n")


def main():
    """Run all tests."""
    print("\n" + "="*60)
    print("RiboStructMapper - Data Processor Test Suite")
    print("="*60 + "\n")
    
    try:
        test_aggregate_density()
        test_apply_offsets()
        test_integration()
        
        print("\n" + "="*60)
        print("ALL TESTS PASSED! ✓")
        print("="*60 + "\n")
        
    except AssertionError as e:
        print(f"\n❌ TEST FAILED: {e}\n")
        raise


if __name__ == "__main__":
    main()
