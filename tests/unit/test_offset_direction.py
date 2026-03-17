
import unittest
import numpy as np
from ribostruct.core.processor import shift_nucleotide_density

class TestOffsetDirection(unittest.TestCase):
    def test_positive_offset_shifts_downstream(self):
        """
        Verify that any positive offset shifts density to HIGHER indices
        (downstream, towards the 3' end / already-translated sequence).

        This is the biologically correct direction for P-site offset correction.

        Test with offset = +1.
        Input:  [0, 1, 0, 0, 0] (Signal at index 1)
        Expect: [0, 0, 1, 0, 0] (Signal moves to index 2)
        """
        density = np.array([0.0, 1.0, 0.0, 0.0, 0.0])
        offset = 1

        # We expect the signal to move RIGHT (downstream, higher index)
        expected = np.array([0.0, 0.0, 1.0, 0.0, 0.0])

        result = shift_nucleotide_density(density, offset)

        print(f"\nInput:    {density}")
        print(f"Offset:   +{offset}")
        print(f"Result:   {result}")
        print(f"Expected: {expected}")

        np.testing.assert_array_equal(result, expected,
                                      err_msg="Positive offset should shift downstream (right)!")

    def test_negative_input_treated_as_absolute_downstream(self):
        """
        Verify that negative inputs are treated via abs() and still shift
        downstream (the one-way property of the offset system).

        The old system required negative numbers, e.g. -12.
        The new system converts them to abs(-12) = 12, and shifts downstream.

        Test with offset = -1 (should behave identically to offset = +1).
        Input:  [0, 1, 0, 0, 0]
        Expect: [0, 0, 1, 0, 0]  (same as +1)
        """
        density = np.array([0.0, 1.0, 0.0, 0.0, 0.0])
        offset = -1

        # abs(-1) = 1, so shift downstream by 1
        expected = np.array([0.0, 0.0, 1.0, 0.0, 0.0])

        result = shift_nucleotide_density(density, offset)

        print(f"\nInput:    {density}")
        print(f"Offset:   {offset} (-> abs = 1, downstream)")
        print(f"Result:   {result}")
        print(f"Expected: {expected}")

        np.testing.assert_array_equal(result, expected,
                                      err_msg="Negative offset should be treated as abs() and shift downstream!")

    def test_psite_offset_12(self):
        """
        Real-world example: P-site offset of 12 nucleotides.

        [1 at index 12] -> after offset 12 -> signal appears at index 24
        """
        density = np.zeros(30)
        density[12] = 1.0  # Signal at nt 12

        result = shift_nucleotide_density(density, 12)

        # Signal should move to index 24 (12 + 12)
        expected = np.zeros(30)
        expected[24] = 1.0

        print(f"\nP-site offset 12 test:")
        print(f"  Non-zero at input: index 12")
        print(f"  Non-zero at result: index {np.nonzero(result)[0]}")

        np.testing.assert_array_equal(result, expected,
                                      err_msg="P-site offset 12 should shift signal from index 12 to index 24!")

if __name__ == '__main__':
    unittest.main()
