
import unittest
import numpy as np
from ribostruct.core.processor import shift_nucleotide_density

class TestOffsetDirection(unittest.TestCase):
    def test_negative_offset_shifts_downstream(self):
        """
        User requirement: "-12 should show the signal 4 AAs down stream"
        This means offset = -12 should result in a shift to higher indices (right).
        
        Test with offset = -1.
        Input:  [0, 1, 0, 0, 0] (Signal at index 1)
        Expect: [0, 0, 1, 0, 0] (Signal at index 2)
        """
        density = np.array([0.0, 1.0, 0.0, 0.0, 0.0])
        offset = -1
        
        # We expect the signal to move RIGHT (downstream)
        expected = np.array([0.0, 0.0, 1.0, 0.0, 0.0])
        
        result = shift_nucleotide_density(density, offset)
        
        print(f"\nInput:    {density}")
        print(f"Offset:   {offset}")
        print(f"Result:   {result}")
        print(f"Expected: {expected}")
        
        np.testing.assert_array_equal(result, expected, 
                                      err_msg="Negative offset should shift downstream (right)!")

    def test_positive_offset_shifts_upstream(self):
        """
        If negative shifts downstream, positive should shift upstream.
        Offset = +1.
        Input:  [0, 1, 0, 0, 0]
        Expect: [1, 0, 0, 0, 0]
        """
        density = np.array([0.0, 1.0, 0.0, 0.0, 0.0])
        offset = 1
        
        expected = np.array([1.0, 0.0, 0.0, 0.0, 0.0])
        
        result = shift_nucleotide_density(density, offset)
        
        np.testing.assert_array_equal(result, expected, 
                                      err_msg="Positive offset should shift upstream (left)!")

if __name__ == '__main__':
    unittest.main()
