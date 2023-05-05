#!/usr/bin/python
import unittest
import numpy as np
from nussinov import initialize, is_paired, fill, traceback, pairs_to_structure, nussinov
import collections

class TestNussinov(unittest.TestCase):
    def test_initialize(self):
        N = 4
        min_dist = 3
        matrix = initialize(N, min_dist)
        expected = np.array([
            [0, np.nan, np.nan, np.nan],
            [0, 0, np.nan, np.nan],
            [np.nan, 0, 0, np.nan],
            [np.nan, np.nan, 0, 0]
        ])
        np.testing.assert_equal(matrix, expected)


    def test_is_paired(self):
        self.assertEqual(is_paired(('A', 'U')), 2)
        self.assertEqual(is_paired(('G', 'C')), 3)
        self.assertEqual(is_paired(('G', 'U')), 1)
        self.assertEqual(is_paired(('A', 'A')), 0)

    def test_fill(self):
        seq = 'UGGGGUU'
        length = len(seq)
        m = initialize(length, 3)
        for k in range(length):
            for i in range(length - k):
                j = i + k
                m[i][j] = fill(i, j, seq)
        self.assertEqual(m[0][6], 1)
        self.assertEqual(m[1][5], 1)

    def test_traceback(self):
        seq = "GGAAGCC"
        weights = [2, 3, 1]
        length = len(seq)
        min_dist = 3
        matrix = initialize(length, min_dist)
        for k in range(length):
            for i in range(length-k):
                j = i + k
                matrix[i][j] = fill(i, j, seq, min_dist, weights)
        pairs = traceback(matrix, seq, weights)
        expected_pairs = [(0, 6), (1, 5)]
        if not collections.Counter(pairs) == collections.Counter(expected_pairs):
            print("Pairs:", pairs)
            print("Expected pairs:", expected_pairs)
        self.assertCountEqual(pairs, expected_pairs)


    def test_pairs_to_structure(self):
        pairs = [(0, 6), (1, 5)]
        seq = 'GGAAGCC'
        struct = pairs_to_structure(pairs, seq)
        self.assertEqual(struct, '((...))')

    def test_nussinov(self):
        seq = 'UGGGGUUUAAGGCCCC'
        weights = [2, 3, 1]
        min_dist = 3
        seq_id = 'test'
        nussinov(seq, weights, min_dist, seq_id)
        with open('test.out', 'r') as f:
            output = f.read()
        expected_output = '>test\nUGGGGUUUAAGGCCCC\n.(((((...)..))))'
        self.assertEqual(output.strip(), expected_output.strip())

if __name__ == '__main__':
    unittest.main()