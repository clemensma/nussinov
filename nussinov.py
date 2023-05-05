#!/usr/bin/python
import numpy as np

def initialize(N, min_dist):
    """
    Initializes an NxN matrix for Nussinov algorithm.
    
    Args:
        N (int): Length of the RNA sequence.
        min_dist (int): Minimum distance allowed between paired bases.
        
    Returns:
        np.ndarray: Initialized NxN matrix.
    """
    matrix = np.empty((N,N))
    matrix[:] = np.nan
    matrix[range(N), range(N)] = 0
    matrix[range(1, N), range(N-1)] = 0
    return matrix

def is_paired(pair, weights = [2, 3, 1]):
    """
    Determines if a base pair is paired and returns the corresponding weight.
    
    Args:
        pair (tuple): Tuple containing two bases to check for pairing.
        weights (list): List of weights for each type of pair (default: [2, 3, 1]).
        
    Returns:
        int: The weight of the pair if paired, 0 otherwise.
    """
    if pair in [('A', 'U'),('U', 'A')]:
        return weights[0]
    if pair in [('G', 'C'), ('C', 'G')]:
        return weights[1]
    if pair in [('G', 'U'), ('U', 'G')]:
        return weights[2]
    else:
        return 0

def fill(i, j, seq, min_dist = 3, weights = [2, 3, 1]):
    """
    Recursive function to fill the NxN matrix for Nussinov algorithm.
    
    Args:
        i (int): Row index in the matrix.
        j (int): Column index in the matrix.
        seq (str): RNA sequence.
        min_dist (int): Minimum distance allowed between paired bases (default: 3).
        weights (list): List of weights for each type of pair (default: [2, 3, 1]).
        
    Returns:
        int: The maximum number of base pairings from i to j.
    """
    # Base case
    if i >= j - min_dist:
        return 0

    score_paired = fill(i + 1, j - 1, seq, min_dist, weights) + is_paired((seq[i],seq[j]), weights)
    score_i_not_paired = fill(i + 1, j, seq, min_dist, weights)
    score_j_not_paired = fill(i, j - 1, seq, min_dist, weights)
    score_all_bifurcation = [fill(i, k, seq, min_dist, weights) + fill(k + 1, j, seq, min_dist, weights) for k in range(i, j - min_dist)]
    score_bifurcation = max(score_all_bifurcation)

    return max(score_paired, score_i_not_paired, score_j_not_paired, score_bifurcation)

def traceback(matrix, seq, weights):
    """
    Traceback through the filled matrix to obtain the list of base pairings.
    
    Args:
        matrix (np.ndarray): Filled NxN matrix from Nussinov algorithm.
        seq (str): RNA sequence.
        weights (list): List of weights for each type of pair.
        
    Returns:
        list: List of tuples containing the paired bases' indices.
    """
    stack = [(0,len(seq)-1)]
    pairs = []
    while(stack):
        i, j = stack.pop()
        if i >= j:
            continue
        elif int(matrix[i+1][j]) == int(matrix[i][j]):
            stack.append((i+1, j))
        elif int(matrix[i][j-1]) == int(matrix[i][j]):
            stack.append((i, j-1))
        elif int(matrix[i+1][j-1]) + int(is_paired((seq[i], seq[j]), weights)) == int(matrix[i][j]):
            pairs.append((i, j))
            stack.append((i+1, j-1))
        else:
            for k in range(i, j):
                if int(matrix[i][k-1]) + int(matrix[k+1][j-1]) + 1 == int(matrix[i][j]):
                    stack.append((k+1, j))
                    stack.append((i, k))
                    break
    return pairs

def pairs_to_structure(pairs, seq):
    """
    Converts the list of base pairings to dot-bracket notation.
    
    Args:
        pairs (list): List of tuples containing the paired bases' indices.
        seq (str): RNA sequence.
        
    Returns:
        str: RNA secondary structure in dot-bracket notation.
    """
    dots = ['.' for _ in seq]
    for pair in pairs:
        i, j = pair
        dots[i], dots[j] = '(', ')'
    return ''.join(dots)

def nussinov(seq, weights, min_dist, seq_id):
    """
    Executes the Nussinov algorithm to predict RNA secondary structure.
    
    Args:
        seq (str): RNA sequence.
        weights (list): List of weights for each type of pair.
        min_dist (int): Minimum distance allowed between paired bases.
        seq_id (str): Sequence identifier for output file.
        
    Returns:
        int: 0 when execution is successful.
    """
    length = len(seq)
    m = initialize(length, min_dist)
    for k in range(length):
        for i in range(length-k):
            j = i + k
            m[i][j] = fill(i, j, seq, min_dist, weights)
    pairs = traceback(m, seq, weights)
    struct = pairs_to_structure(pairs, seq)
    with open(seq_id + '.out', 'w') as f:
        f.write('>' + seq_id + '\n' + seq + '\n' + struct)
    
    return 0

def main():
    """
    Main function to run the Nussinov algorithm on a sample RNA sequence.
    """
    WEIGHTS = [2, 3, 1]
    MIN_DIST = 3
    SEQ = 'UGGGGUUUAAGGCCCC'
    SEQ_ID = 'test'
    nussinov(SEQ, WEIGHTS, MIN_DIST, SEQ_ID)

    return 0

if __name__ == '__main__':
    main()
