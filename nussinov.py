#/usr/bin/python
import numpy as np

def initialize(N, min_dist):
    matrix = np.empty((N,N))
    matrix[:] = np.nan
    matrix[range(N), range(N)] = 0
    matrix[range(1, N), range(N-1)] = 0
    return matrix

def is_paired(pair, weights = [2, 3, 1]):
    if pair in [('A', 'U'),('U', 'A')]:
        return weights[0]
    if pair in [('G', 'C'), ('C', 'G')]:
        return weights[1]
    if pair in [('G', 'U'), ('U', 'G')]:
        return weights[2]
    else:
        return 0

def fill(i, j, seq, min_dist = 3, weights = [2, 3, 1]):
    
    # Base case
    if i >= j - min_dist:
        return 0

    score_paired = fill(i + 1, j - 1, seq, min_dist, weights) + is_paired((seq[i],seq[j]), weights)
    score_i_not_paired = fill(i + 1, j, seq, min_dist, weights)
    score_j_not_paired = fill(i, j - 1, seq, min_dist, weights)
    score_all_bifurcation = [fill(i, k, seq, min_dist, weights) + fill(k + 1, j, seq, min_dist, weights) for k in range(i, j - min_dist)]
    score_bifurcation = max(score_all_bifurcation)

    return max(score_paired, score_i_not_paired, score_j_not_paired, score_bifurcation)

def main():
    WEIGHTS = [1, 1, 1]
    MIN_DIST = 0
    SEQ = 'GGGAAAUCC'
    LEN = len(SEQ)
    m = initialize(LEN, MIN_DIST)
    print(m)
    secundary_structure = []
    for k in range(LEN):
        for i in range(LEN-k):
            j = i + k
            m[i][j] = fill(i, j, SEQ, MIN_DIST, WEIGHTS)
    print(m)
    return 0


if __name__ == '__main__':
    main()