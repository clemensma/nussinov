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

def traceback(matrix, seq, weights):
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
    dots = ['.' for _ in seq]
    for pair in pairs:
        i, j = pair
        dots[i], dots[j] = '(', ')'
    return ''.join(dots)

def nussinov(seq, weights, min_dist, seq_id):
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
    WEIGHTS = [2, 3, 1]
    MIN_DIST = 3
    SEQ = 'UGGGGUUUAAGGCCCC'
    SEQ_ID = 'test'
    nussinov(SEQ, WEIGHTS, MIN_DIST, SEQ_ID)

    return 0


if __name__ == '__main__':
    main()