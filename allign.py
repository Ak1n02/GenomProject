def read_sequences_from_file(filename):
    with open(filename, 'r') as file:
        lines = [line.strip() for line in file if line.strip()]
    if len(lines) < 2:
        raise ValueError("File must contain at least two sequences.")
    return lines[0], lines[1]


def affine_gap_alignment(seq1, seq2, match=3, mismatch=-2, gap_open=-3, gap_extend=-1):
    len1, len2 = len(seq1), len(seq2)

    # Initialize matrices
    M = [[float('-inf')] * (len2 + 1) for _ in range(len1 + 1)]
    X = [[float('-inf')] * (len2 + 1) for _ in range(len1 + 1)]
    Y = [[float('-inf')] * (len2 + 1) for _ in range(len1 + 1)]
    backtrack = [[None] * (len2 + 1) for _ in range(len1 + 1)]

    # Initialization
    M[0][0] = 0
    for i in range(1, len1 + 1):
        X[i][0] = gap_open + (i - 1) * gap_extend
        M[i][0] = X[i][0]
    for j in range(1, len2 + 1):
        Y[0][j] = gap_open + (j - 1) * gap_extend
        M[0][j] = Y[0][j]

    # Fill matrices
    for i in range(1, len1 + 1):
        for j in range(1, len2 + 1):
            a, b = seq1[i - 1], seq2[j - 1]
            score = match if a == b else mismatch

            X[i][j] = max(X[i - 1][j] + gap_extend, M[i - 1][j] + gap_open + gap_extend)
            Y[i][j] = max(Y[i][j - 1] + gap_extend, M[i][j - 1] + gap_open + gap_extend)
            M[i][j] = max(M[i - 1][j - 1] + score, X[i][j], Y[i][j])

            if M[i][j] == M[i - 1][j - 1] + score:
                backtrack[i][j] = 'D'
            elif M[i][j] == X[i][j]:
                backtrack[i][j] = 'U'
            else:
                backtrack[i][j] = 'L'

    # Traceback
    i, j = len1, len2
    aligned_seq1 = []
    aligned_seq2 = []

    while i > 0 or j > 0:
        if backtrack[i][j] == 'D':
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append(seq2[j - 1])
            i -= 1
            j -= 1
        elif backtrack[i][j] == 'U':
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append('-')
            i -= 1
        elif backtrack[i][j] == 'L':
            aligned_seq1.append('-')
            aligned_seq2.append(seq2[j - 1])
            j -= 1
        else:
            break

    return ''.join(reversed(aligned_seq1)), ''.join(reversed(aligned_seq2)), M[len1][len2]


if __name__ == "__main__":
    import sys

    if len(sys.argv) < 2:
        print("Usage: allign.py <filename.seq>")
        sys.exit(1)

    filename = sys.argv[1]

    try:
        seq1, seq2 = read_sequences_from_file(filename)
        aligned1, aligned2, score = affine_gap_alignment(seq1, seq2)
        print("Aligned Sequence 1:\n", aligned1)
        print("Aligned Sequence 2:\n", aligned2)
        print("Alignment Score:", score)
    except Exception as e:
        print("Error:", str(e))
