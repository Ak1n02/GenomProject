# Function to read sequences from a file
def read_sequences_from_file(filename):
    """
    Reads two sequences from the given file.

    Args:
        filename (str): The name of the file containing the sequences.

    Returns:
        tuple: A tuple containing two sequences (seq1, seq2).

    Raises:
        ValueError: If the file contains less than two sequences.
    """
    with open(filename, 'r') as file:
        # Read and strip lines, ignoring blank ones
        lines = [line.strip() for line in file if line.strip()]

    # Ensure there are at least two sequences in the file
    if len(lines) < 2:
        raise ValueError("File must contain at least two sequences.")

    return lines[0], lines[1]


# Function to perform affine gap alignment
def affine_gap_alignment(seq1, seq2, match=3, mismatch=-2, gap_open=-3, gap_extend=-1):
    """
    Performs needleman-wunsch affine gap alignment on two sequences.

    Args:
        seq1 (str): The first sequence.
        seq2 (str): The second sequence.
        match (int): Score for a match. Default is 3.
        mismatch (int): Penalty for a mismatch. Default is -2.
        gap_open (int): Penalty for opening a gap. Default is -3.
        gap_extend (int): Penalty for extending a gap. Default is -1.

    Returns:
        tuple: Aligned sequences (aligned_seq1, aligned_seq2) and the alignment score.
    """
    # Get lengths of the sequences
    len1, len2 = len(seq1), len(seq2)

    # Initialize scoring matrices for affine gap penalties
    # M - Main scoring matrix
    M = [[float('-inf')] * (len2 + 1) for _ in range(len1 + 1)]
    # X - Scoring matrix for gaps in sequence 2
    X = [[float('-inf')] * (len2 + 1) for _ in range(len1 + 1)]
    # Y - Scoring matrix for gaps in sequence 1
    Y = [[float('-inf')] * (len2 + 1) for _ in range(len1 + 1)]
    # Backtracking matrix to reconstruct the alignment
    backtrack = [[None] * (len2 + 1) for _ in range(len1 + 1)]

    # Initialization of scores for starting gaps
    M[0][0] = 0  # Starting point of alignment (no penalty)
    for i in range(1, len1 + 1):
        X[i][0] = gap_open + (i - 1) * gap_extend  # Penalize for gaps in seq2
        M[i][0] = X[i][0]  # Set the main matrix to this score
    for j in range(1, len2 + 1):
        Y[0][j] = gap_open + (j - 1) * gap_extend  # Penalize for gaps in seq1
        M[0][j] = Y[0][j]  # Set the main matrix to this score

    # Fill in the scoring matrices
    for i in range(1, len1 + 1):
        for j in range(1, len2 + 1):
            a, b = seq1[i - 1], seq2[j - 1]
            # Assign score for match/mismatch
            score = match if a == b else mismatch

            # Calculate score for extending or starting a gap (vertical - seq2 gaps)
            X[i][j] = max(X[i - 1][j] + gap_extend,  # Extend gap
                          M[i - 1][j] + gap_open + gap_extend)  # Open new gap

            # Calculate score for extending or starting a gap (horizontal - seq1 gaps)
            Y[i][j] = max(Y[i][j - 1] + gap_extend,  # Extend gap
                          M[i][j - 1] + gap_open + gap_extend)  # Open new gap

            # Calculate score for aligning residues (match or mismatch)
            M[i][j] = max(M[i - 1][j - 1] + score,  # Match/Mismatch
                          X[i][j],  # From seq2 gap
                          Y[i][j])  # From seq1 gap

            # Determine the backtracking direction
            if M[i][j] == M[i - 1][j - 1] + score:
                backtrack[i][j] = 'D'  # Diagonal (match/mismatch)
            elif M[i][j] == X[i][j]:
                backtrack[i][j] = 'U'  # Up (seq2 gap)
            else:
                backtrack[i][j] = 'L'  # Left (seq1 gap)

    # Traceback to recover the aligned sequences
    i, j = len1, len2
    aligned_seq1 = []  # Aligned version of seq1
    aligned_seq2 = []  # Aligned version of seq2

    # Start backtracking
    while i > 0 or j > 0:
        if backtrack[i][j] == 'D':
            # Diagonal move: residues from both sequences are aligned
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append(seq2[j - 1])
            i -= 1
            j -= 1
        elif backtrack[i][j] == 'U':
            # Up move: Residue in seq1 is aligned to a gap
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append('-')
            i -= 1
        elif backtrack[i][j] == 'L':
            # Left move: Residue in seq2 is aligned to a gap
            aligned_seq1.append('-')
            aligned_seq2.append(seq2[j - 1])
            j -= 1
        else:
            # Break conditions for robustness
            break

    # Reverse the aligned sequences (since traceback constructs them backward)
    return ''.join(reversed(aligned_seq1)), ''.join(reversed(aligned_seq2)), M[len1][len2]


# Main script execution block
if __name__ == "__main__":
    import sys

    # Check if a filename is provided as a command-line argument
    if len(sys.argv) < 2:
        print("Usage:python allign.py <filename.seq>")
        sys.exit(1)

    filename = sys.argv[1]  # File containing sequences

    try:
        # Read sequences from the file
        seq1, seq2 = read_sequences_from_file(filename)
        # Perform affine gap alignment
        aligned1, aligned2, score = affine_gap_alignment(seq1, seq2)
        # Print the results
        print("Aligned Sequence 1:\n", aligned1)
        print("Aligned Sequence 2:\n", aligned2)
        print("Alignment Score:", score)
    except Exception as e:
        # Handle any errors during execution
        print("Error:", str(e))
