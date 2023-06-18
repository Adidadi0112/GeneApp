import io
import matplotlib.pyplot as plt
from numpy import full
import numpy as np

def dotplot(first_seq, second_seq, window, threshold):
    fig, ax = plt.subplots()
    fig.set_size_inches(10, 10)
    for i in range(len(first_seq)):
        for j in range(len(second_seq)):
            if i >= window and j >= window:
                score = sum([1 for k in range(window) if first_seq[i - k] == second_seq[j - k]])
                if score >= threshold:
                    ax.plot(i, j, ".", color='black')
    ax.set_xlabel('Sequence 1')
    ax.set_ylabel('Sequence 2')
    ax.set_title('Dot Plot window: ' + str(window) + ', threshold: ' + str(threshold))

    # Zapisanie wykresu w formacie PNG w buforze pamięci
    buffer = io.BytesIO()
    plt.savefig(buffer, format='png')
    buffer.seek(0)

    # Zwrócenie danych wykresu jako bajtów
    return buffer.getvalue()

def upload_sequence(file):
  seq = ""
  with open(file) as file:
    lines = file.readlines()[1:]
    for i in lines:
      seq += i
  seq = seq.replace('\n',"")
  return seq


def needleman_wunsch(seq1, seq2, gap_penalty=-1, match_score=1,
                     mismatch_score=-1):
    up_arrow = "\u2191"
    left_arrow = "\u2190"
    up_left_arrow = "\u2196"

    n_rows = len(seq1) + 1  # dodatkowy wiersz
    n_columns = len(seq2) + 1  # dodatkowa kolumna

    scoring_array = full([n_rows, n_columns], 0)
    traceback_array = full([n_rows, n_columns], "-")

    # Inicjalizacja macierzy podobieństwa

    for i in range(n_rows):
        scoring_array[i, 0] = i * gap_penalty
    for j in range(len(seq2) + 1):
        scoring_array[0, j] = j * gap_penalty

    # Obliczanie wartości macierzy podobieństwa
    for i in range(1, n_rows):
        for j in range(1, n_columns):
            match = scoring_array[i - 1][j - 1] + (match_score if seq1[i - 1] ==
                                                                  seq2[j - 1] else mismatch_score)
            delete = scoring_array[i - 1][j] + gap_penalty
            insert = scoring_array[i][j - 1] + gap_penalty
            scoring_array[i][j] = max(match, delete, insert)

    # Obliczenie wartości dopasowania i zwrócenie wyniku
    aligned_seq1, aligned_seq2 = "", ""
    i, j = len(seq1), len(seq2)
    while i > 0 or j > 0:
        if i > 0 and j > 0 and scoring_array[i][j] == scoring_array[i - 1][j - 1] + (
        match_score if seq1[i - 1] == seq2[j - 1] else mismatch_score):
            aligned_seq1 = seq1[i - 1] + aligned_seq1
            aligned_seq2 = seq2[j - 1] + aligned_seq2
            traceback_array[i, j] = up_left_arrow
            i -= 1
            j -= 1
        elif i > 0 and scoring_array[i][j] == scoring_array[i - 1][j] + gap_penalty:
            aligned_seq1 = seq1[i - 1] + aligned_seq1
            aligned_seq2 = "-" + aligned_seq2
            traceback_array[i, j] = left_arrow
            i -= 1
        else:
            aligned_seq1 = "-" + aligned_seq1
            aligned_seq2 = seq2[j - 1] + aligned_seq2
            traceback_array[i, j] = up_arrow
            j -= 1
    return scoring_array[len(seq1)][len(seq2)], aligned_seq1, aligned_seq2, scoring_array, traceback_array

def central_sequence(seqs):
    M = full([len(seqs),len(seqs)],0)
    for i in range(len(seqs)):
      for j in range(len(seqs)):
        result = needleman_wunsch(seqs[i],seqs[j],-1,1,-1)
        M[i,j] = result[0]
    sum_row = np.sum(M, axis=1)
    return M

def star_msa(seqs):
    M = central_sequence(seqs)
    n = len(seqs)

    alignment = []

    central_seq_index = np.argmax(np.sum(M, axis=1))
    central_seq = seqs[central_seq_index]
    alignment.append(central_seq)

    for i in range(n):
        if i != central_seq_index:
            aligned_seq = ''
            seq = seqs[i]
            central_len = len(central_seq)
            seq_len = len(seq)
            j = 0
            k = 0
            while j < central_len or k < seq_len:
                if j < central_len and k < seq_len:
                    if central_seq[j] == seq[k]:
                        aligned_seq += central_seq[j]
                        j += 1
                        k += 1
                    else:
                        aligned_seq += '-'
                        j += 1
                        k += 1
                elif j < central_len:
                    aligned_seq += '-'
                    j += 1
                elif k < seq_len:
                    aligned_seq += seq[k]
                    k += 1
            alignment.append(aligned_seq)
    return alignment


