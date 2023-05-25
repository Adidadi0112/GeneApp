import io
import matplotlib.pyplot as plt

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

def upload_sequence(first_file, second_file):
  first_seq = ""
  second_seq = ""
  with open(first_file) as first_file:
    lines = first_file.readlines()[1:]
    for i in lines:
      first_seq += i
  with open(second_file) as second_file:
    lines = second_file.readlines()[1:]
    for j in lines:
      second_seq += j
  first_seq = first_seq.replace('\n',"")
  second_seq = second_seq.replace('\n',"")
  return first_seq, second_seq