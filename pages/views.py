import base64
from django.shortcuts import render, redirect
from .forms import DataForm
from .utils.utils import validate_sequence, find_sequence
from .utils.align_utils import dotplot, upload_sequence, needleman_wunsch


def home(request):
    return render(request, "home.html", {})


def input_sequence(request):
    if request.method == 'POST':
        form = DataForm(request.POST)
        if form.is_valid():
            input_data = form.cleaned_data['input_data']
            if validate_sequence(input_data):
                result = find_sequence(input_data)
                context = {
                    'form': form,
                    'input_data': input_data,
                    'result': result
                }
                return render(request, 'input_sequence.html', context)
            else:
                form.add_error('input_data', 'input data is not correct')
    else:
        form = DataForm()

    return render(request, 'input_sequence.html', {'form': form})


def input_protein(request):
    return render(request, "input_protein.html", {})


def align_menu(request):
    return render(request, "align_menu.html", {})


def dotplot_view(request):
    if request.method == 'POST':
        first_seq = request.POST.get('first_seq', '')
        second_seq = request.POST.get('second_seq', '')
        window = int(request.POST.get('window', '0'))
        threshold = int(request.POST.get('threshold', '0'))

        first_file = request.FILES.get('first_seq_file')
        second_file = request.FILES.get('second_seq_file')
        if first_file and second_file:
            first_seq, second_seq = upload_sequence(first_file, second_file)

        plot_data = dotplot(first_seq, second_seq, window, threshold)

        # Konwersja danych wykresu na base64
        plot_base64 = base64.b64encode(plot_data).decode('utf-8')

        context = {
            'plot_base64': plot_base64,
            'first_seq': first_seq,
            'second_seq': second_seq,
            'window': window,
            'threshold': threshold
        }

        return render(request, "dotplot.html", context)
    else:
        return render(request, "dotplot.html", {})


def needleman_wunsch_view(request):
    if request.method == 'POST':
        seq1 = request.POST.get('seq1', '')
        seq2 = request.POST.get('seq2', '')

        if seq1 and seq2:
            score, aligned_seq1, aligned_seq2, scoring_array, traceback_array = needleman_wunsch(seq1, seq2)

            context = {
                'seq1': seq1,
                'seq2': seq2,
                'score': score,
                'aligned_seq1': aligned_seq1,
                'aligned_seq2': aligned_seq2,
                'scoring_array': scoring_array,
                'traceback_array': traceback_array,
                'show_result': True  # Dodaj nowy klucz kontekstu 'show_result' i ustaw go na True
            }

            return render(request, "needleman_wunsch.html", context)

    return render(request, "needleman_wunsch.html", {})


def msa_view(request):
    if request.method == 'POST':
        sequences = request.session.get('sequences', [])
        action = request.POST.get('action')

        if action == 'add':
            sequence = request.POST.get('sequence')
            sequences.append(sequence)
            request.session['sequences'] = sequences
        elif action == 'delete':

            if sequences:
                sequences.pop()

                request.session['sequences'] = sequences

        return redirect('msa')

    else:
        sequences = request.session.get('sequences', [])

    context = {
        'sequences': sequences
    }

    return render(request, 'msa.html', context)


