import base64
from django.shortcuts import render
from .forms import DataForm
from .utils.utils import validate_sequence, find_sequence
from .utils.align_utils import dotplot, upload_sequence

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
                print(input_data)
            else:
                form.add_error('input_data', 'input data is not correct')
    else:
        form = DataForm()

    return render(request, 'input_sequence.html', {'form': form})

def input_protein(request):
    return render(request, "input_protein.html", {})


def align_sequences(request):
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

        return render(request, "align_sequences.html", context)
    else:
        return render(request, "align_sequences.html", {})

