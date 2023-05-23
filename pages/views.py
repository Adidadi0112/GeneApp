
from django.shortcuts import render
from .forms import DataForm
from .utils import validate_sequence, find_sequence

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

def align_sequences(request):
    return render(request, "align_sequences.html", {})
