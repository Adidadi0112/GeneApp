from django import forms

class DataForm(forms.Form):
    input_data = forms.CharField(label='Wprowadź dane', max_length=100)
