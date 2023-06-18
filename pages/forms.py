from django import forms


class DataForm(forms.Form):
    input_data = forms.CharField(label='Input:', widget=forms.Textarea(attrs={'rows': 1, 'cols': 100}))

