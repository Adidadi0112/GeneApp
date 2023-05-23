from django.urls import path
from . import views
urlpatterns = [
    path('', views.home, name='home'),
    path('input_sequence/', views.input_sequence, name='input_sequence'),
    path('input_protein/', views.input_protein, name='input_protein'),
    path('align_sequences/', views.align_sequences, name='align_sequences'),
]
