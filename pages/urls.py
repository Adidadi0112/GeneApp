from django.urls import path
from . import views
urlpatterns = [
    path('', views.home, name='home'),
    path('input_sequence/', views.input_sequence, name='input_sequence'),
    path('input_protein/', views.input_protein, name='input_protein'),
    path('dotplot/', views.dotplot_view, name='dotplot'),
    path('align_menu/', views.align_menu, name='align_menu'),
    path('needleman_wunsch/', views.needleman_wunsch_view, name='needleman_wunsch'),
    path('msa/', views.msa_view, name='msa'),
]
