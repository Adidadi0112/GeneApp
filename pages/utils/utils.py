from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Seq import Seq
from bs4 import BeautifulSoup
def validate_sequence(sequence):
    # Logika walidacji sekwencji
    valid_characters = 'ACDEFGHIKLMNPQRSTVWY'
    for char in sequence:
        if char not in valid_characters:
            return False
    return True

def find_sequence(sequence):
    protein_sequence = Seq(sequence)
    result_handle = NCBIWWW.qblast('blastp', 'nr', protein_sequence)

    # Parsowanie wyników BLAST
    blast_records = NCBIXML.read(result_handle)

    # Przygotowanie tabeli HTML
    table_html = "<table>"
    table_html += "<tr><th>Alignment Title</th><th>HSP</th></tr>"

    # Dodawanie danych do tabeli HTML
    alignments = []
    for alignment in blast_records.alignments:
        for hsp in alignment.hsps:
            alignment_data = {
                'title': alignment.title[:80],
                'hsp': str(hsp)
            }
            alignments.append(alignment_data)

    for alignment_data in alignments[:5]:  # Wyświetlaj tylko 5 najlepszych dopasowań
        table_html += "<tr>"
        table_html += f"<td>{alignment_data['title']}</td>"
        table_html += f"<td>{alignment_data['hsp']}</td>"
        table_html += "</tr>"

    table_html += "</table>"

    # Utworzenie obiektu BeautifulSoup
    soup = BeautifulSoup(table_html, 'html.parser')

    # Zwrócenie tabeli HTML jako kodu HTML
    return str(soup)
