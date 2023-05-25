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
    count = 0  # Licznik dopasowań
    for alignment in blast_records.alignments:
        for hsp in alignment.hsps:
            table_html += "<tr>"
            table_html += f"<td>{alignment.title[:80]}</td>"
            table_html += f"<td>{str(hsp)}</td>"
            table_html += "</tr>"
            count += 1

            if count == 5:  # Przerwij pętlę po znalezieniu 5 rekordów
                break

        if count == 5:  # Przerwij pętlę po znalezieniu 5 rekordów
            break

    table_html += "</table>"

    # Utworzenie obiektu BeautifulSoup
    soup = BeautifulSoup(table_html, 'html.parser')

    # Zwrócenie tabeli HTML jako kodu HTML
    return str(soup)
