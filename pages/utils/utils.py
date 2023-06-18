from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Seq import Seq
from bs4 import BeautifulSoup
from Bio import Entrez, SeqIO
def validate_sequence(sequence):
    # Logika walidacji sekwencji
    valid_characters = 'ACDEFGHIKLMNPQRSTVWY'
    for char in sequence:
        if char not in valid_characters:
            return False
    return True

def find_sequence_id(sequence):
    protein_sequence = Seq(sequence)
    result_handle = NCBIWWW.qblast('blastp', 'nr', protein_sequence)

    # Parsowanie wyników BLAST
    blast_records = NCBIXML.read(result_handle)

    # Przygotowanie tabeli HTML z obramowaniem
    table_html = '<table style="border-collapse: collapse; border: 1px solid black;">'
    table_html += '<tr><th style="border: 1px solid black;">Alignment Title</th>'
    table_html += '<th style="border: 1px solid black;">HSP Details</th></tr>'

    # Dodawanie danych do tabeli HTML
    alignments = []
    for alignment in blast_records.alignments:
        for hsp in alignment.hsps:
            # Tworzenie bardziej informacyjnego tytułu HSP
            hsp_title = f'Score: {hsp.score}, E-value: {hsp.expect}'

            alignment_data = {
                'title': alignment.title[:80],
                'hsp': str(hsp_title)
            }
            alignments.append(alignment_data)

    # Generowanie wierszy tabeli HTML
    for alignment_data in alignments:
        table_html += '<tr>'
        table_html += f'<td style="border: 1px solid black;">{alignment_data["title"]}</td>'
        table_html += f'<td style="border: 1px solid black;">{alignment_data["hsp"]}</td>'
        table_html += '</tr>'

    table_html += '</table>'

    return table_html


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
def get_sequence_from_id(sequence_id):
    Entrez.email = 'adanos1312@wp.pl'

    try:
        handle = Entrez.efetch(db='protein', id=sequence_id, rettype='fasta', retmode='text')
        record = SeqIO.read(handle, 'fasta')
        handle.close()
        return str(record.seq)
    except Exception as e:
        print(f'Wystąpił błąd: {e}')
        return None
