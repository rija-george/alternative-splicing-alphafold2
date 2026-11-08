"""
Step 2 - Global pairwise alignment (Needleman-Wunsch) of isoform sequences
against PDB chain sequences. BLOSUM62, gap open -10, gap extend -0.5.

Code excerpt as published in Chapter 4.2 of the thesis.
Usage: python 02_pairwise_alignment.py input_file.fasta
(first record = query isoform, remaining records = targets)
"""
import sys

from Bio import Align, SeqIO
from Bio.Align import substitution_matrices
from Bio.SeqRecord import SeqRecord

matrix = substitution_matrices.load("BLOSUM62")
aligner = Align.PairwiseAligner()
aligner.mode = "global"
aligner.substitution_matrix = matrix
aligner.open_gap_score = -10
aligner.extend_gap_score = -0.5

fasta = sys.argv[1] if len(sys.argv) > 1 else "input_file.fasta"
records = list(SeqIO.parse(fasta, "fasta"))
query_seq = str(records[0].seq).strip()

for record in records[1:]:
    target_seq = str(record.seq).strip()
    alignments = aligner.align(query_seq, target_seq)
    if alignments:
        alignment = alignments[0]
        score = alignment.score
        aligned_records = [
            SeqRecord(alignment[0], id=records[0].id, description=f"Score={score}"),
            SeqRecord(alignment[1], id=record.id, description=f"Score={score}"),
        ]
        output_file = f"{records[0].id}_{record.id}_aligned.fasta"
        SeqIO.write(aligned_records, output_file, "fasta")
        print(f"{records[0].id} vs {record.id}: score {score}")
