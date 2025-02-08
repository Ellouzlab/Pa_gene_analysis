from utils import *
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def translate_fasta(input_file, output_file):
    # Read the FASTA file
    records = SeqIO.parse(input_file, 'fasta')
    
    # Translate each record and update description
    translated_records = []
    for record in records:
        # Translate the sequence, assuming it's a valid DNA sequence for translation
        translated_seq = record.seq.translate(to_stop=False)
        
        # Create a new SeqRecord with the same ID, empty description and translated sequence
        new_record = SeqRecord(translated_seq, id=record.id, description="")
        translated_records.append(new_record)
    
    # Write the translated sequences to an output file
    SeqIO.write(translated_records, output_file, 'fasta')

if __name__ == '__main__':
    args = arguments("input", "output")
    translate_fasta(args.input, args.output)

