import argparse
from Bio import SeqIO
import os

def main():
    parser = argparse.ArgumentParser(description="Split a FASTA file into separate files for each record.")
    parser.add_argument("input_file", type=str, help="Input FASTA file containing DNA records.")
    parser.add_argument("output_dir", type=str, help="Directory where the output FASTA files will be stored.")
    
    args = parser.parse_args()
    
    # Ensure the output directory exists
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    
    # Read and process each record
    record_iterator = SeqIO.parse(args.input_file, "fasta")
    for index, record in enumerate(record_iterator, start=1):
        output_path = os.path.join(args.output_dir, f"virus_{index}.fasta")
        SeqIO.write(record, output_path, "fasta")
        print(f"Written {output_path}")

if __name__ == "__main__":
    main()

