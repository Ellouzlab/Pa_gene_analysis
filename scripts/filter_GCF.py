import argparse
from Bio import SeqIO

def filter_fasta(input_file, output_file):
    # Read the input FASTA file and filter out records that start with 'GCF'
    with open(input_file, "r") as infile, open(output_file, "w") as outfile:
        filtered_records = (record for record in SeqIO.parse(infile, "fasta") if not record.id.startswith("GCF"))
        SeqIO.write(filtered_records, outfile, "fasta")
    print(f"Filtered FASTA file saved as {output_file}")

def main():
    parser = argparse.ArgumentParser(description="Filter out FASTA records with headers starting with 'GCF'")
    parser.add_argument("input", help="Input FASTA file")
    parser.add_argument("output", help="Output FASTA file")
    args = parser.parse_args()

    filter_fasta(args.input, args.output)

if __name__ == "__main__":
    main()
