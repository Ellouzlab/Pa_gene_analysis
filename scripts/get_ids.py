import argparse
from Bio import SeqIO

def arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="Input file in fasta format", required=True)
    parser.add_argument("-o", "--output", help="Output file with ids", required=True)
    args = parser.parse_args()
    return args

def read_fasta(input_file):
    with open(input_file, "r") as handle:
        records = list(SeqIO.parse(handle, "fasta"))
    return records

def main():
    args = arguments()
    records = read_fasta(args.input)
    txt='\n'.join([record.id.split('|')[0] for record in records])
    
    with open(args.output, "w") as handle:
        handle.write(txt)

if __name__ == "__main__":
    main()