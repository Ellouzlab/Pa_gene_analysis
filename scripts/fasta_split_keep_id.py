import argparse, os
from utils import *

def arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--fasta', required=True, type=str, help="Path to fasta file")
    parser.add_argument('--output', required=True, type=str, help="output directory")
    return parser.parse_args()

def main():
    args = arguments()

    records = read_fasta(args.fasta)

    os.makedirs(args.output, exist_ok=True)
    for record in records:
        fasta_name = f"{args.output}/{record.id}.fasta".replace("||","__")
        write_fasta(fasta_name, [record])
    print("Done")

if __name__ == "__main__":
    main()