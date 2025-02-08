import argparse
import pandas
from utils import *


def arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--phage_csv', required=True, type = str, help="Path to csv file containing phage data")
    parser.add_argument('--fasta', required=True, type = str, help="Path to fasta file containing phage sequences")
    parser.add_argument('--output', required=True, type = str, help="output fasta name")
    return parser.parse_args()


def main():
    args=arguments()
    initial_records=read_fasta(args.fasta)
    
    df=pd.read_csv(args.phage_csv)

    for record in initial_records:
        try:
            record.id = f"{record.id}||{df.loc[df['contig_name']==record.id][' prediction'].values[0]}"
        except:
            record.id = f"{record.id}||unknown"

    write_fasta(args.output, initial_records)
    


if __name__ == "__main__":
    main()