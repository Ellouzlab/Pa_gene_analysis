import argparse
import os
import matplotlib.pyplot as plt
import pandas as pd
from utils import *


def arguments():
    parser = argparse.ArgumentParser(description="Analyze results from checkv")
    parser.add_argument('-i', '--indir', required=True, type=str, help="input directory, that is the output from checkv")
    parser.add_argument('-o', '--outdir', required=True, type=str, help="output directory")
    parser.add_argument('-t', '--threshold', required=False, type=float, default=90.0, help='threshold for completeness')

    args = parser.parse_args()
    return args

def main():
    args=arguments()

    completeness_file = f"{args.indir}/quality_summary.tsv"
    completeness_data = pd.read_csv(completeness_file, sep='\t')
    print(completeness_data.head())

    virus_file=f"{args.indir}/viruses.fna"
    virus_records=read_fasta(virus_file)

    provirus_file=f"{args.indir}/proviruses.fna"
    provirus_records=read_fasta(provirus_file)

    all_records = virus_records + provirus_records
    
    filtered_records = []
    for record in all_records:
        try:
            if record.id.endswith("_1"):
                record.id = record.id[:-2]
            completeness=completeness_data.loc[completeness_data['contig_id']==record.id]['completeness'].item()
            viral_gene_count=completeness_data.loc[completeness_data['contig_id']==record.id]['viral_genes'].item()
            print(record.id, completeness)
            if completeness>args.threshold and viral_gene_count>0:
                record.id = f"{record.id}||complete={completeness}"
                record.description = ""
                filtered_records.append(record)
        except:
            print(f"Error! record {record.id} does not have a completeness value")

    outdir = makedir(args.outdir)
    write_fasta(f"{outdir}/complete_decontaminated_viruses.fna", filtered_records)

if __name__ == '__main__':
    main()