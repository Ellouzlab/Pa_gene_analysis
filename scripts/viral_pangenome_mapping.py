from utils import *
from mmseqs_utils import *
import os
import pandas as pd

@main_logging
def main():
    args=arguments("genes", "viral_db", "outdir")
    

    os.makedirs(args.outdir, exist_ok=True)
    init_logging(f"{args.outdir}/viral_pangenome_mapping.log")
    
    viral_db_path = f"{args.outdir}/viral_db"
    os.makedirs(viral_db_path, exist_ok=True)
    viral_db = mmseqs_makedb(args.viral_db, viral_db_path)

    genes_db_path = f"{args.outdir}/genes"
    os.makedirs(genes_db_path, exist_ok=True)
    genes_db = mmseqs_makedb(args.genes, genes_db_path)

    m8_file = mmseqs_search(genes_db, viral_db, args.outdir, f"{args.outdir}/tmp")
    columns = ["query", "target", "pident", "alnlen", "mismatch", "numgapopen", "qstart", "qend", "tstart", "tend", "evalue", "bitscore"]

    df = pd.read_csv(m8_file, sep="\t", header=None, names=columns)
    df['query'] = df.query.map(lambda x: x.split('_')[0])
    df = df[df['bitscore'] > 50].groupby("query").apply(lambda df: df.loc[df['bitscore'].idxmax()]).to_csv(f"{args.outdir}/viral_pangenome_mapping.tsv", sep="\t", index=False)
    print(df)

if __name__ == "__main__":
    main()
