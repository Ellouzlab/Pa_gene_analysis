from utils import *
from mmseqs_utils import *
import pandas as pd
import glob

def searching():
    pass

def run_vfdb_search(query, vfdb, vfdb_desc, outdir):
    if not os.path.exists(f"{outdir}/vfdb_search.tsv"):
        os.makedirs(outdir, exist_ok=True)
        init_logging(f"{outdir}/vfdb_search.log")

        tmp = f"{outdir}/tmp"
        os.makedirs(tmp, exist_ok=True)
        
        vfdb_tmp_path = f"{tmp}/vfdb"
        os.makedirs(vfdb_tmp_path, exist_ok=True)
        vfdb_db = mmseqs_makedb(vfdb, vfdb_tmp_path)
        
        query_tmp_path = f"{tmp}/query"
        os.makedirs(query_tmp_path, exist_ok=True)
        query_db = mmseqs_makedb(query, query_tmp_path)
        
        search_tmp_path = f"~/tmp/search"
        os.makedirs(search_tmp_path, exist_ok=True)
        
        search_output_dir = f"{tmp}/search_output"
        os.makedirs(search_output_dir, exist_ok=True)
        search_output = mmseqs_search(query_db, vfdb_db, search_output_dir, search_tmp_path)
        
        
        columns = ["query", "target", "pident", "alnlen", "mismatch", "numgapopen", "qstart", "qend", "tstart", "tend", "evalue", "bitscore"]
        df = pd.read_csv(search_output, sep="\t", header=None, names=columns)
        df = df[(df['bitscore'] > 50) & (df['pident'] > .70)]
        df = df.loc[df.groupby('query')['bitscore'].idxmax()]
        df['target'] = df['target'].map(lambda x: x.split("(")[0])
        print(df)
        
        fasta_records = read_fasta(vfdb)
        fasta_dict = {record.id.split("(")[0]: f"VFC{record.description.split('(VFC')[-1].split(')')[0]}" for record in fasta_records}
        df.target = df.target.map(fasta_dict)
        
        descript_df = pd.read_csv(vfdb_desc)
        descript_df = descript_df.drop_duplicates(subset=["VFCID"])
        print(descript_df.columns)
        
        descript_df.index = descript_df.VFCID
        df["VF_cat"] = df.target.map(descript_df["VFcategory"])
        df.to_csv(f"{outdir}/vfdb_search.tsv", sep="\t", index=False)
    return f"{outdir}/vfdb_search.tsv"
    
    
    
    

def main(args):
    queries = [entry for entry in glob.glob(f"{args.bakta}/*/*.faa") if not "hypothetical" in entry]
    result_paths = []
    for query in queries:
        outdir = f"{args.outdir}/{query.split('/')[-1].split('.')[0]}"
        
        result_path = run_vfdb_search(query = query, vfdb = args.vfdb, vfdb_desc = args.vfdb_descs, outdir = outdir)
        result_paths.append(result_path)
    
    df_list = []
    for result_path in result_paths:
        df = pd.read_csv(result_path, sep="\t")
        df_list.append(df)
    df = pd.concat(df_list)
    df.to_csv(f"{args.outdir}/vfdb_search.tsv", sep="\t", index=False)
    
    
if __name__ == "__main__":
    args = arguments("bakta", "vfdb", "vfdb_descs", "outdir")
    main(args)