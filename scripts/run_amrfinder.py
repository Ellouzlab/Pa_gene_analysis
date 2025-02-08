from utils import *
import glob, os, subprocess
import pandas as pd

@running_message
def amrfinder(prot_file, out_file):
    if not os.path.exists(out_file):
        cmd = f"amrfinder -p {prot_file} -o {out_file} --plus --threads 32"
        run_command(cmd)

def main(args):
    protein_files = glob.glob(f"{args.bakta}/*/*.faa")
    
    outfiles=[]
    for prot_file in protein_files:
        if not "hypothetical" in prot_file:
            basename = prot_file.split("/")[-1].split(".")[0]
            out_file = f"{args.outdir}/{basename}.tsv"
            amrfinder(prot_file, out_file)
            outfiles.append(out_file)
    
    df_list = []
    for file in outfiles:
        df = pd.read_csv(file, sep="\t")
        name = file.split("/")[-1].split(".")[0]
        df["name"] = name
        df_list.append(df)
    df = pd.concat(df_list)
    df.to_csv(f"{args.outdir}/amrfinder.tsv", sep="\t", index=False)
        
    

if __name__ == '__main__':
    args = arguments("bakta", "outdir")
    os.makedirs(args.outdir, exist_ok=True)
    main(args)