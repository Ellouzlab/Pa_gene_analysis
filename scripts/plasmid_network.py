from utils import *

@running_message
def makedb(input, database):
    if not os.path.exists(f"{database}.dmnd"):
        cmd=f"diamond makedb --in {input} -d {database}"
        run_command(cmd, shell=True)
    else:
        print("Database already exists, using existing database")
        
@running_message
def diamond(input, database, tsv_path, bitscore, threads, sensitivity=1):
    sensitivity_types = {
        0: "--faster",
        1: "",
        2: "--sensitive",
        3: "--more-sensitive",
        4: "--very-sensitive",
        5: "--ultra-sensitive"
    }
    if not sensitivity in sensitivity_types:
        raise ValueError("Sensitivity should be between 0 and 4")
    
    if not os.path.exists(tsv_path):
        cmd = f"diamond blastp {sensitivity_types[sensitivity]} -q {input} -d {database} --log --min-score {bitscore} -o {tsv_path} --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovhsp ppos --threads {threads}"
        run_command(cmd, shell=True)
    else: 
        print("TSV file already exists, using existing file")
    
def plasmid_network(plsdb, entero_db, outdir):
    os.makedirs(f"{outdir}/plsdb", exist_ok=True)
    plsdb_database = f"{outdir}/plsdb/plsdb"
    makedb(plsdb, plsdb_database)
    
    os.makedirs(f"{outdir}/enterobacterales", exist_ok=True)
    entero_database = f"{outdir}/enterobacterales/enterobacterales"
    makedb(entero_db, entero_database)
    
    tsv_path = f"{outdir}/plsdb_vs_enterobacterales.tsv"
    diamond(plsdb, entero_database, tsv_path, 50, 32, sensitivity=1)
    
    

def main():
    args = arguments("data", "outdir")
    entero_db = f"{args.data}/enterobacterales_db/enterobacterales_db.faa"
    plsdb = f"{args.data}/plsdb_reps.fasta"
    
    os.makedirs(args.outdir, exist_ok=True)
    init_logging(f"{args.outdir}/plasmid_network.log")
    plasmid_network(plsdb, entero_db, args.outdir)

if __name__ == "__main__":
    main()