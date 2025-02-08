import os
from utils import running_message, run_command

@running_message
def mmseqs_makedb(input, outdir):
    db = f"{outdir}/mmseqs_db"
    if not os.path.exists(f"{outdir}/mmseqs_db"):
        cmd = f"mmseqs createdb {input} {outdir}/mmseqs_db"
        run_command(cmd, shell=True)
    else:
        print("Database already exists, using existing database")
    return db

@running_message
def mmseqs_cluster_cmd(db, outdir, threads, sensitivity, coverage=0.5, min_id = 0.3 ,e_value=1e5):
    tmp = f"{outdir}/tmp"
    os.makedirs(tmp, exist_ok=True)
    
    cluster_output = f"{outdir}/cluster_output"
    if not os.path.exists(f"{cluster_output}.index"):
        cmd = f"mmseqs cluster -s {sensitivity} -c {coverage} --min-seq-id {min_id} -e {e_value} --threads {threads} {db} {cluster_output} {tmp}"
        run_command(cmd, shell=True)
    else:
        print("Output file already exists, using existing file")
    return cluster_output

@running_message
def mmseqs_createtsv(db, cluster_output, outdir):
    tsv_path = f"{outdir}/cluster_output.tsv"
    if not os.path.exists(tsv_path):
        cmd = f"mmseqs createtsv {db} {db} {cluster_output} {tsv_path}"
        run_command(cmd, shell=True)
    else:
        print("Output file already exists, using existing file")
    return tsv_path

@running_message
def mmseqs_search(query_db, ref_db, outdir, tmp_path):
    if not os.path.exists(f"{outdir}/network.m8"):
        cmd = f"mmseqs search {query_db} {ref_db} {outdir}/network_int {tmp_path}"
        run_command(cmd, shell=True, check=True)
        cmd2 = f"mmseqs convertalis {query_db} {ref_db} {outdir}/network_int {outdir}/network.m8"
        run_command(cmd2, shell=True, check=True)
    return f"{outdir}/network.m8"
