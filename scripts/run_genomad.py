import argparse
import os
import subprocess as sp
from utils import *


def arguments():
    parser = argparse.ArgumentParser(description="runs genomad on fasta files")
    parser.add_argument('-i', '--indir', required=True, type=str, help="input directory")
    parser.add_argument('--gdb', default='/media/sulman/EllouzSSD1/pantoea_runs/pangenome_analysis/genomad_db', type=str, help='path to genomad db')
    parser.add_argument('-o', '--outdir', required=True, type=str, help="output directory")
    args = parser.parse_args()
    return args

def mkdir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

@running_message
def genomad(fasta_file, gdb, output):
    cmd = f"genomad end-to-end {fasta_file} {output}  {gdb}"
    sp.run(cmd, shell=True)

def main():
    args=arguments()
    for fasta_file in os.listdir(args.indir):
        if fasta_file.endswith('.fasta') or fasta_file.endswith('.fna'):
            mkdir(args.outdir)
            genomad(os.path.join(args.indir, fasta_file), args.gdb, f"{args.outdir}/{fasta_file.replace('.fasta', '')}")

if __name__ == '__main__':
    main()
