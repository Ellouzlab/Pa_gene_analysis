import argparse
import os
import subprocess as sp
import threading


def arguments():
    parser = argparse.ArgumentParser(description="runs bakta on fasta files")
    parser.add_argument('-i', '--indir', required=True, type=str, help="input directory")
    parser.add_argument('-o', '--outdir', required=True, type=str, help="output directory")
    parser.add_argument('-t', '--threads', default=32, type=int, help='number of threads to use')
    args = parser.parse_args()
    return args

def mkdir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

def bakta(fasta_file, threads,output):
    cmd = f"prokka --cpus {threads} --force --compliant --outdir {output} {fasta_file}"
    sp.run(cmd, shell=True)

def main():
    args=arguments()
    for fasta_file in os.listdir(args.indir):
        if fasta_file.endswith('.fasta'):
            mkdir(args.outdir)
            bakta(os.path.join(args.indir, fasta_file), args.threads,f"{args.outdir}/{fasta_file.replace('.fasta', '')}")

if __name__ == '__main__':
    main()
