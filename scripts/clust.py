import subprocess
import pandas as pd
from sklearn.manifold import TSNE
import plotly.express as px
import argparse
import os

def run_blast(fasta_file, num_cores):
    db_name = "protein_db"
    subprocess.run(['makeblastdb', '-in', fasta_file, '-dbtype', 'prot', '-out', db_name], text=True)
    output_file = "blast_output.txt"
    subprocess.run(['blastp', '-query', fasta_file, '-db', db_name, '-out', output_file, '-outfmt', '6 qseqid sseqid pident', '-num_threads', str(num_cores)], text=True)
    return output_file

def construct_distance_matrix(blast_output):
    data = pd.read_csv(blast_output, sep='\t', names=['query', 'subject', 'identity'])
    data['distance'] = 100 - data['identity']
    data = data.groupby(['query', 'subject']).agg({'distance': 'min'}).reset_index()
    matrix = data.pivot(index='query', columns='subject', values='distance').fillna(0)
    matrix.to_csv("distance_matrix.csv")
    return matrix


def main():
    parser = argparse.ArgumentParser(description="Run BLAST to create a distance matrix and perform 3D t-SNE.")
    parser.add_argument('fasta_file', type=str, help="The path to the FASTA file.")
    parser.add_argument('num_cores', type=int, help="Number of processor cores to use.")
    args = parser.parse_args()

    distance_matrix_file = "distance_matrix.csv"
    
    # Check if the distance matrix file already exists
    if os.path.exists(distance_matrix_file):
        print("Using existing distance matrix.")
        distance_matrix = pd.read_csv(distance_matrix_file, index_col=0)
    else:
        print("Running BLAST and constructing distance matrix.")
        blast_output = run_blast(args.fasta_file, args.num_cores)
        distance_matrix = construct_distance_matrix(blast_output)


if __name__ == "__main__":
    main()
