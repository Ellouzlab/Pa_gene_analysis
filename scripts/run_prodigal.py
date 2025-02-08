import argparse
from Bio import SeqIO
import os
import concurrent.futures
import subprocess
from tqdm import tqdm

def parse_args():
    parser = argparse.ArgumentParser(description="Run Prodigal on multiple threads by splitting the input FASTA file into chunks.")
    parser.add_argument('-i', '--input', required=True, help="Input FASTA file")
    parser.add_argument('-o', '--output_dir', required=True, help="Output directory for Prodigal results")
    parser.add_argument('-c', '--chunk_size', type=int, default=1000, help="Number of sequences per chunk (default: 100)")
    parser.add_argument('-t', '--threads', type=int, default=32, help="Number of threads to use (default: 4)")
    parser.add_argument('-f', '--final_output', required=True, help="Final output file combining all Prodigal results")
    return parser.parse_args()

def split_fasta(input_file, output_dir, chunk_size):
    os.makedirs(output_dir, exist_ok=True)
    chunk = []
    chunk_count = 0

    with open(input_file, "r") as infile:
        for record in tqdm(SeqIO.parse(infile, "fasta"), desc="Splitting FASTA", unit="chunks"):
            chunk.append(record)
            if len(chunk) == chunk_size:
                chunk_file = os.path.join(output_dir, f"chunk_{chunk_count}.fasta")
                with open(chunk_file, "w") as outfile:
                    SeqIO.write(chunk, outfile, "fasta")
                chunk = []
                chunk_count += 1

        # Write the remaining chunk if it's not empty
        if chunk:
            chunk_file = os.path.join(output_dir, f"chunk_{chunk_count}.fasta")
            with open(chunk_file, "w") as outfile:
                SeqIO.write(chunk, outfile, "fasta")

def run_prodigal(chunk_file, output_dir):
    output_file = os.path.join(output_dir, f"{os.path.basename(chunk_file)}.genes")
    cmd = f"prodigal -i {chunk_file} -a {output_file} -p meta"
    try:
        with open(os.devnull, 'w') as devnull:
            subprocess.run(cmd, shell=True, check=True, stdout=devnull, stderr=devnull)
    except subprocess.CalledProcessError as e:
        print(f"Error running Prodigal on {chunk_file}: {e}")
    return output_file

def run_prodigal_parallel(input_dir, output_dir, threads):
    os.makedirs(output_dir, exist_ok=True)
    chunk_files = [os.path.join(input_dir, f) for f in os.listdir(input_dir) if f.endswith(".fasta")]
    
    with tqdm(total=len(chunk_files), desc="Running Prodigal", unit="chunks") as pbar:
        with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
            futures = [executor.submit(run_prodigal, chunk_file, output_dir) for chunk_file in chunk_files]
            for future in concurrent.futures.as_completed(futures):
                try:
                    future.result()
                    pbar.update(1)
                except Exception as e:
                    print(f"Error: {e}")

def combine_results(input_dir, output_file):
    with open(output_file, 'w') as outfile:
        for filename in tqdm(sorted(os.listdir(input_dir)), desc="Combining results", unit="files"):
            if filename.endswith(".genes"):
                chunk_file_path = os.path.join(input_dir, filename)
                with open(chunk_file_path, 'r') as infile:
                    for line in infile:
                        outfile.write(line)

if __name__ == "__main__":
    args = parse_args()
    
    split_fasta(args.input, args.output_dir, args.chunk_size)
    run_prodigal_parallel(args.output_dir, args.output_dir, args.threads)
    combine_results(args.output_dir, args.final_output)
