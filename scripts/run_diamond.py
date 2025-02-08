import argparse
import subprocess
import csv
import tempfile
import os
import sys
import multiprocessing

def main():
    parser = argparse.ArgumentParser(description='Run DIAMOND and output results in CSV format.')
    parser.add_argument('--query', required=True, help='Path to the query FASTA file.')
    parser.add_argument('--db', required=True, help='Path to the DIAMOND database file.')
    parser.add_argument('--out', required=True, help='Path to the output CSV file.')
    parser.add_argument('--diamond_mode', default='blastp', choices=['blastp', 'blastx'],
                        help='DIAMOND mode to use (blastp or blastx).')
    parser.add_argument('--diamond_path', default='diamond', help='Path to the DIAMOND executable.')
    parser.add_argument('--threads', default=str(multiprocessing.cpu_count()), help='Number of threads to use.')
    parser.add_argument('--write_fasta', help='Path to the output FASTA file.')
    parser.add_argument('--bitscore', type=float, default=50, help='Bitscore threshold for writing FASTA sequences.')
    args = parser.parse_args()

    # Check if bitscore is provided when write_fasta is used
    if args.write_fasta and args.bitscore is None:
        parser.error('--bitscore is required when --write_fasta is specified.')

    # Create a temporary file to store DIAMOND output
    with tempfile.NamedTemporaryFile(delete=False) as tmpfile:
        tmp_output = tmpfile.name

    try:
        # Construct the DIAMOND command with fields as separate arguments
        diamond_command = [
            args.diamond_path, args.diamond_mode,
            '--query', args.query,
            '--db', args.db,
            '--out', tmp_output,
            '--outfmt', '6',
            'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
            'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore',
            'qlen', 'slen', 'qcovhsp', 'ppos',
            '--threads', args.threads,
        ]

        print('Running DIAMOND command:', ' '.join(diamond_command))
        subprocess.run(diamond_command, check=True)

        # Initialize set to collect query IDs meeting bitscore threshold
        query_ids = set()

        # Read DIAMOND output and write to CSV with headers
        with open(tmp_output, 'r') as infile, open(args.out, 'w', newline='') as outfile:
            writer = csv.writer(outfile)
            writer.writerow([
                'query_id', 'subject_id', 'identity', 'alignment_length', 'mismatches',
                'gap_opens', 'query_start', 'query_end', 'subject_start', 'subject_end',
                'evalue', 'bitscore', 'query_length', 'subject_length', 'query_coverage', 'positive_matches'
            ])
            for line in infile:
                row = line.strip().split('\t')
                writer.writerow(row)

                # Collect query IDs if bitscore meets threshold
                if args.write_fasta:
                    bitscore = float(row[11])  # bitscore is at index 11
                    if bitscore >= args.bitscore:
                        query_ids.add(row[0])

        # Write sequences to FASTA file if requested
        if args.write_fasta and query_ids:
            with open(args.query, 'r') as fasta_in, open(args.write_fasta, 'w') as fasta_out:
                write_seq = False
                seq_id = ''
                for line in fasta_in:
                    if line.startswith('>'):
                        seq_id = line[1:].split()[0]  # Extract ID without '>'
                        write_seq = seq_id in query_ids
                    if write_seq:
                        fasta_out.write(line)

    except subprocess.CalledProcessError as e:
        print('DIAMOND failed with error code:', e.returncode)
        sys.exit(1)
    finally:
        # Clean up temporary file
        os.remove(tmp_output)

if __name__ == '__main__':
    main()
