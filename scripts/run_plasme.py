from utils import *
import glob
from tqdm import tqdm
import os

def main():
    args = arguments("indir", "outdir", "plasme_path", "db")
    os.makedirs(args.outdir, exist_ok=True)
    init_logging(f"{args.outdir}/run_plasme.log")
    input_files = glob.glob(f"{args.indir}/*.fna")

    for input_file in tqdm(input_files):
        try:
            os.makedirs("~/temp", exist_ok=True)
            output_path = f"{args.outdir}/{os.path.basename(input_file).replace('.fna', '.tsv')}"
            os.makedirs(output_path, exist_ok=True)
            cmd=f"python {args.plasme_path} {input_file} {output_path}/plasmids -t 47 --temp ~temp -d {args.db}"
            run_command(cmd)
        except Exception as e:
            logging.error(f"An error occurred: {e}", exc_info=True)
            print(f"An error occurred: {e}")

if __name__ == "__main__":
    main()