import argparse
import os
from utils import *

def run_phagcn(file_path):
    pass

def arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--genomad_folders', type=str, required=True, help="Path to folder containing all the results from genomad")
    parser.add_argument('--output', type=str, required=True, help="output file name")
    args = parser.parse_args()
    return args

def grab_folder(folder_path, query, exclude):
    for folder in os.listdir(folder_path):
        if query in folder and exclude not in folder:
            return folder

def prep_records(viral_records, bacteria_name):
    counter=1
    for record in viral_records:
        length = len(record.seq)
        record.id = f"{bacteria_name}||phage_{counter}||{length}bp"
        counter+=1
    return viral_records

def main():
    args = arguments()
    all_viruses = []
    for bacteria_folder in os.listdir(args.genomad_folders):
        bacteria_name= bacteria_folder.split(".")[0].split("/")[-1]
        bacteria_folder_path =f"{args.genomad_folders}/{bacteria_folder}"
        
        #Grab summary folder
        summary_folder = grab_folder(bacteria_folder_path, "summary", ".log")
        summary_folder_path = f"{bacteria_folder_path}/{summary_folder}"

        virus_fna= grab_folder(summary_folder_path, "virus.fna", ".log")
        virus_fna_path = f"{summary_folder_path}/{virus_fna}"

        viral_records=read_fasta(virus_fna_path)
        viral_records_mod = prep_records(viral_records, bacteria_name)
        all_viruses.extend(viral_records_mod)
    
    write_fasta(args.output, all_viruses)
        

if __name__ == "__main__":
    main()
    



