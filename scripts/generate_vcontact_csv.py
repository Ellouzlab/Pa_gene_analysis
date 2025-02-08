import argparse
import pandas
from utils import *

def arguments():
    parser = argparse.ArgumentParser(description="generate vcontact csv")
    parser.add_argument('-i', '--indir', required=True, type=str, help="input directory, that is the output from genomad")
    parser.add_argument('-o', '--outfile', required=True, type=str, help="output file")
    args = parser.parse_args()
    return args

        
def grab_file_by_ending(indir, ending, exclude=None):
    for file in os.listdir(indir):
        if file.endswith(ending) and (exclude is None or not file.endswith(exclude)):
            return f"{indir}/{file}"

def main():
    args=arguments()
    
    #get path to relevant files
    summary_folder = grab_file_by_ending(args.indir, "summary", exclude='.log')
    virus_fnas = grab_file_by_ending(summary_folder, "virus.fna")
    virus_faa = grab_file_by_ending(summary_folder, "virus_proteins.faa")
    
    #read fasta files
    protein_records=read_fasta(virus_faa)
    contig_records=read_fasta(virus_fnas)
    
    #create a dictionary to store the matching proteins for each contig
    matching_dict={}
    for contig in contig_records:
        matching_dict[contig.id]=[]
        
    #fill the dictionary
    for protein in protein_records:
        contig_id='_'.join(protein.id.split('_')[:-1])
        
        if contig_id in matching_dict:
            matching_dict[contig_id].append(protein.id)
        else:
            print(f"Error! protein {protein.id} does not match any contig")
    
    # Create the DataFrame
    data = []
    for contig_id, protein_ids in matching_dict.items():
        for protein_id in protein_ids:
            data.append({"protein_id": protein_id, "contig_id": contig_id, "keywords": "None"})
    df = pd.DataFrame(data)
    
    df.to_csv(args.outfile, index=False)
    
    
if __name__ == '__main__':
    main()