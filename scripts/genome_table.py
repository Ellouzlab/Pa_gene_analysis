from utils import *
import pandas as pd
import os, json
from collections import defaultdict

if __name__ == "__main__":
    args = arguments("data", "output")
    
    assemblies = [strain.split("/")[-1].split(".")[0].replace('Pantoea_agglomerans_', '') for strain in os.listdir(f"{args.data}/assemblies")]
    
    df = pd.read_csv(f"{args.data}/ncbi_dataset.tsv", sep="\t")
    df["Organism Infraspecific Names Strain"] = df["Organism Infraspecific Names Strain"].map(lambda x: str(x).replace(" ", "_"))
    df = df.drop_duplicates(subset=['Organism Infraspecific Names Strain'])
    df = df[df['Organism Infraspecific Names Strain'].isin(assemblies)]
    df["Organism Infraspecific Names Strain"] = df["Organism Infraspecific Names Strain"].map(lambda x: x.replace("_", " ") if not "MM2021" in x else x)
    
    
    JSONL_path = f"{args.data}/assembly_data_report.jsonl"

    # Open the file and read each line as a JSON object
    with open(JSONL_path, 'r') as f:
        json_lines = [json.loads(line) for line in f]

    json_dict = {}
    for json_line in json_lines:
        try:
            strain_name=json_line['organism']['infraspecificNames']['strain']
        except:
            try:
                strain_name=json_line['organism']['infraspecificNames']['isolate']
            except:
                print(json_line)
                
        if not strain_name in json_dict:
            json_dict[strain_name]=json_line
        else:
            print(f"Duplicate strain name found: {strain_name}")
    
    print(json_dict.keys())
    def get_host_info(strain_name):
        try:
            return json_dict[strain_name]['assemblyInfo']['biosample'].get('host', '')
        except KeyError:
            return ''
    
    def get_location_info(strain_name):
        try:
            return json_dict[strain_name]['assemblyInfo']['biosample'].get('geoLocName', '')
        except KeyError:
            return ''

    # Apply the mapping to the 'Host' column, using an empty string for missing data
    df["Host"] = df["Organism Infraspecific Names Strain"].map(get_host_info)
    df["Location"] = df["Organism Infraspecific Names Strain"].map(get_location_info)
    df['Country']=df.Location.map(lambda x : x.split(":")[0])
    df['FullName']="Pantoea agglomerans "+df["Organism Infraspecific Names Strain"]
    
    df.to_csv(f"genome_table.tsv", sep="\t", index=False)