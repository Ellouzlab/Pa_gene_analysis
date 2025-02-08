import pandas as pd
import plotly.express as px
import os
from glob import glob
from tqdm import tqdm
from utils import *

def read_fasta_in_chunks(file_path, chunk_size=1000):
    """
    Generator function that reads a FASTA file and yields records in chunks.
    """
    with open(file_path, 'r') as file:
        chunk = []
        for line in file:
            if line.startswith(">"):
                if len(chunk) >= chunk_size:
                    yield chunk
                    chunk = []
            chunk.append(line.strip())
        if chunk:
            yield chunk

def mod_csv(csv_path, mod_csv_path):
    if not os.path.exists(mod_csv_path):
        # Load the CSV file
        df = pd.read_csv(csv_path)
        
        # Remove rows where GCA starts with "GCF"
        df = df[~df['GCA'].str.startswith("GCF")]
        
        # Group by GCA and count the number of rows for each GCA
        grouped_df = df.groupby('GCA').agg({
            'ID': 'first',         # Keep the first ID in each group (or change as needed)
            'Family': 'first',     # Keep the first Family in each group
            'Genus': 'first',      # Keep the first Genus in each group
            'Species': 'first',    # Keep the first Species in each group
            'Number': 'size'       # Count the number of rows per GCA
        }).reset_index()

        # Rename the 'Number' column to represent the count of proteins
        grouped_df.rename(columns={'Number': 'Protein_Count'}, inplace=True)
        
        # Save the modified DataFrame to the new CSV file
        grouped_df.to_csv(mod_csv_path, index=False)
    else:
        print("File already exists. Skipping processing.")
        grouped_df = pd.read_csv(mod_csv_path)
    
    return grouped_df

    
def main():
    args = arguments("data", "output")
    
    path_to_db = glob(f"{args.data}/entero*db/*.faa")
    csv_path = f"{args.output}/entero_db.csv"
    chunk_size = 1000  # Define the number of records to process per chunk

    # Only create directory if needed and CSV file does not already exist
    if not os.path.exists(csv_path):
        os.makedirs(args.output, exist_ok=True)
        
        # Open the CSV file in write mode initially to write headers
        with open(csv_path, 'w') as f:
            # Define the headers and write them
            f.write("ID,GCA,Family,Genus,Species,Number\n")

        # Process each file in the path_to_db list
        with open(csv_path, 'a') as f:
            for path in path_to_db:
                for chunk in tqdm(read_fasta_in_chunks(path, chunk_size=chunk_size), desc="Processing chunks", unit="chunk"):
                    
                    # Prepare data for the current chunk
                    chunk_data = []
                    for line in chunk:
                        if line.startswith(">"):
                            # Extract ID and parse it
                            id = line[1:]  # Removing '>'
                            split_id = id.split("___")
                            gca = split_id[0]
                            
                            try:
                                family = split_id[1]
                                genus = split_id[2]
                                species = split_id[3]
                                protein = split_id[4].split("_")[1].split(' ')[0]
                            except:
                                print(f"Error parsing ID: {id}")
                                
                            chunk_data.append(f"{id},{gca},{family},{genus},{species},{protein}")

                    # Write chunk data to CSV
                    f.write("\n".join(chunk_data) + "\n")
                    
        print(f"Data saved to {csv_path}")
    else:
        print(f"File {csv_path} already exists. Skipping processing.")
    
    mod_csv_path = f"{args.output}/entero_db_mod.csv"
    
    df = mod_csv(csv_path, mod_csv_path)
    
    # 1. Distribution of Families
    family_counts = df['Family'].value_counts().reset_index()
    family_counts.columns = ['Family', 'Genome_Count']  # Rename columns for clarity
    fig1 = px.bar(family_counts, 
                  x='Family', y='Genome_Count', 
                  labels={'Family': 'Family', 'Genome_Count': 'Genome Count'},
                  color='Family',
                  template='plotly_dark')
    fig1.update_layout(
        xaxis={'categoryorder':'total descending'},
        showlegend=False,
        font=dict(size=10),
        xaxis_title_font=dict(size=12),
        yaxis_title_font=dict(size=12),
        legend_font=dict(size=10),
        title_font=dict(size=14)
    )
    fig1.show()
    fig1.write_image(f"{args.output}/fig1_family_distribution.svg")

    # 2. Genera per Family
    genera_per_family = df.groupby(['Family', 'Genus']).size().reset_index(name='Counts')
    fig2 = px.bar(genera_per_family, x='Family', y='Counts', color='Genus', 
                  labels={'Counts': 'Number of Genera'},
                  template='plotly_dark')
    fig2.update_layout(
        xaxis={'categoryorder':'total descending'},
        font=dict(size=10),
        xaxis_title_font=dict(size=12),
        yaxis_title_font=dict(size=12),
        legend_font=dict(size=10),
        title_font=dict(size=14)
    )
    fig2.show()
    fig2.write_image(f"{args.output}/fig2_genera_per_family.svg")

    # 3. Species per Genus (Box Plot)
    species_per_genus = df.groupby('Genus')['Species'].nunique().reset_index(name='Species_Count')
    fig3 = px.box(species_per_genus, y='Species_Count', 
                  labels={'Species_Count': 'Species Count per Genus'},
                  template='plotly_dark')
    fig3.update_layout(
        font=dict(size=10),
        xaxis_title_font=dict(size=12),
        yaxis_title_font=dict(size=12),
        legend_font=dict(size=10),
        title_font=dict(size=14)
    )
    fig3.show()
    fig3.write_image(f"{args.output}/fig3_species_per_genus.svg")

    MIN_SPECIES_THRESHOLD = 10

    # Filter out species with fewer than the specified minimum threshold
    species_counts = df['Species'].value_counts()
    large_species = species_counts[species_counts >= MIN_SPECIES_THRESHOLD].index
    df_filtered = df[df['Species'].isin(large_species)]

    # Remove genera that no longer contain any species
    genera_counts = df_filtered['Genus'].value_counts()
    large_genera = genera_counts[genera_counts > 1].index
    df_filtered = df_filtered[df_filtered['Genus'].isin(large_genera)]

    # Remove families that no longer contain any genera
    family_counts = df_filtered['Family'].value_counts()
    large_families = family_counts[family_counts > 1].index
    df_filtered = df_filtered[df_filtered['Family'].isin(large_families)]

    # Create a new column that combines Genus and Species for better representation in the sunburst chart
    df_filtered['Genus_Species'] = df_filtered['Genus'] + ' ' + df_filtered['Species']

    # 4. Generate the hierarchical summary using the filtered DataFrame
    fig4 = px.sunburst(df_filtered, path=['Family', 'Genus', 'Genus_Species'], 
                       template='plotly_dark')
    fig4.update_layout(
        font=dict(size=10),
        margin=dict(t=10, l=10, r=10, b=10)
    )
    fig4.show()
    fig4.write_image(f"{args.output}/fig4_hierarchical_summary.svg")

    # 5. Generate the hierarchical summary using the filtered DataFrame (Family and Genus only)
    fig5 = px.sunburst(df_filtered, path=['Family', 'Genus'],
                       template='plotly_dark')
    fig5.update_layout(
        font=dict(size=10),
        margin=dict(t=10, l=10, r=10, b=10)
    )
    fig5.show()
    fig5.write_image(f"{args.output}/fig5_hierarchical_summary_family_genus.svg")

if __name__ == "__main__":
    main()
