from utils import *
import pandas as pd
from tqdm import tqdm
import os
import logging
from umap import UMAP
import plotly.express as px
from mmseqs_utils import *
import glob
from sklearn.preprocessing import StandardScaler
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.cluster.hierarchy import dendrogram
import matplotlib.pyplot as plt
from sklearn.cluster import DBSCAN

@running_message
def add_rgp(gff, rgp_tsv, gene_fam, output):
    output_file = f"{output}/rgp_gff.tsv"
    if not os.path.exists(output_file) or os.path.getsize(output_file) == 0:
        rgp_df = pd.read_csv(rgp_tsv, sep="\t")
        
        columns = ["contig", "type", "source", "start", "end", "score", "strand", "frame", "attributes"]
        gff_df = pd.read_csv(gff, header=None, names=columns, index_col=False)
        print(gff_df)
        gff_df = gff_df[gff_df.source == "CDS"]
        gff_df["protein_id"] = gff_df.attributes.map(lambda x: x.split("ID=")[1].split(";")[0] if "ID=" in x else x.split("gene_id ")[1].split(";")[0])

        gene_fam = pd.read_csv(gene_fam, sep="\t", names=["rep", "protein_id", "something"])
        gene_fam_dict = gene_fam.set_index("protein_id")["rep"].to_dict()

        # Map 'rep' column in gff_df
        gff_df["rep"] = gff_df["protein_id"].map(gene_fam_dict)
        gff_df["rgp_id"] = None

        # Find start and end indices
        start_index = gff_df[gff_df["protein_id"].isin(rgp_df["start"])].index
        end_index = gff_df[gff_df["protein_id"].isin(rgp_df["stop"])].index

        coords = zip(start_index, end_index)
        coords = list(coords)  # Convert to list to reuse multiple times

        # Initialize variables before loop
        all_protein_indices = []
        all_rgp_ids = []
        rgp_id = 1  # Initialize rgp_id here

        with tqdm(total=len(coords), desc="Processing coordinates", unit="RGP") as pbar:
            for coord_pair in coords:
                start, end = coord_pair
                # Determine protein indices within the range
                protein_indices = gff_df.index[start:end + 1].intersection(gff_df.index)
                all_protein_indices.extend(protein_indices)
                all_rgp_ids.extend([rgp_id] * len(protein_indices))
                rgp_id += 1
                pbar.update(1)

        # Assign rgp_ids in bulk
        gff_df.loc[all_protein_indices, "rgp_id"] = all_rgp_ids

        # Save the final DataFrame
        gff_df.to_csv(output_file, sep="\t", index=False)
        return gff_df
    
    else:
        logging.info("Output file already exists, using existing file")
        return pd.read_csv(output_file, sep="\t")
    
def classify_rgp(gff_df, genomad, assembly_dir, bakta_dir, plasme_results, vfdb_res, amr_res):
    gff_df["genome_locus"] = gff_df.attributes.map(lambda x: x.split("ID=")[-1].split('_')[0])

    def read_first_line(filename):
        with open(filename, 'r') as file:
            first_line = file.readline().strip()
        return first_line
    
    class genome:
        def __init__(self, fasta_path):
            self.fasta_path = fasta_path
            self.records = read_fasta_no_progress(fasta_path)
            self.id_list = [record.id for record in self.records]
            self.genome_id = fasta_path.split("/")[-1].split(".")[0]
            i=1
            self.contig_convert_dict = {}
            for record in self.records:
                self.contig_convert_dict[f"contig_{i}"] = record.id
                i+=1
    
    genome_locus_dict = {}
    for file in tqdm(glob.glob(f"{bakta_dir}/*/*.faa"), desc = "Getting genome locus", unit="Genome"):
        first_line = read_first_line(file)
        if not ">" in first_line:
            raise ValueError(f"First line of {file} is not a fasta header")
        else:
            genome_locus = first_line.split(">")[1].split("_")[0]
            genome_locus_dict[genome_locus] = file.split("/")[-1].split(".")[0]

    genome_obj_dict = {}
    for fasta_path in tqdm(os.listdir(assembly_dir), desc="Processing genomes", unit="Genome"):
        genome_obj = genome(f"{assembly_dir}/{fasta_path}")
        genome_obj_dict[genome_obj.genome_id] = genome_obj

    def get_contig_true(row):
        genome_id = row.genome_id
        genome_obj = genome_obj_dict[genome_id]
        contig = genome_obj.contig_convert_dict[row.contig]
        return contig

    gff_df["genome_id"] = gff_df.genome_locus.map(lambda x: genome_locus_dict[x])
    gff_df["contig_true"] = gff_df.apply(get_contig_true, axis=1)

    
    rgp_df = gff_df.groupby("rgp_id").agg(
        start=("start", "min"),
        end=("end", "max"),
        contig=("contig", "first"),
        genome_locus=("genome_locus", "first"),
        contig_true=("contig_true", "first"),
        genome_id=("genome_id", "first")
    ).reset_index()

    # columns are: seq_name, length, topology, n_genes, genetic_code, plasmid_score, fdr, n_hallmarks, marker_enrichment, conjugation, genes
    genomad_files = glob.glob(f"{genomad}/*/*_summary/*virus_summary.tsv")

    genomad_df_list = []
    for file in genomad_files:
        df = pd.read_csv(file, sep="\t")
        df["genome_id"] = file.split("/")[-1].split("_virus_summary")[0]
        genomad_df_list.append(df)
    genomad_df = pd.concat(genomad_df_list)
    genomad_df["contig_true"] = genomad_df.seq_name.map(lambda x: x.split("|")[0])
    
    def process_coords(row):
        if "provirus" in row.seq_name:
            row["start"]=int(row.seq_name.split("provirus_")[-1].split("_")[0])
            row["end"]=int(row.seq_name.split("provirus_")[-1].split("_")[1])
        else:
            row["start"]=1
            row["end"]=int(row.length)
        return row
    genomad_df = genomad_df.apply(process_coords, axis=1)

    def viral_checker(row):
        filtered_df = genomad_df[
            (row["contig_true"] == genomad_df["contig_true"]) &  # same contig
            (
                (row["start"] >= genomad_df["start"]) & (row["start"] <= genomad_df["end"]) |  # start of RGP is within viral region
                (row["end"] >= genomad_df["start"]) & (row["end"] <= genomad_df["end"])  # end of RGP is within viral region
            )
        ]
        return len(filtered_df) > 0

    rgp_df["viral"] = rgp_df.apply(viral_checker, axis=1)

    viral_rgp_dict = dict(zip(rgp_df.rgp_id, rgp_df.viral))
    gff_df["viral"] = gff_df.rgp_id.map(lambda x: viral_rgp_dict[x])
    
    
    class plasmid:
        def __init__(self, genome_id, record):
            self.genome_id = genome_id
            self.record = record
            self.contig = record.id
            self.length = len(record.seq)

    manual_plasmid_dict = {}
    plasme_df_list = []

    # Add plasmid results
    for file in glob.glob(f"{plasme_results}/*/*_report.csv"):
        genome_id = file.split("/")[-2].split('.')[0]

        if os.path.getsize(file) == 0:
            assembly_file = glob.glob(f"{assembly_dir}/{genome_id}.fna")
            if len(assembly_file) == 0 or len(assembly_file) > 1:
                raise ValueError(f"ERROR! {len(assembly_file)} assemblies found for {genome_id}: {assembly_file}")
            else:
                assembly_records = read_fasta_no_progress(assembly_file[0])
                assembly_records = sorted(assembly_records, key=lambda record: len(record.seq), reverse=True)
                plasmids = assembly_records[1:]  # Exclude the first (presumably chromosomal) record
                for plasmid_record in plasmids:
                    plasmid_obj = plasmid(genome_id, plasmid_record)
                    manual_plasmid_dict[plasmid_obj.contig] = plasmid_obj
        else:
            df = pd.read_csv(file, sep="\t")
            df["genome_id"] = genome_id
            plasme_df_list.append(df)

    # Concatenate all plasme DataFrames
    plasme_df = pd.concat(plasme_df_list, ignore_index=True)

    # Define the columns for plasme_df
    plasme_df_cols = ["contig", "length", "reference", "order", "evidence", "score", "amb_region"]

    # Append manual plasmid data to plasme_df
    for contig, plasmid_obj in manual_plasmid_dict.items():
        new_row = pd.DataFrame({
            "contig": [plasmid_obj.contig],
            "length": [plasmid_obj.length],
            "reference": ["manual"],
            "order": ["manual"],
            "evidence": ["manual"],
            "score": [1],
            "amb_region": ["manual"]
        }, columns=plasme_df_cols)
        
        plasme_df = pd.concat([plasme_df, new_row], ignore_index=True)
        
    
    gff_df["plasmid"] = gff_df.contig_true.map(lambda x: x in plasme_df.contig.values)
    print("RGP not plasmid associated",len(gff_df[gff_df.plasmid==False]))
    print("RGP plasmid associated",len(gff_df[gff_df.plasmid==True]))
    
    #analyze vfdb and amr
    vfdb_df = pd.read_csv(vfdb_res, sep="\t")
    total_vfdb = len(vfdb_df)
    gff_prot_list = gff_df.protein_id.tolist()
    
    def analyze_vfdb():
        vfdb_df_RGP = vfdb_df[vfdb_df["query"].isin(gff_prot_list)]
        rgp_vfdb_len = len(vfdb_df_RGP)
        
        print(f"Total VFDB: {total_vfdb}")
        print(f"RGP VFDB: {rgp_vfdb_len}")
        
        vfdb_df = pd.read_csv(vfdb_res, sep="\t")
        total_vfdb = len(vfdb_df)

        # Convert gff_df protein IDs to list
        gff_prot_list = gff_df.protein_id.tolist()

        # Filter RGP entries
        vfdb_df_RGP = vfdb_df[vfdb_df["query"].isin(gff_prot_list)]
        rgp_vfdb_len = len(vfdb_df_RGP)

        # Filter non-RGP entries
        vfdb_df_non_RGP = vfdb_df[~vfdb_df["query"].isin(gff_prot_list)]

        # Calculate counts for each VF_cat in both RGP and non-RGP datasets
        rgp_vfdb_counts = vfdb_df_RGP["VF_cat"].value_counts(normalize=True) * 100
        non_rgp_vfdb_counts = vfdb_df_non_RGP["VF_cat"].value_counts(normalize=True) * 100

        # Align both RGP and non-RGP VF_cat counts to have the same categories
        vf_cats = sorted(set(vfdb_df["VF_cat"]))  # Get all VF categories

        rgp_vfdb_counts = rgp_vfdb_counts.reindex(vf_cats, fill_value=0)
        non_rgp_vfdb_counts = non_rgp_vfdb_counts.reindex(vf_cats, fill_value=0)

        # Create double bar plot
        fig, ax = plt.subplots()

        bar_width = 0.35
        index = range(len(vf_cats))

        # Plot bars
        bars_rgp = ax.bar(index, rgp_vfdb_counts, bar_width, label='RGP')
        bars_non_rgp = ax.bar([i + bar_width for i in index], non_rgp_vfdb_counts, bar_width, label='Non-RGP')

        # Labeling
        ax.set_xlabel('VF Category')
        ax.set_ylabel('Percentage of Total Hits (%)')
        ax.set_title('RGP vs Non-RGP Virulence Factor Composition')
        ax.set_xticks([i + bar_width / 2 for i in index])
        ax.set_xticklabels(vf_cats, rotation=45, ha="right")
        ax.legend()

        # Show plot
        plt.tight_layout()
        plt.show()
    #analyze_vfdb()
    
    def analyze_amr():
        amr_df = pd.read_csv(amr_res, sep="\t")

        # Filter for AMR element type only
        amr_df = amr_df[amr_df["Element type"] == "AMR"]
        all_amr = len(amr_df)

        # Convert gff_df protein IDs to list
        gff_prot_list = gff_df.protein_id.tolist()

        # Filter RGP AMR entries
        amr_df_RGP = amr_df[amr_df["Protein identifier"].isin(gff_prot_list)]
        rgp_amr_len = len(amr_df_RGP)

        # Filter non-RGP AMR entries
        amr_df_non_RGP = amr_df[~amr_df["Protein identifier"].isin(gff_prot_list)]

        # Calculate counts for each "Class" in both RGP and non-RGP datasets
        rgp_amr_counts = amr_df_RGP["Class"].value_counts(normalize=True) * 100
        non_rgp_amr_counts = amr_df_non_RGP["Class"].value_counts(normalize=True) * 100

        # Align both RGP and non-RGP Class counts to have the same categories
        amr_classes = sorted(set(amr_df["Class"]))  # Get all unique AMR classes

        rgp_amr_counts = rgp_amr_counts.reindex(amr_classes, fill_value=0)
        non_rgp_amr_counts = non_rgp_amr_counts.reindex(amr_classes, fill_value=0)

        # Create double bar plot
        fig, ax = plt.subplots()

        bar_width = 0.35
        index = range(len(amr_classes))

        # Plot bars
        bars_rgp = ax.bar(index, rgp_amr_counts, bar_width, label='RGP')
        bars_non_rgp = ax.bar([i + bar_width for i in index], non_rgp_amr_counts, bar_width, label='Non-RGP')

        # Labeling
        ax.set_xlabel('AMR Class')
        ax.set_ylabel('Percentage of Total AMR Genes (%)')
        ax.set_title('RGP vs Non-RGP AMR Gene Composition')
        ax.set_xticks([i + bar_width / 2 for i in index])
        ax.set_xticklabels(amr_classes, rotation=45, ha="right")
        ax.legend()

        # Show plot
        plt.tight_layout()
        plt.show()
        
        
    analyze_amr()
    return gff_df
    



@running_message
def UMAP_visualization(gff_df, output):
    colors = [
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", 
    "#e377c2", "#7f7f7f", "#bcbd22", "#17becf", "#1b9e77", "#d95f02", 
    "#7570b3", "#e7298a", "#66a61e", "#e6ab02", "#a6761d", "#666666",
    "#ff6f61", "#6b5b95"
    ]
    umap_csv = f"{output}/umap.csv"
    
    if not os.path.exists(umap_csv):
        rgp_abs_matrix = gff_df.groupby(["rep", "rgp_id"]).size().unstack(fill_value=0).transpose()
        logging.info(rgp_abs_matrix)

        logging.info("UMAP 2D processing")
        proj_2d = UMAP(n_components=2).fit_transform(rgp_abs_matrix)

        umap_df = pd.DataFrame({
            "rgp_id": rgp_abs_matrix.index,
            "x": proj_2d[:,0],
            "y": proj_2d[:,1]
        })
        logging.info("Saving UMAP to csv")
        umap_df.to_csv(umap_csv, index=False)
    else:
        logging.info("UMAP csv already exists, using existing file")
        umap_df = pd.read_csv(umap_csv)

    # Always generate and show the dendrogram
    proj_2d = umap_df[["x", "y"]].values
    Z = linkage(proj_2d, method='ward')
    
    plt.figure(figsize=(10, 7))
    dendrogram(Z)
    plt.title("Dendrogram")
    plt.show()

    # Perform DBSCAN clustering
    logging.info("Using DBSCAN for clustering")
    dbscan = DBSCAN(eps=0.3, min_samples=10)  # Adjust parameters as needed
    cluster_labels = dbscan.fit_predict(proj_2d)

    umap_df["cluster"] = cluster_labels
    umap_df.to_csv(umap_csv, index=False)  # Save the updated DataFrame with clusters

    # Assuming 'plasmid' is already in gff_df, use it to color the plot
    umap_df["plasmid"] = gff_df["plasmid"].reset_index(drop=True)
    umap_df["cluster_color"] = umap_df["cluster"].astype(str)  # Convert cluster labels to strings for color mapping

    # Generate the bar chart showing cluster sizes and plasmid percentages
    cluster_summary = umap_df.groupby('cluster').agg(
        total_count=('plasmid', 'size'),
        plasmid_count=('plasmid', 'sum')
    ).reset_index()

    cluster_summary['plasmid_percentage'] = (cluster_summary['plasmid_count'] / cluster_summary['total_count']) * 100

    # Plot the bar chart
    fig, ax1 = plt.subplots(figsize=(10, 6))

    # Bar chart for total count per cluster
    ax1.bar(cluster_summary['cluster'], cluster_summary['total_count'], color='lightblue', alpha=0.6, label='Total Count')
    ax1.set_xlabel('Cluster')
    ax1.set_ylabel('Total Count')
    ax1.set_title('Cluster Sizes and Plasmid Percentages')
    
    # Line plot for plasmid percentage per cluster
    ax2 = ax1.twinx()
    ax2.plot(cluster_summary['cluster'], cluster_summary['plasmid_percentage'], color='orange', marker='o', linestyle='-', linewidth=2, label='Plasmid Percentage')
    ax2.set_ylabel('Plasmid Percentage (%)')

    # Combine legends from both axes
    ax1.legend(loc='upper left')
    ax2.legend(loc='upper right')

    plt.show()

    # 2D Scatter plot using plotly
    fig2d = px.scatter(
        umap_df, 
        x="x", 
        y="y", 
        opacity=0.7, 
        color="cluster_color",  # Use cluster color for distinct coloring
        symbol="plasmid",       # Differentiate by plasmid/non-plasmid status
        color_discrete_sequence=colors  # Apply custom color palette
    )

    fig2d.show()

def main():
    args = arguments("rgp_tsv", "gff", "gene_fam","output", "genomad", "assembly_dir", "bakta_dir", "plasme_results", "vfdb_res", "amr_res")
    os.makedirs(args.output, exist_ok=True)
    init_logging(f"{args.output}/cluster_rgp.log")
                 
    gff_df = add_rgp(args.gff, args.rgp_tsv, args.gene_fam, args.output)
    logging.info("RGP added to GFF file")
    print(gff_df)
    gff_df = gff_df[gff_df.rgp_id.notnull()]
    logging.info(gff_df)

    gff_df = classify_rgp(gff_df, args.genomad, args.assembly_dir, args.bakta_dir, args.plasme_results, args.vfdb_res, args.amr_res)

    #UMAP_visualization(gff_df, args.output)

if __name__ == "__main__":
    main()