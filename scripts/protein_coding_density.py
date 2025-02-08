import pandas as pd
from utils import *
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
import numpy as np
import seaborn as sns
from scipy.stats import f_oneway
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from statsmodels.stats.multicomp import MultiComparison


def fasta_list_to_dict(fasta_list, type):
    return {record.id: type for record in fasta_list}

def make_histo_obj(ax, df, type, colour):
    data = df[df["gene_type"] == type].length
    ax.hist(data, weights=np.ones(len(data)) / len(data), bins=500, alpha=0.25, label=type, 
            color=colour, histtype='stepfilled', edgecolor=colour, linewidth=1.5)

def plot_histo(df):
    _, ax = plt.subplots()

    make_histo_obj(ax, df, "persistent", "blue")
    make_histo_obj(ax, df, "shell", "green")
    make_histo_obj(ax, df, "cloud", "red")

    plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
    plt.show()

def plot_boxplot(df):
    plt.figure(figsize=(10, 6))
    sns.boxplot(
        x="gene_type", 
        y="length", 
        data=df, 
        palette={"persistent": "blue", "shell": "green", "cloud": "red"},
        showfliers=False  # This hides the outliers in the boxplot
    )
    plt.title("Comparison of Gene Lengths by Gene Type (Outliers Removed)")
    plt.xlabel("Gene Type")
    plt.ylabel("Gene Length")
    plt.show()

def perform_statistical_tests(df):
    # Extract gene lengths by type
    persistent_lengths = df[df["gene_type"] == "persistent"].length
    shell_lengths = df[df["gene_type"] == "shell"].length
    cloud_lengths = df[df["gene_type"] == "cloud"].length

    # Calculate mean lengths
    mean_persistent = persistent_lengths.mean()
    mean_shell = shell_lengths.mean()
    mean_cloud = cloud_lengths.mean()

    # ANOVA to test all groups together
    anova_result = f_oneway(persistent_lengths, shell_lengths, cloud_lengths)

    print(f"ANOVA result: p-value = {anova_result.pvalue}\n")

    # Perform Tukey's HSD test
    mc = MultiComparison(df['length'], df['gene_type'])
    tukey_result = mc.tukeyhsd()

    # Print the Tukey HSD test results
    print(tukey_result)

    # Plot Tukey's HSD results
    plt.figure(figsize=(10, 6))
    tukey_result.plot_simultaneous(comparison_name='persistent', xlabel="Mean Difference", ylabel="Gene Type")
    plt.title("Tukey's HSD Test Results")
    plt.show()

def main():
    args = arguments("gff", "persistent", "cloud", "shell")
    columns = ["contig", "type", "source", "start", "end", "score", "strand", "frame", "attributes"]
    df = pd.read_csv(args.gff, header=None, names=columns, index_col=False)

    df["length"] = df.end - df.start
    df["contig_length"] = pd.NA
    df.loc[df["source"] == "region", "contig_length"] = df["length"]

    df["contig_length"] = df["contig_length"].ffill()
    df = df[df.source == "CDS"]
    
    persistent_fasta = fasta_list_to_dict(read_fasta(args.persistent), "persistent")
    cloud_fasta = fasta_list_to_dict(read_fasta(args.cloud), "cloud")
    shell_fasta = fasta_list_to_dict(read_fasta(args.shell), "shell")

    gene_type_dict = {**persistent_fasta, **cloud_fasta, **shell_fasta}

    df["gene_id"] = df.attributes.map(
        lambda x: x.split("ID=")[1].split(";")[0] if "ID=" in x else x.split("gene_id ")[1].split(";")[0]
    )
    
    df["gene_type"] = df.gene_id.map(lambda x: gene_type_dict[x])

    plot_histo(df)
    plot_boxplot(df)
    perform_statistical_tests(df)

if __name__ == "__main__":
    main()
