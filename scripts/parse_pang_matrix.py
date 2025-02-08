import pandas as pd
import matplotlib.pyplot as plt
from utils import arguments

def plot_gene_distribution(percents_dict, genome_columns):
    import matplotlib.pyplot as plt

    # Assuming 'percents_dict' is already calculated and contains the percentages
    # for 'persistent', 'cloud', and 'shell' categories across different genomes.

    fig, ax = plt.subplots()

    # Plotting each category in a different color
    ax.hist(percents_dict['persistent'], bins=25, alpha=0.5, label='Persistent', color='blue')
    ax.hist(percents_dict['cloud'], bins=25, alpha=0.5, label='Cloud', color='red')
    ax.hist(percents_dict['shell'], bins=25, alpha=0.5, label='Shell', color='green')

    # Adding labels and title
    ax.set_xlabel('Percentage of Genes')
    ax.set_ylabel('Frequency')
    ax.set_title('Percentage Composition of Genes Across Genomes')
    ax.legend()

    # Show the plot
    plt.show()

if __name__ == '__main__':
    args = arguments("input")
    
    df = pd.read_csv(args.input)
    print(df.columns.to_list())
    
    non_genome_column = ['Gene', 'Non-unique Gene name', 'Annotation', 'No. isolates', 'No. sequences', 'Avg sequences per isolate', 
                         'Accessory Fragment', 'Genome Fragment', 'Order within Fragment', 'Accessory Order with Fragment', 
                         'QC', 'Min group size nuc', 'Max group size nuc', 'Avg group size nuc']
    
    gene_names = df['Gene']
    all_columns = df.columns
    genome_columns = all_columns[df.columns.map(lambda x: x not in non_genome_column)]
    
    for col in genome_columns:
        df[col] = df[col].map(lambda x: 0 if pd.isna(x) else 1)
    
    df_genes = df[genome_columns]
    gene_presence_percentage = (df_genes.sum(axis=1) / df_genes.shape[1]) * 100  # Percentage
    
    types_ = df['Non-unique Gene name']
    
    # Identify genes present in all genomes
    persistent_index = df[types_=='persistent'].index
    cloud_index = df[types_=='cloud'].index
    shell_index = df[types_=='shell'].index
    
    indexes = {'persistent': persistent_index, "cloud": cloud_index, "shell": shell_index}
    
    # Calculate the percentage of each genome composed of fully present genes
    percents_dict = {}
    for key, index in indexes.items():
        percents_dict[key] = []
        for column in genome_columns:
            num_genes = df.loc[index, column].sum()
            total_genes_in_genome = df[column].sum()
            percent = (num_genes / total_genes_in_genome) * 100
            percents_dict[key].append(percent)
    
    plot_gene_distribution(percents_dict, genome_columns)
    
    for key, value in percents_dict.items():
        print(f"{key} - mean: {sum(value) / len(value):.2f}, min: {min(value):.2f}, max: {max(value):.2f}, std: {pd.Series(value).std():.2f}")
