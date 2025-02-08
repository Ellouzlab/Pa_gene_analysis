import argparse
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from collections import defaultdict, Counter
import sys
from matplotlib.patches import Patch
from matplotlib import colormaps
import numpy as np  # For numerical operations
import os  # For directory operations
from pyvis.network import Network  # For interactive network visualization


def arguments():
    parser = argparse.ArgumentParser(description="Draw a network from a CSV file.")
    parser.add_argument("--taxonomy", required=True, help="Database taxonomy file")
    parser.add_argument("--hits_csv", required=True, help="CSV file with hits")
    parser.add_argument("--id", type=float, default=1, help="Min identity for hit (0-1)")
    parser.add_argument("--cov", type=float, default=1, help="Min coverage for hit (0-1)")
    parser.add_argument("--output", required=True, help="Output folder for images")
    parser.add_argument("--aggregate", help="Aggregation level", default="species",
                        choices=["genus", "family", "species"])
    parser.add_argument("--color", help="Color by", choices=["family", "genus", "species"], default="genus")
    parser.add_argument("--agglomerans", action='store_true',
                        help="Treat Pantoea agglomerans as a separate node and do not aggregate with other Pantoea or Erwinia species")
    parser.add_argument("--inclusive", action='store_true',
                        help="Use all genomes in the taxonomy data for edge weight calculations")
    parser.add_argument("--min_genomes", type=int, default=1,
                        help="Minimum number of genomes to include a node in the network")
    parser.add_argument("--max_nodes", type=int, default=50,
                        help="Maximum number of nodes to include in the network")
    # New arguments for simplification
    parser.add_argument("--simplify", action='store_true',
                        help="Output a simplified version of the graph as well")
    parser.add_argument("--simplify_top_n", type=int, default=50,
                        help="Number of top nodes to include in the simplified graph based on centrality")
    parser.add_argument("--simplify_min_edge_weight", type=float, default=0.0,
                        help="Minimum edge weight to include in the simplified graph")
    args = parser.parse_args()

    if args.id < 0 or args.id > 1 or args.cov < 0 or args.cov > 1:
        parser.error("Identity and coverage thresholds must be between 0 and 1.")

    return args

def get_taxonomy_counts(taxonomy_path):
    # Read the taxonomy CSV
    taxonomy_df = pd.read_csv(taxonomy_path)
    # Create a mapping from GCA to taxonomy
    taxonomy_dict = {}
    for idx, row in taxonomy_df.iterrows():
        gca = row['GCA']
        taxonomy_dict[gca] = {
            'Family': row['Family'],
            'Genus': row['Genus'],
            'Species': row['Species'],
            'Protein_Count': row.get('Protein_Count', 0)
        }
    return taxonomy_df, taxonomy_dict


def filter_hits(hits_csv, id_threshold, cov_threshold):
    # Read the hits CSV
    hits_df = pd.read_csv(hits_csv)
    # Filter hits based on identity and coverage thresholds
    hits_filtered = hits_df[
        (hits_df['identity'] >= id_threshold * 100) &
        (hits_df['query_coverage'] >= cov_threshold * 100)
    ]
    # Remove duplicate query-to-subject pairs
    hits_filtered = hits_filtered.drop_duplicates(subset=['query_id', 'subject_id'])
    return hits_filtered

def aggregate_graph(hits_df, aggregation_level, taxonomy_dict, agglomerans_option, inclusive_option, min_genomes=3):
    query_hits = defaultdict(set)
    node_genomes = defaultdict(set)  # Map node to set of genomes with hits
    all_node_genomes = defaultdict(set)  # Map node to all genomes in taxonomy
    node_labels = {}  # Map node to label for display

    # Build all_node_genomes mapping from taxonomy data
    for gca, taxonomy in taxonomy_dict.items():
        if agglomerans_option and taxonomy['Genus'] == 'Pantoea' and taxonomy['Species'] == 'agglomerans':
            node = 'Pantoea agglomerans'
        else:
            if aggregation_level == 'species':
                if taxonomy['Species'] == 'sp.':
                    continue  # Exclude 'sp.' species
                node = taxonomy['Species']
            elif aggregation_level == 'genus':
                node = taxonomy['Genus']
            elif aggregation_level == 'family':
                node = taxonomy['Family']
            else:
                print(f"Unknown aggregation level: {aggregation_level}")
                sys.exit(1)
        all_node_genomes[node].add(gca)

    # Process hits
    for idx, row in hits_df.iterrows():
        query_id = row['query_id']
        subject_id = row['subject_id']

        # Parse the subject_id to extract GCA and taxonomy
        parts = subject_id.split('___')
        if len(parts) >= 4:
            gca = parts[0]
            family = parts[1]
            genus = parts[2]
            species = parts[3]
        else:
            print(f"Warning: subject_id '{subject_id}' does not have the expected format.")
            continue

        # Get the taxonomy info from the taxonomy_dict
        taxonomy = taxonomy_dict.get(gca)
        if not taxonomy:
            print(f"Warning: GCA '{gca}' not found in taxonomy data.")
            continue

        # Determine the node based on aggregation level and agglomerans_option
        if agglomerans_option and taxonomy['Genus'] == 'Pantoea' and taxonomy['Species'] == 'agglomerans':
            node = 'Pantoea agglomerans'
            label = 'P. agglomerans'
        else:
            if aggregation_level == 'species':
                if taxonomy['Species'] == 'sp.':
                    continue  # Exclude 'sp.' species
                node = taxonomy['Species']
                label = f"{taxonomy['Genus']} {taxonomy['Species']}"
            elif aggregation_level == 'genus':
                node = taxonomy['Genus']
                label = taxonomy['Genus']
            elif aggregation_level == 'family':
                node = taxonomy['Family']
                label = taxonomy['Family']
            else:
                print(f"Unknown aggregation level: {aggregation_level}")
                sys.exit(1)

        query_hits[query_id].add(node)
        node_genomes[node].add(gca)  # Only genomes with hits are added
        if node not in node_labels:
            node_labels[node] = label

    if aggregation_level == 'family':
        # For family level, always use all_node_genomes for filtering
        nodes_to_keep = {node for node, genomes in all_node_genomes.items() if len(genomes) >= min_genomes}
    else:
        if inclusive_option:
            # Use all_node_genomes for filtering
            nodes_to_keep = {node for node, genomes in all_node_genomes.items() if len(genomes) >= min_genomes}
        else:
            # Use node_genomes (genomes with hits) for filtering
            nodes_to_keep = {node for node, genomes in node_genomes.items() if len(genomes) >= min_genomes}

    # Filter node_genomes and node_labels
    node_genomes = {node: genomes for node, genomes in node_genomes.items() if node in nodes_to_keep}
    node_labels = {node: label for node, label in node_labels.items() if node in nodes_to_keep}

    # Build the graph with filtered nodes
    G = nx.Graph()
    G.add_nodes_from(node_genomes.keys())

    # Build edges based on shared hits to the same query_id
    edge_weights = defaultdict(int)
    for query_id, nodes in query_hits.items():
        nodes = [node for node in nodes if node in node_genomes]  # Only consider filtered nodes
        for i in range(len(nodes)):
            for j in range(i + 1, len(nodes)):
                node1 = nodes[i]
                node2 = nodes[j]
                if node1 > node2:
                    node1, node2 = node2, node1
                edge_weights[(node1, node2)] += 1

    # Calculate edge weights
    for (node1, node2), shared_hits in edge_weights.items():
        if inclusive_option:
            # Use all genomes from taxonomy data
            num_genomes_node1 = len(all_node_genomes[node1])
            num_genomes_node2 = len(all_node_genomes[node2])
        else:
            # Use only genomes with hits
            num_genomes_node1 = len(node_genomes[node1])
            num_genomes_node2 = len(node_genomes[node2])
        total_genomes = num_genomes_node1 + num_genomes_node2
        weight = shared_hits / total_genomes if total_genomes > 0 else 0
        G.add_edge(node1, node2, weight=weight)

    return G, node_genomes, node_labels



def color_nodes(G, color_by, node_genomes, taxonomy_dict, agglomerans_option):
    # Assign colors to nodes based on the 'color_by' parameter
    node_attributes = {}
    for node in G.nodes():
        genomes = node_genomes.get(node, set())
        attrs = []
        for gca in genomes:
            taxonomy = taxonomy_dict.get(gca)
            if taxonomy:
                # Handle Pantoea agglomerans separately
                if agglomerans_option and node == 'Pantoea agglomerans':
                    attrs.append('Pantoea agglomerans')
                else:
                    attrs.append(taxonomy[color_by.capitalize()])

        if attrs:
            # Assign the most common attribute
            most_common_attr = Counter(attrs).most_common(1)[0][0]
            node_attributes[node] = most_common_attr

    # Assign colors
    unique_attrs = set(node_attributes.values())
    unique_attrs.discard('Unknown')
    unique_attrs = sorted(unique_attrs)
    attr_to_color = {}

    # Use a discrete colormap
    base_cmap = colormaps.get_cmap('hsv')
    num_colors = len(unique_attrs)
    if num_colors == 0:
        print("No unique attributes found for coloring.")
        return [], {}
    color_list = [base_cmap(i / num_colors) for i in range(num_colors)]
    for idx, attr in enumerate(unique_attrs):
        # Convert RGBA to HEX
        rgba_color = color_list[idx]
        hex_color = '#%02x%02x%02x' % tuple(int(c * 255) for c in rgba_color[:3])
        attr_to_color[attr] = hex_color

    node_colors = []
    for node in G.nodes():
        attr = node_attributes.get(node, 'Unknown')
        color = attr_to_color.get(attr, '#000000')  # Default to black
        node_colors.append(color)

    return node_colors, attr_to_color


def draw_network(G, node_colors, attr_to_color, output_base, color_attr, aggregation_level, node_labels):
    plt.figure(figsize=(6, 6))
    try:
        # Use weights in the layout
        pos = nx.spring_layout(G, k=0.15, iterations=20, weight='weight')
    except Exception as e:
        print(f"Error computing layout: {e}")
        sys.exit(1)

    # Invert weights for betweenness centrality
    for u, v, data in G.edges(data=True):
        data['inv_weight'] = 1 / data['weight'] if data['weight'] != 0 else float('inf')

    # Compute eigenvector centrality with weights
    eigen_centrality = nx.eigenvector_centrality(G, max_iter=1000, weight='weight')

    # Compute betweenness centrality with inverted weights
    betweenness_centrality = nx.betweenness_centrality(G, weight='inv_weight')

    # Normalize centrality measures for visualization
    ec_values = np.array(list(eigen_centrality.values()))
    ec_min, ec_max = ec_values.min(), ec_values.max()
    if ec_max > ec_min:
        ec_norm_values = (ec_values - ec_min) / (ec_max - ec_min)
    else:
        ec_norm_values = ec_values

    # Create a dict mapping node to normalized centrality
    ec_norm_dict = {node: ec_norm_values[idx] for idx, node in enumerate(eigen_centrality.keys())}

    # Adjust node sizes based on centrality
    node_sizes = [5 + 100 * ec_norm_dict[node] for node in G.nodes()]  # Scale node sizes

    # Draw nodes with sizes based on eigenvector centrality
    nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=node_sizes)

    # Draw edges
    edges = G.edges()
    weights = [G[u][v]['weight'] for u, v in edges]
    max_weight = max(weights) if weights else 1
    edge_widths = [(w / max_weight) * 2 for w in weights]  # Edge width scaling
    nx.draw_networkx_edges(G, pos, edge_color='grey', alpha=0.5, width=edge_widths)

    # Add node labels (names on nodes)
    if aggregation_level == 'species':
        labels = node_labels
        font_size = 5  # Decreased font size for species level
    else:
        labels = {node: node_labels.get(node, node) for node in G.nodes()}
        font_size = 7  # Decreased font size
    nx.draw_networkx_labels(G, pos, labels, font_size=font_size)

    # Create color legend
    legend_handles = [Patch(facecolor=color, label=attr) for attr, color in attr_to_color.items()]
    plt.legend(handles=legend_handles, title=color_attr.capitalize(), loc='best', fontsize='small')

    plt.axis('off')
    plt.tight_layout()
    plt.savefig(output_base + '.svg', format='svg')
    plt.close()
    print(f"Network graph saved to {output_base}.svg")

    # Create centrality bar charts
    create_centrality_bar_charts(G, eigen_centrality, betweenness_centrality, node_labels, output_base)

    # Create interactive network graph using pyvis
    create_interactive_network(G, node_labels, output_base, node_colors, attr_to_color, ec_norm_dict)


def select_top_nodes_by_edge_weight(G, node_labels, top_n=50):
    # Calculate summed edge weights for each node
    node_edge_weight_sum = {}
    for node in G.nodes():
        total_weight = sum([G[node][neighbor]['weight'] for neighbor in G.neighbors(node)])
        node_edge_weight_sum[node] = total_weight

    # Select top N nodes by summed edge weights
    top_nodes = sorted(node_edge_weight_sum.items(), key=lambda x: x[1], reverse=True)[:top_n]
    top_nodes_set = set([node for node, weight in top_nodes])

    # Create subgraph with top N nodes
    G_sub = G.subgraph(top_nodes_set).copy()

    # Update node_labels to include only labels of nodes in subgraph
    node_labels_sub = {node: label for node, label in node_labels.items() if node in top_nodes_set}

    return G_sub, node_labels_sub


def create_centrality_bar_charts(G, eigen_centrality, betweenness_centrality, node_labels, output_base):
    # Create DataFrame for centrality measures
    centrality_df = pd.DataFrame({
        'Node': list(G.nodes()),
        'Eigenvector Centrality': list(eigen_centrality.values()),
        'Betweenness Centrality': list(betweenness_centrality.values())
    })

    # Map node labels
    centrality_df['Label'] = centrality_df['Node'].map(node_labels)

    # Sort by Eigenvector Centrality and select top 20
    centrality_df_ec = centrality_df.sort_values('Eigenvector Centrality', ascending=False).head(20)

    # Plot Eigenvector Centrality
    plt.figure(figsize=(10, 6))
    plt.bar(centrality_df_ec['Label'], centrality_df_ec['Eigenvector Centrality'], color='skyblue')
    plt.xticks(rotation=90, fontsize=8)  # Decreased font size
    plt.yticks(fontsize=8)  # Decreased font size
    plt.xlabel('Node', fontsize=10)
    plt.ylabel('Eigenvector Centrality', fontsize=10)
    plt.title('Top 20 Nodes by Eigenvector Centrality', fontsize=12)
    plt.tight_layout()
    eigen_output = output_base + '_eigenvector_centrality.svg'
    plt.savefig(eigen_output, format='svg')
    plt.close()
    print(f"Eigenvector centrality chart saved to {eigen_output}")

    # Sort by Betweenness Centrality and select top 20
    centrality_df_bc = centrality_df.sort_values('Betweenness Centrality', ascending=False).head(20)

    # Plot Betweenness Centrality
    plt.figure(figsize=(10, 6))
    plt.bar(centrality_df_bc['Label'], centrality_df_bc['Betweenness Centrality'], color='salmon')
    plt.xticks(rotation=90, fontsize=8)  # Decreased font size
    plt.yticks(fontsize=8)  # Decreased font size
    plt.xlabel('Node', fontsize=10)
    plt.ylabel('Betweenness Centrality', fontsize=10)
    plt.title('Top 20 Nodes by Betweenness Centrality', fontsize=12)
    plt.tight_layout()
    betweenness_output = output_base + '_betweenness_centrality.svg'
    plt.savefig(betweenness_output, format='svg')
    plt.close()
    print(f"Betweenness centrality chart saved to {betweenness_output}")


def simplify_graph(G, top_n, min_edge_weight):
    # Simplify the graph by keeping top N nodes based on centrality and edges above a weight threshold
    eigen_centrality = nx.eigenvector_centrality(G, max_iter=1000, weight='weight')

    # Select top N nodes
    top_nodes = sorted(eigen_centrality.items(), key=lambda x: x[1], reverse=True)[:top_n]
    top_nodes_set = set(node for node, centrality in top_nodes)

    # Create subgraph with top N nodes
    G_top = G.subgraph(top_nodes_set).copy()

    # Remove edges below the minimum edge weight
    edges_to_remove = [(u, v) for u, v, data in G_top.edges(data=True) if data['weight'] < min_edge_weight]
    G_top.remove_edges_from(edges_to_remove)

    # Optionally, remove isolated nodes
    isolated_nodes = list(nx.isolates(G_top))
    G_top.remove_nodes_from(isolated_nodes)

    return G_top

def plot_agglomerans_edges(G, output_base, node_labels):
    # Check if Pantoea agglomerans is in the graph
    if 'Pantoea agglomerans' in G.nodes:
        agglomerans_node = 'Pantoea agglomerans'
    elif 'agglomerans' in G.nodes:
        agglomerans_node = 'agglomerans'
    else:
        print("Pantoea agglomerans is not present in the network graph.")
        return

    agglomerans_edges = []
    for neighbor in G.neighbors(agglomerans_node):
        edge_data = G.get_edge_data(agglomerans_node, neighbor)
        weight = edge_data.get('weight', None)
        if weight is not None:
            agglomerans_edges.append((neighbor, weight))
        else:
            print(f"Warning: Edge between '{agglomerans_node}' and '{neighbor}' does not have a 'weight' attribute.")

    if not agglomerans_edges:
        print(f"No edges with 'weight' attribute found for '{agglomerans_node}'.")
        return

    # Sort by weight in descending order
    agglomerans_edges = sorted(agglomerans_edges, key=lambda x: x[1], reverse=True)

    # Extract nodes and weights for the plot
    nodes, weights = zip(*agglomerans_edges)

    # Map nodes to full Genus species names using node_labels
    labels = [node_labels.get(node, node) for node in nodes]

    # Plotting the bar graph
    plt.figure(figsize=(10, 6))
    plt.bar(labels, weights, color='teal')
    plt.xticks(rotation=90, fontsize=8)
    plt.yticks(fontsize=8)
    plt.xlabel('Node', fontsize=10)
    plt.ylabel('Edge Weight', fontsize=10)
    plt.title('Edge Weights from Pantoea agglomerans to Other Nodes', fontsize=12)
    plt.tight_layout()

    # Save the plot
    agglomerans_output = output_base + '_agglomerans_edges.svg'
    plt.savefig(agglomerans_output, format='svg')
    plt.close()
    print(f"Bar graph of edge weights from Pantoea agglomerans saved to {agglomerans_output}")

def create_interactive_network(G, node_labels, output_base, node_colors, attr_to_color, ec_norm_dict):
    # Create a pyvis Network instance
    net = Network(height="750px", width="100%", bgcolor="#ffffff", font_color="black")
    net.barnes_hut()  # Set physics model for layout
    net.show_buttons(filter_=['nodes', 'edges', 'physics'])

    # Update node attributes
    for idx, node in enumerate(G.nodes()):
        node_id = node
        label = node_labels.get(node_id, str(node_id))
        color = node_colors[idx]
        # Scale node size by eigenvector centrality
        # Let's scale sizes between min_size and max_size
        min_size = 10
        max_size = 50
        ec_norm_value = ec_norm_dict[node_id]
        size = min_size + (max_size - min_size) * ec_norm_value
        net.add_node(node_id, label=label, title=label, color=color, size=size)

    # Get edge weights
    weights = [data.get('weight', 1) for u, v, data in G.edges(data=True)]
    if weights:
        min_weight = min(weights)
        max_weight = max(weights)
    else:
        min_weight = max_weight = 1

    # Add edges
    for u, v, data in G.edges(data=True):
        weight = data.get('weight', 1)
        # Scale edge width and opacity based on weight
        # Edge width between min_width and max_width
        min_width = 0.2
        max_width = 10
        # Edge opacity between min_opacity and max_opacity
        min_opacity = 0.3
        max_opacity = 1.0
        if max_weight > min_weight:
            normalized_weight = (weight - min_weight) / (max_weight - min_weight)
        else:
            normalized_weight = 1.0
        width = min_width + (max_width - min_width) * normalized_weight
        opacity = min_opacity + (max_opacity - min_opacity) * normalized_weight
        # Set edge color with opacity
        gray_color = f'rgba(128, 128, 128, {opacity})'
        net.add_edge(u, v, value=weight, width=width, color=gray_color, smooth=False)

    # Set options to have straight edges
    net.set_edge_smooth('straight')

    # Save the interactive network graph
    html_output = output_base + '.html'
    net.show(html_output, notebook=False)  # Set notebook=False here
    print(f"Interactive network graph saved to {html_output}")

def main():
    args = arguments()

    # Check if output folder exists, if not, create it
    if not os.path.exists(args.output):
        os.makedirs(args.output)

    taxonomy_df, taxonomy_dict = get_taxonomy_counts(args.taxonomy)
    hits_df = filter_hits(args.hits_csv, args.id, args.cov)
    G_full, node_genomes, node_labels_full = aggregate_graph(
        hits_df, args.aggregate, taxonomy_dict, args.agglomerans, args.inclusive, args.min_genomes
    )

    # Select top nodes by summed edge weights
    G, node_labels = select_top_nodes_by_edge_weight(G_full, node_labels_full, args.max_nodes)

    # Update node_genomes to include only nodes in the subgraph
    node_genomes_sub = {node: node_genomes[node] for node in G.nodes() if node in node_genomes}

    node_colors, attr_to_color = color_nodes(G, args.color, node_genomes_sub, taxonomy_dict, args.agglomerans)

    output_base = os.path.join(args.output, 'network')
    draw_network(G, node_colors, attr_to_color, output_base, args.color, args.aggregate, node_labels)

    # Pass node_labels to plot_agglomerans_edges
    plot_agglomerans_edges(G, output_base, node_labels)

    # Output simplified graph if requested
    if args.simplify:
        print("Generating simplified graph...")
        G_simplified = simplify_graph(G, args.simplify_top_n, args.simplify_min_edge_weight)
        simplified_output_base = os.path.join(args.output, 'network_simplified')

        # Update node_labels and node_colors for the simplified graph
        node_labels_simplified = {node: node_labels[node] for node in G_simplified.nodes()}
        node_genomes_simplified = {node: node_genomes_sub[node] for node in G_simplified.nodes()}
        node_colors_simplified, attr_to_color_simplified = color_nodes(
            G_simplified, args.color, node_genomes_simplified, taxonomy_dict, args.agglomerans)

        draw_network(G_simplified, node_colors_simplified, attr_to_color_simplified,
                     simplified_output_base, args.color, args.aggregate, node_labels_simplified)
        print(f"Simplified network graph saved to {simplified_output_base}.svg")

        # Create interactive network for simplified graph
        create_interactive_network(G_simplified, node_labels_simplified, simplified_output_base, node_colors_simplified, attr_to_color_simplified, ec_norm_dict)


if __name__ == '__main__':
    main()
