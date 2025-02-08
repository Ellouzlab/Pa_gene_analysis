import argparse
import csv
import networkx as nx
import matplotlib.pyplot as plt
from collections import defaultdict
import matplotlib.cm as cm
import numpy as np
from matplotlib.patches import Patch
import sys

def main():
    parser = argparse.ArgumentParser(description='Create a network graph from DIAMOND output CSV.')
    parser.add_argument('--input_csv', required=True, help='Input CSV file from DIAMOND output.')
    parser.add_argument('--output_image', required=True, help='Output image file for the graph (e.g., graph.png).')
    parser.add_argument('--identity_threshold', type=float, default=80.0, help='Minimum identity percentage to include a hit.')
    parser.add_argument('--coverage_threshold', type=float, default=80.0, help='Minimum coverage percentage to include a hit.')
    parser.add_argument('--aggregate_by_genus', action='store_true', help='Aggregate nodes by genus instead of genome.')
    args = parser.parse_args()

    # Dictionaries to store genome information and query hits
    genome_info = {}  # Maps genome_id to its family for coloring
    query_hits = defaultdict(set)  # Maps query_id to a set of genome_ids or genera
    genus_to_genomes = defaultdict(set)  # Maps genus to set of genomes

    # Read the CSV file and parse the data
    with open(args.input_csv, 'r') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            try:
                # Apply filtering based on identity and coverage
                identity = float(row['identity'])
                coverage = float(row['query_coverage'])
                if identity < args.identity_threshold or coverage < args.coverage_threshold:
                    continue  # Skip this hit

                query_id = row['query_id']
                subject_id = row['subject_id']
                # Parse the subject_id to extract genome_id, family, genus, species
                parts = subject_id.split('___')
                if len(parts) >= 4:
                    genome_id = parts[0]
                    family = parts[1]
                    genus = parts[2]
                    species = parts[3]
                    # Store genome information
                    genome_info[genome_id] = {'family': family, 'genus': genus, 'species': species}
                    genus_to_genomes[genus].add(genome_id)
                    # Add genus to the set of genera that hit the query_id
                    if args.aggregate_by_genus:
                        query_hits[query_id].add(genus)
                    else:
                        query_hits[query_id].add(genome_id)
                else:
                    print(f"Warning: subject_id '{subject_id}' does not have the expected format.")
            except ValueError:
                print(f"Warning: Unable to parse identity or coverage in row: {row}")
                continue

    # Build the graph
    G = nx.Graph()
    if args.aggregate_by_genus:
        # Nodes are genera
        all_genera = set()
        for genus_set in query_hits.values():
            all_genera.update(genus_set)
        for genus in all_genera:
            G.add_node(genus)
        # Keep track of pairwise counts
        pairwise_counts = defaultdict(int)
        # Count co-occurrences and build edge weights
        for genera in query_hits.values():
            genera = list(genera)
            for i in range(len(genera)):
                for j in range(i + 1, len(genera)):
                    genus1 = genera[i]
                    genus2 = genera[j]
                    # Ensure consistent ordering
                    if genus1 > genus2:
                        genus1, genus2 = genus2, genus1
                    pairwise_counts[(genus1, genus2)] += 1
        # Add edges with weights adjusted by number of unique genomes
        for (genus1, genus2), shared_hits in pairwise_counts.items():
            num_genomes_genus1 = len(genus_to_genomes[genus1])
            num_genomes_genus2 = len(genus_to_genomes[genus2])
            total_genomes = num_genomes_genus1 + num_genomes_genus2
            edge_weight = shared_hits / total_genomes if total_genomes > 0 else 0
            G.add_edge(genus1, genus2, weight=edge_weight)
    else:
        # Nodes are genome_ids with family as an attribute
        for genome_id, info in genome_info.items():
            G.add_node(genome_id, family=info['family'])
        # Add edges between genomes that have hits to the same query protein
        for genomes in query_hits.values():
            genomes = list(genomes)
            for i in range(len(genomes)):
                for j in range(i + 1, len(genomes)):
                    G.add_edge(genomes[i], genomes[j])

    # Assign colors to nodes based on family
    if args.aggregate_by_genus:
        # Since nodes are genera, map genus to family
        genus_to_family = {}
        for info in genome_info.values():
            genus_to_family[info['genus']] = info['family']
        families = set(genus_to_family[genus] for genus in G.nodes() if genus in genus_to_family)
        family_list = sorted(families)
        num_families = len(family_list)
        colormap = cm.get_cmap('hsv', num_families)
        family_to_color = {family: colormap(idx / num_families) for idx, family in enumerate(family_list)}
        node_colors = []
        for node in G.nodes():
            family = genus_to_family.get(node, 'Unknown')
            color = family_to_color.get(family, (0, 0, 0, 1))  # Default to black if family not found
            node_colors.append(color)
    else:
        families = set(info['family'] for info in genome_info.values())
        family_list = sorted(families)
        family_color_map = {family: idx for idx, family in enumerate(family_list)}
        num_families = len(family_list)
        colormap = cm.get_cmap('hsv', num_families)
        family_to_color = {family: colormap(idx / num_families) for family, idx in family_color_map.items()}
        node_colors = [family_to_color[G.nodes[node]['family']] for node in G.nodes()]

    # Draw the graph
    plt.figure(figsize=(12, 12))
    try:
        pos = nx.spring_layout(G, k=0.15, iterations=20)
    except Exception as e:
        print(f"Error computing layout: {e}")
        sys.exit(1)
    # Increase node size
    node_size = 100 if args.aggregate_by_genus else 50
    nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=node_size)
    # Edge widths based on weights
    if args.aggregate_by_genus:
        edges = G.edges()
        weights = [G[u][v]['weight'] for u, v in edges]
        max_weight = max(weights) if weights else 1
        edge_widths = [(w / max_weight) * 5 for w in weights]  # Scale to max width 5
        nx.draw_networkx_edges(G, pos, edge_color='grey', alpha=0.5, width=edge_widths)
    else:
        nx.draw_networkx_edges(G, pos, alpha=0.5, width=0.5)
    # Display genus names on nodes
    if args.aggregate_by_genus:
        labels = {node: node for node in G.nodes()}
        nx.draw_networkx_labels(G, pos, labels, font_size=8)
    # Optionally, add labels for genome nodes (can be cluttered for large graphs)
    # else:
    #     nx.draw_networkx_labels(G, pos, font_size=8)

    # Create legend for families
    legend_elements = [Patch(facecolor=family_to_color[family], label=family) for family in family_list]
    plt.legend(handles=legend_elements, title='Family', loc='best', fontsize='small')

    plt.axis('off')
    plt.tight_layout()
    plt.savefig(args.output_image, dpi=300)
    plt.close()

    print(f"Graph saved to {args.output_image}")

if __name__ == '__main__':
    main()
