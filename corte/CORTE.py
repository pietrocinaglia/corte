# -*- coding: utf-8 -*-

import os
import requests
import networkx as nx
import pandas as pd
import scipy
import statistics
from itertools import combinations
import matplotlib.pyplot as plt
from joblib import Parallel, delayed

import csv

class CORTE:
    WORKSPACE = os.path.dirname(os.path.realpath(__file__)) + "/"
    AGES = ['20-29', '30-39', '40-49', '50-59', '60-69', '70-79']
    UNIT = 'TPM'
    DATASETID = "gtex_v10"
    LIMIT = 100000
    metadata = []
    temporal_network = list()
    genes_of_interest = list()
    id_type = 'gene_code'
    tissues_of_interest = list()
    threshold = 0.05
    verbose = False

    def __init__(self, genes_of_interest:list, tissues_of_interest:list, threshold:float=0.05, verbose:bool=False):
        if len(genes_of_interest) == 0:
            raise Exception("Genes of interest ('genes_of_interest:list') are mandatory.")

        self.genes_of_interest = genes_of_interest
        self.tissues_of_interest = tissues_of_interest
        self.threshold = threshold
        self.verbose = verbose

        print("[INFO] Loading metadata...")
        df = pd.read_csv(self.WORKSPACE + 'metadata.csv', header=0, sep=" ")
        df = df[df['gene_symbol'].isin(self.genes_of_interest)]
        self.metadata = dict(zip(df.gene_symbol, df.ensembl_id))
        print("- Metadata [OK] ")

    def __retrieve_data(self, action='geneExpression', genes_of_interest=None) -> pd.DataFrame:
        if genes_of_interest is None:
            genes_of_interest = self.genes_of_interest

        if action == 'geneExpression':
            api = "https://gtexportal.org/api/v2/expression/geneExpression"
            params = {
                'itemsPerPage': self.LIMIT,
                'datasetId': self.DATASETID,
                'gencodeId': genes_of_interest,
                'tissueSiteDetailId': self.tissues_of_interest,
                'attributeSubset': 'ageBracket',
                'format': 'json'
            }
            result = requests.get(api, params=params).json()
            return pd.DataFrame(result['data'])
        else:
            raise Exception("The action in Data Retrieving is not supported.")

    def construct_temporal_network(self) -> list:
        if self.verbose:
            print("[INFO] Data Retrieving...")

        gene_symbol_id = self.metadata
        gencodeIds = list(gene_symbol_id.values())
        df = self.__retrieve_data(action='geneExpression', genes_of_interest=gencodeIds)

        gene_pairs = list(combinations(gene_symbol_id, 2))
        temporal_network = []

        if self.verbose:
            print("[INFO] Temporal Network Construction...")

        for i in range(len(self.AGES)):
            age_group = self.AGES[i]
            timepoint = nx.Graph()
            timepoint.add_nodes_from(gene_symbol_id.keys())

            def compute_edge(u, v):
                u_data = df[(df['gencodeId'] == gene_symbol_id[u]) & 
                            (df['tissueSiteDetailId'].isin(self.tissues_of_interest)) &
                            (df['unit'] == self.UNIT) & 
                            (df['subsetGroup'] == age_group)]

                v_data = df[(df['gencodeId'] == gene_symbol_id[v]) & 
                            (df['tissueSiteDetailId'].isin(self.tissues_of_interest)) &
                            (df['unit'] == self.UNIT) & 
                            (df['subsetGroup'] == age_group)]

                u_gexp = list()
                if len(u_data.data) == 1:
                    u_gexp = u_data.data.values[0]
                else:
                    u_gexp = [ statistics.median(gexp_i) if len(gexp_i) > 0 else 0 for gexp_i in u_data.data ]
                
                v_gexp = list()
                if len(v_data.data) == 1:
                    v_gexp = v_data.data.values[0]
                else:
                    v_gexp = [ statistics.median(gexp_i) if len(gexp_i) > 0 else 0 for gexp_i in v_data.data ]

                if ( len(u_gexp) < 3 or len(v_gexp) < 3):
                    if self.verbose:
                        raise Warning("Size of sample less than 3; skipped.")
                    return None

                s, p = scipy.stats.pearsonr(u_gexp, v_gexp)

                return (u, v, p) if p < self.threshold else None

            edges = Parallel(n_jobs=-1)(delayed(compute_edge)(u, v) for u, v in gene_pairs)

            for e in edges:
                if e:
                    u, v, p = e
                    timepoint.add_edge(u, v, pvalue=p)

            temporal_network.append(timepoint)

            if self.verbose:
                print(f"-- Timepoint {age_group} built with {timepoint.number_of_edges()} edges.")

        self.temporal_network = temporal_network
        return temporal_network

    def analyze_temporal_network(self, temporal_network: list) -> pd.DataFrame:
        all_stats = []
        previous_stats = None  # For computing deltas

        for i, G in enumerate(temporal_network):
            stats = {"timepoint": self.AGES[i]}
            n_nodes = G.number_of_nodes()
            n_edges = G.number_of_edges()

            stats["nodes"] = n_nodes
            stats["edges"] = n_edges
            stats["density"] = nx.density(G)
            stats["avg_degree"] = sum(dict(G.degree()).values()) / n_nodes if n_nodes > 0 else 0
            stats["avg_clustering"] = nx.average_clustering(G) if n_nodes > 0 else 0
            stats["components"] = nx.number_connected_components(G) if n_nodes > 0 else 0
            stats["transitivity"] = nx.transitivity(G) if n_nodes > 0 else 0
            stats["assortativity"] = nx.degree_pearson_correlation_coefficient(G) if n_edges > 0 and n_nodes > 1 else None

            # Centralities
            if n_nodes > 1 and n_edges > 0:
                deg_cent = nx.degree_centrality(G)
                bet_cent = nx.betweenness_centrality(G)
                clo_cent = nx.closeness_centrality(G)

                stats["avg_degree_centrality"] = statistics.mean(deg_cent.values())
                stats["avg_betweenness"] = statistics.mean(bet_cent.values())
                stats["avg_closeness"] = statistics.mean(clo_cent.values())

                stats["top_node_by_degree"] = max(deg_cent, key=deg_cent.get)
                stats["top_node_by_betweenness"] = max(bet_cent, key=bet_cent.get)
                stats["top_node_by_closeness"] = max(clo_cent, key=clo_cent.get)
            else:
                stats.update({
                    "avg_degree_centrality": 0,
                    "avg_betweenness": 0,
                    "avg_closeness": 0,
                    "top_node_by_degree": None,
                    "top_node_by_betweenness": None,
                    "top_node_by_closeness": None,
                })

            # Largest connected component stats
            if n_nodes > 0 and n_edges > 0:
                components = sorted(nx.connected_components(G), key=len, reverse=True)
                largest_cc = G.subgraph(components[0])
                if largest_cc.number_of_nodes() > 1:
                    try:
                        stats["diameter"] = nx.diameter(largest_cc)
                        stats["avg_shortest_path_length"] = nx.average_shortest_path_length(largest_cc)
                    except:
                        stats["diameter"] = None
                        stats["avg_shortest_path_length"] = None
                else:
                    stats["diameter"] = None
                    stats["avg_shortest_path_length"] = None
            else:
                stats["diameter"] = None
                stats["avg_shortest_path_length"] = None

            # Temporal deltas from previous timepoint
            if previous_stats:
                for key in ["edges", "density", "avg_degree", "avg_clustering",
                            "avg_degree_centrality", "avg_betweenness", "avg_closeness"]:
                    prev_val = previous_stats.get(key, 0)
                    curr_val = stats.get(key, 0)
                    stats[f"delta_{key}"] = curr_val - prev_val
            else:
                for key in ["edges", "density", "avg_degree", "avg_clustering",
                            "avg_degree_centrality", "avg_betweenness", "avg_closeness"]:
                    stats[f"delta_{key}"] = 0  # First timepoint has no delta

            previous_stats = stats.copy()
            all_stats.append(stats)

        return pd.DataFrame(all_stats)


    def extract_high_degree_genes(self, temporal_network: list, top_n: int = 10) -> dict:
        """
        Extracts the top-N genes with the highest degree for each timepoint.

        Returns:
            dict: { 'age_group': [(gene, degree), ...] }
        """
        top_genes = {}

        for i, G in enumerate(temporal_network):
            age_group = self.AGES[i]

            if G.number_of_nodes() == 0:
                top_genes[age_group] = []
                continue

            degrees = sorted(G.degree, key=lambda x: x[1], reverse=True)
            top_genes[age_group] = degrees[:top_n]

        return top_genes


    def export_adjacency_matrices(self, temporal_network: list, output_path: str):#, nan: int= 1):
        for i, G in enumerate(temporal_network):
            nodes = list(G.nodes())
            # Create adjacency matrix weighted by 'pvalue'
            adj = nx.to_pandas_adjacency(G, nodelist=nodes, weight='pvalue')#.fillna(nan) # if missing edges (NaN) have to replaced,e.g., with nan (input value)
            filename = os.path.join(output_path, f"adjacency_timepoint_{i+1}.csv")
            adj.to_csv(filename)


    def plot(self, temporal_network:list, with_labels:bool=True, output_path:str=None):
        if not temporal_network:
            raise Exception("Temporal Network is empty.")

        for i, tp in enumerate(temporal_network):
            plt.figure(figsize=(8,6))
            plt.title(f"Timepoint {i+1} - {self.AGES[i]}")
            pos = nx.spring_layout(tp)
            degrees = dict(tp.degree())
            node_color = [degrees[n] for n in tp.nodes()]
            nx.draw(tp, pos, node_color=node_color, cmap=plt.cm.viridis, with_labels=with_labels, node_size=300)

            if output_path:
                plt.savefig(output_path + f"timepoint_{i+1}.png", dpi=500)
                plt.clf()
            else:
                plt.show()

    def save_as_files(self, temporal_network:list, output_path:str):
        if not temporal_network:
            raise Exception("Temporal Network is empty.")
        if not output_path:
            raise Exception("You must provide an output path.")

        for i, tp in enumerate(temporal_network):
            tp_title = "timepoint" + str(i+1)
            nx.write_edgelist(tp, output_path + tp_title + '.txt', data=True)
