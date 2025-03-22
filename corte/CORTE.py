# -*- coding: utf-8 -*-

import requests
import networkx as nx
import pandas as pd
import scipy
import statistics
from itertools import combinations
import matplotlib.pyplot as plt

###
# @author: Pietro Cinaglia
# @mail: cinaglia@unicz.it
# @description: An open source and user-friendly tool for COnstructing Real-world TEmporal networks from genotype-tissue expression data (CoRTE)
# @url: https://github.com/pietrocinaglia/corte
###

class CORTE:
    #
    AGES = ['20-29', '30-39', '40-49', '50-59', '60-69', '70-79']
    UNIT = 'TPM'
    DATASETID = "gtex_v10"
    LIMIT = 100000 # Max is 100000
    #
    temporal_network = list()
    genes_of_interest = list()
    id_type = 'gene_code' # 'gene_code' or 'gene_symbol'
    tissues_of_interest = list()
    threshold = 0.05
    verbose = False
    #

    def __init__(self, genes_of_interest:list, tissues_of_interest:list, threshold:float=0.05, id_type:str='gene_code', verbose:bool=False):
        if len(genes_of_interest) == 0:
            raise Exception("Genes of interest ('genes_of_interest:list') are mandatory; this list cannot be empty or none.")

        if id_type == 'gene_code' or id_type == 'gene_symbol':
            pass
        else:
            raise Exception("id_type value is not supported, you can only use 'gene_code' or 'gene_symbol'.")
    
        self.genes_of_interest = genes_of_interest
        self.id_type = id_type
        self.tissues_of_interest = tissues_of_interest
        self.threshold = threshold
        self.verbose = verbose

        if self.verbose == True:
            print( "***********************************************************" )
            print( " -----------" )
            print( " > CoRTE < " )
            print( " -----------" )
            print( "  Version: 0.1" )
            print( "  GitHub:\t" + "https://github.com/pietrocinaglia/corte" )
            print( "  Contact:\tPietro Cinaglia (cinaglia@unicz.it)" )
            print( "***********************************************************" )
        pass
  
    def __retrieve_data(self, action:str='geneExpression', genes_of_interest=None) -> pd.DataFrame:
        result = None

        if genes_of_interest is None:
            genes_of_interest = self.genes_of_interest

        if action == 'geneExpression':
            api = "https://gtexportal.org/api/v2/expression/geneExpression"
            params = {'itemsPerPage': self.LIMIT, 'datasetId': self.DATASETID, 'gencodeId': genes_of_interest, 'tissueSiteDetailId': self.tissues_of_interest, 'attributeSubset': 'ageBracket', 'format':'json'}
            result =  requests.get( api, params ).json()

        elif action == 'metasoft':
            api = "https://gtexportal.org/api/v2/association/metasoft"
            params = {'datasetId': self.DATASETID, 'gencodeId': genes_of_interest, 'format':'json'}
            result =  requests.get( api, params ).json()
            result = pd.DataFrame(result)

        elif action == 'metadata':
            api = "https://gtexportal.org/api/v2/reference/gene"
            params = {'itemsPerPage': self.LIMIT, 'geneId': genes_of_interest}
            result =  requests.get( api, params ).json()

        else:
            raise Exception("The action in Data Retrieving is not supported.")
        
        # Checking result
        if result is None:
            raise Exception("An error occurring in data retrieving.")

        return pd.DataFrame(result['data'])
    
    def construct_temporal_network(self) -> list:
        if self.verbose == True:
            print("[INFO] Data Retrieving...")
        
        gencodeIds = None
        if self.id_type == 'gene_symbol':
            genes_metadata = self.__retrieve_data(action='metadata', genes_of_interest=self.genes_of_interest)
            gencodeIds = genes_metadata['gencodeId'].values
        elif self.id_type == 'gene_code':
            gencodeIds = self.genes_of_interest # translation is not necessary
        else:
            raise Exception("id_type value is not supported, you can only use 'gene_code' or 'gene_symbol'.")
        
        gene_symbol_id = {self.genes_of_interest[i]: gencodeIds[i] for i in range(len(self.genes_of_interest))}

        if self.verbose == True:
            print("- Metadata [OK] ")
        
        df = self.__retrieve_data(action='geneExpression', genes_of_interest=gencodeIds)
        
        gene_pairs = list(combinations(gene_symbol_id, 2))

        if self.verbose == True:
            print("- Gene Expression Data [OK] ")
            print()

        if self.verbose == True:
            print("[INFO] Temporal Network Construction...")

        temporal_network = list()
        for i in range(len(self.AGES)):

            if self.verbose == True:
                print("- Time Point #" + str(i) + " (" + str(self.AGES[i]) + ") ...")

            timepoint = nx.Graph()
            timepoint.add_nodes_from( list(gene_symbol_id.keys()) )

            for u, v in gene_pairs:
                
                u_data = df.loc[ (df['gencodeId'] == gene_symbol_id[u]) & (df['tissueSiteDetailId'].isin(self.tissues_of_interest)) & (df['unit'] == self.UNIT) & (df['subsetGroup'] == self.AGES[i]) ]
                v_data = df.loc[ (df['gencodeId'] == gene_symbol_id[v]) & (df['tissueSiteDetailId'].isin(self.tissues_of_interest)) & (df['unit'] == self.UNIT) & (df['subsetGroup'] == self.AGES[i]) ]

                if (u_data.empty or v_data.empty):
                    continue

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
                    if self.verbose == True:
                        raise Warning("Size of sample less than 3; skipped.")
                    continue

                s, p = scipy.stats.pearsonr(u_gexp, v_gexp)
                
                p = round(p, 5)

                if p < self.threshold:
                    timepoint.add_edge(u,v, pvalue=p)

            temporal_network.append(timepoint)
            if self.verbose == True:
                print("-- [OK]")
        
        if self.verbose == True:
            print("- Temporal Network built.")

        return temporal_network

    def plot(self, temporal_network:list, with_labels:bool=True, output_path:str=None):
        if len(temporal_network) == 0:
            raise Exception("Temporal Network is empty; firtsly, you must call 'construct_temporal_network'.")

        # Print output
        for i in range(len(temporal_network)):
            tp = temporal_network[i]
            tp_title = "Timepoint #" + str(i)
            if self.verbose == True:
                print(tp_title)
                print(tp.nodes)
                print(tp.edges(data=True))
                print()
            
            plt.figure()  
            plt.title(tp_title)

            nx.draw(tp, with_labels=with_labels)

            if output_path is not None:
                plt.savefig(output_path + tp_title + '.png')
                plt.clf()
            else:
                plt.show()

    def save_as_files(self, temporal_network:list, output_path:str):
        if len(temporal_network) == 0:
            raise Exception("Temporal Network is empty; firtsly, you must call 'construct_temporal_network'.")
        
        if output_path is None:
            raise Exception("You must provide an output path for storing data.")

        # Store timepoints as file
        for i in range(len(temporal_network)):
            tp = temporal_network[i]
            tp_title = "timepoint" + str(i)
            nx.write_edgelist(tp, output_path+tp_title+'.txt', data=True)