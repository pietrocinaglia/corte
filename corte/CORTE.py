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
    DATASETID = "gtex_v8"
    #
    temporal_network = list()
    genes_of_interest = list()
    tissues_of_interest = list()
    threshold = 0.05
    verbose = False
    #

    def __init__(self, genes_of_interest:list, tissues_of_interest:list, threshold:float=0.05, verbose:bool=False):
        if len(genes_of_interest) == 0:
            raise Exception("Genes of interest ('genes_of_interest:list') are mandatory; this list cannot be empty or none.")

        self.genes_of_interest = genes_of_interest
        self.tissues_of_interest = tissues_of_interest
        self.threshold = threshold
        self.verbose = verbose

        if self.verbose == True:
            print( "***********************************************************" )
            print( " -----------" )
            print( " > CoRTE < " )
            print( " -----------" )
            print( "  GitHub:\t" + "https://github.com/pietrocinaglia/corte" )
            print( "  Contact:\tPietro Cinaglia (cinaglia@unicz.it)" )
            print( "***********************************************************" )
        pass
  
    def __retrieve_data(self, action:str='geneExpression') -> pd.DataFrame:
        result = None
        if action == 'geneExpression':
            api = "https://gtexportal.org/api/v2/expression/geneExpression"
            params = {'itemsPerPage': 100000, 'datasetId': self.DATASETID, 'gencodeId': self.genes_of_interest, 'tissueSiteDetailId': self.tissues_of_interest, 'attributeSubset': 'ageBracket', 'format':'json'}
            result =  requests.get( api, params ).json()
        elif action == 'metasoft':
            api = "https://gtexportal.org/api/v2/association/metasoft"
            params = {'datasetId': self.DATASETID, 'gencodeId': self.genes_of_interest, 'format':'json'}
            result =  requests.get( api, params ).json()
            result = pd.DataFrame(result)
        elif action == 'metadata':
            api = "https://gtexportal.org/api/v2/reference/gene"
            params = {'itemsPerPage': 100000, 'geneId': self.genes_of_interest}
            result =  requests.get( api, params ).json()
        else:
            raise Exception("The action in Data Retrieving is not supported.")
        
        if result is not None:
            result = pd.DataFrame(result["data"])

        return result

    def __savefig(self, output_path:str = None):
        # Print output
        for i in range(len(self.temporal_network)):
            tp = self.temporal_network[i]
            if self.verbose == True:
                print("Timepoint #" + str(i))
                print(tp.nodes)
                print(tp.edges(data=True))
                print()
            plt.savefig(output_path + 'tp' + str(i) + '.png')
            plt.clf()
    
    def construct_temporal_network(self) -> list:
        if self.verbose == True:
            print("[INFO] Data Retrieving...")
        
        genes_metadata = self.__retrieve_data(self.genes_of_interest, tissues_of_interest=self.tissues_of_interest, action='metadata')
        gencodeIds = genes_metadata['gencodeId'].values
        gene_symbol_id = {self.genes_of_interest[i]: gencodeIds[i] for i in range(len(self.genes_of_interest))}
        if self.verbose == True:
            print("- Metadata [OK] ")

        df = self.__retrieve_data(gencodeIds, action='geneExpression')
        gene_pairs = list(combinations(gene_symbol_id, 2))
        if self.verbose == True:
            print("- Gene Expression Data [OK] ")

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
                    timepoint.add_edge(u,v, p=p)

            temporal_network.append(timepoint)
            if self.verbose == True:
                print("-- [OK]")
        
        if self.verbose == True:
            print("- Temporal Network built.")

        return temporal_network
    
    def plot(self, with_labels:bool=True, savefig:str=None):
        if len(self.temporal_network) == 0:
            raise Exception("Temporal Network was not constructed; firtsly, you must call 'construct_temporal_network'.")

        # Print output
        for i in range(len(self.temporal_network)):
            tp = self.temporal_network[i]
            if self.verbose == True:
                print("Timepoint #" + str(i))
                print(tp.nodes)
                print(tp.edges(data=True))
                print()
                        
            if savefig is not None:
                self.__savefig(self, savefig)
            
            nx.draw(tp, with_labels)
            plt.show()

    #############
    #############
    #############
    
    def test(self, plot:bool=True, savefig:str=None):
        self.genes_symbols = ['APBA1','APC','APH1B','BACE1','NOTCH1','PSEN2']#['ENSG00000141510.16', 'ENSG00000146648.17']
        self.tissues_of_interest = ['Brain_Amygdala','Brain_Anterior_cingulate_cortex_BA24','Brain_Caudate_basal_ganglia','Brain_Cerebellar_Hemisphere','Brain_Cerebellum','Brain_Cortex','Brain_Frontal_Cortex_BA9','Brain_Hippocampus','Brain_Hypothalamus','Brain_Nucleus_accumbens_basal_ganglia','Brain_Putamen_basal_ganglia','Brain_Spinal_cord_cervical_c-1','Brain_Substantia_nigra']
        print("################################################")
        print("> TESTING <")
        print("- Genes: " + str(self.genes_symbols))
        print("- Tissues: " + str(self.tissues_of_interest))
        print("################################################")
        print()
        # Data processing
        temporal_network = self.construct_temporal_network()
        # Print output
        for i in range(len(temporal_network)):
            tp = temporal_network[i]
            print("Timepoint #" + str(i))
            print(tp.nodes)
            print(tp.edges(data=True))
            print()

            if savefig is not None:
                self.__savefig(self, savefig+'tp'+str(i)+'.png')
            
            if plot == True:
                nx.draw(tp, with_labels=True)
                plt.show()            

        print(">>> DONE <<<")