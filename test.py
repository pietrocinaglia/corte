# -*- coding: utf-8 -*-

import os

from tgeconet import TGECONET

###
# @author: Pietro Cinaglia
# @mail: cinaglia@unicz.it
# @description: tGeCoNet: a framework for constructing (t)emporal (Ge)ne (Co)-expression (Net)works
# @url: https://github.com/pietrocinaglia/tgeconet
###

WORKSPACE = os.path.dirname(os.path.realpath(__file__)) + "/"

genes_of_interest = ['EWSR1','SMARCA4','DDB2','YAP1','PSMD14','PEBP1','ITPKB','ATF7IP']
tissues_of_interest = ['Brain_Amygdala','Brain_Anterior_cingulate_cortex_BA24','Brain_Caudate_basal_ganglia','Brain_Cerebellar_Hemisphere','Brain_Cerebellum','Brain_Cortex','Brain_Frontal_Cortex_BA9','Brain_Hippocampus','Brain_Hypothalamus','Brain_Nucleus_accumbens_basal_ganglia','Brain_Putamen_basal_ganglia','Brain_Spinal_cord_cervical_c-1','Brain_Substantia_nigra']

print("################################################")
print("> TESTING <")
print("- Genes: " + str(genes_of_interest))
print("- Tissues: " + str(tissues_of_interest))
print("################################################")
print()

# Instantiate
tgeconet = TGECONET(
    genes_of_interest=genes_of_interest,
    tissues_of_interest=tissues_of_interest,
    threshold=0.05,
    verbose=True
)

# Build temporal network
temporal_network = tgeconet.construct_temporal_network()

# Statistics concerning the temporal network
stats = tgeconet.analyze_temporal_network(temporal_network, output_path=WORKSPACE+'results_from_test/')

# Extract top (n) genes
top_genes = tgeconet.extract_high_degree_genes(temporal_network, top_n=10)
print(top_genes)

# Export temporal network as adjacency_matrices with pvalues (layer by layer)
#tgeconet.export_adjacency_matrices(temporal_network, output_path=WORKSPACE+'results_from_test/')

# Export temporal network by using snapshot-based representation as set of edgelist
# tgeconet.save_as_files(temporal_network,output_path=WORKSPACE+'results_from_test/')

# Plotting temporal network and storing plots as image
# (if 'output_path' is defined, then each snapshot will be stored as image; it is not mandatory)
tgeconet.plot(temporal_network, with_labels=True) #, output_path=WORKSPACE)

print(">>> DONE <<<")
