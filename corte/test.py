# -*- coding: utf-8 -*-

import os
from CORTE import CORTE

###
# @author: Pietro Cinaglia
# @mail: cinaglia@unicz.it
# @description: An open source and user-friendly tool for COnstructing Real-world TEmporal networks from genotype-tissue expression data (CoRTE)
# @url: https://github.com/pietrocinaglia/corte
###

WORKSPACE = os.path.dirname(os.path.realpath(__file__)) + "/"

genes_of_interest = ['PSMD14', 'PEBP1', 'ITPKB', 'ATF7IP'] # + ['CCDC92','CCDC74B','ARL14EPL','CCDC74A','ARL14EP','ZNF883','MINPP1','ITPKC','ITPKA','ITPKB','IPMK']
tissues_of_interest = ['Brain_Amygdala','Brain_Anterior_cingulate_cortex_BA24','Brain_Caudate_basal_ganglia','Brain_Cerebellar_Hemisphere','Brain_Cerebellum','Brain_Cortex','Brain_Frontal_Cortex_BA9','Brain_Hippocampus','Brain_Hypothalamus','Brain_Nucleus_accumbens_basal_ganglia','Brain_Putamen_basal_ganglia','Brain_Spinal_cord_cervical_c-1','Brain_Substantia_nigra']

print("################################################")
print("> TESTING <")
print("- Genes: " + str(genes_of_interest))
print("- Tissues: " + str(tissues_of_interest))
print("################################################")
print()

corte = CORTE(genes_of_interest, tissues_of_interest, verbose=True)

# Data processing
temporal_network = corte.construct_temporal_network()

# Plotting temporal network
corte.plot(temporal_network, with_labels=True)
# (or) Plotting temporal network and storing plots as image
# corte.plot(temporal_network, with_labels=True, output_path=WORKSPACE))

# Storing temporal network by using snapshot-based representation as set of edgelist
corte.save_as_files(temporal_network,output_path=WORKSPACE)

print(">>> DONE <<<")
