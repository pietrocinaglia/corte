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

#genes_symbols = ['PSMD14', 'PEBP1', 'ITPKB', 'ATF7IP']
genes_symbols = ['APBA1','APC','APH1B','BACE1','NOTCH1','PSEN2']# ['ENSG00000143801.16', 'ENSG00000134982.16', 'ENSG00000107282.7', 'ENSG00000148400.9', 'ENSG00000186318.16', 'ENSG00000138613.13']
tissues_of_interest = ['Brain_Amygdala','Brain_Anterior_cingulate_cortex_BA24','Brain_Caudate_basal_ganglia','Brain_Cerebellar_Hemisphere','Brain_Cerebellum','Brain_Cortex','Brain_Frontal_Cortex_BA9','Brain_Hippocampus','Brain_Hypothalamus','Brain_Nucleus_accumbens_basal_ganglia','Brain_Putamen_basal_ganglia','Brain_Spinal_cord_cervical_c-1','Brain_Substantia_nigra']

print("################################################")
print("> TESTING <")
print("- Genes: " + str(genes_symbols))
print("- Tissues: " + str(tissues_of_interest))
print("################################################")
print()

corte = CORTE(genes_symbols, tissues_of_interest, threshold=1, verbose=True)

# Data processing
temporal_network = corte.construct_temporal_network()
#corte.plot(temporal_network, with_labels=True)#, output_path=WORKSPACE)

print(">>> DONE <<<")
