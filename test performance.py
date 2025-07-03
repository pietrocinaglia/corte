# -*- coding: utf-8 -*-

import os
import time
import random
import pandas as pd
import matplotlib.pyplot as plt
import warnings
warnings.simplefilter('ignore')
###
# @author: Pietro Cinaglia
# @mail: cinaglia@unicz.it
# @description: An open source and user-friendly tool for COnstructing Real-world TEmporal networks from genotype-tissue expression data (CoRTE)
# @url: https://github.com/pietrocinaglia/corte
###

from corte import CORTE

WORKSPACE = os.path.dirname(os.path.realpath(__file__)) + "/"

# Define test parameters
gene_sizes = [50, 100, 250, 350, 450, 500]  # Increase gene list size
tissue_sizes = [10, 25, 50] # (54) tissues

# Load gene symbols
metadata_df = pd.read_csv(WORKSPACE+'data/metadata.csv', sep=" ")
ALL_GENES = metadata_df['gene_symbol'].unique().tolist()
ALL_GENES = [g for g in ALL_GENES if not g.startswith("ENS")]

# Load tissues
tissues_df = pd.read_csv(WORKSPACE+'data/tissues.csv', sep=",")
ALL_TISSUES = tissues_df['tissue_key'].unique().tolist()

# Record performance data
results = []
print("[INFO] Measuring performance...")

for g_size in gene_sizes:
    for t_size in tissue_sizes:
        genes = ALL_GENES
        tissues = ALL_TISSUES
        if g_size < len(ALL_GENES):
            genes = random.sample(ALL_GENES, g_size)
        if t_size < len(ALL_TISSUES):
            tissues = random.sample(ALL_TISSUES, t_size)

        print(f" - Genes: {g_size}, Tissues: {t_size}...", end=" ", flush=True)

        runtime = time.time()
        corte = CORTE(genes_of_interest=genes, tissues_of_interest=tissues, threshold=0.05, verbose=False)
        network = corte.construct_temporal_network()
        runtime = round(time.time() - runtime, 4)
        results.append({'genes': g_size, 'tissues': t_size, 'runtime': runtime})
        print(f" -- Runtime: {runtime}", end=" ", flush=True)

print("[INFO] Saving results...")

df_results = pd.DataFrame(results)
df_results.to_csv(WORKSPACE+"performance_results.csv", index=False)

print("[INFO] Plotting Heatmap...")

pivot = df_results.pivot("genes", "tissues", "time")
plt.figure(figsize=(8, 6))
plt.title("Runtime (s)")
heatmap = plt.imshow(pivot, cmap="viridis", aspect="auto", interpolation="nearest")
plt.colorbar(heatmap, label="seconds")
plt.xticks(ticks=range(len(pivot.columns)), labels=pivot.columns)
plt.yticks(ticks=range(len(pivot.index)), labels=pivot.index)
plt.xlabel("Number of Tissues")
plt.ylabel("Number of Genes")
plt.tight_layout()
plt.savefig(WORKSPACE+"performance_heatmap.png", dpi=300)
plt.show()

print("[ DONE ]")
