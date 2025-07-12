import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from memory_profiler import memory_usage # pip install memory-profiler

#############

from corte import CORTE

WORKSPACE = os.path.dirname(os.path.realpath(__file__)) + "/"

# --- Fake data generator (simulate REST API invocation)--- 
def generate_fake_expression_data(num_genes, num_tissues, num_samples=10):
    import numpy as np
    import pandas as pd

    # Fake arguments
    gene_symbols = [f"gene_{i}" for i in range(num_genes)]
    tissues = [f"tissue_{i}" for i in range(num_tissues)]

    # Create a dataframe in accordance with real-world data provided by GTEx
    # -- columns: geneSymbol, tissueSiteDetailId, subsetGroup (age), data
    data_rows = []
    age_groups = ['20-29', '30-39', '40-49']

    for tissue in tissues:
        for gene in gene_symbols:
            for age in age_groups:
                # Fake expression data: list of random floats for each sample
                expr = list(np.random.rand(num_samples))
                data_rows.append({
                    'geneSymbol': gene,
                    'tissueSiteDetailId': tissue,
                    'subsetGroup': age,
                    'data': expr
                })

    df = pd.DataFrame(data_rows)
    return df, gene_symbols, tissues


# --- Benchmark ---
def benchmark_construct_temporal_network(gene_counts, num_tissues, threshold=0.05, max_genes=10000):
    results = []

    for num_genes in gene_counts:
        if num_genes > max_genes:
            print(f"[SKIP] Skipping {num_genes} genes (exceeds max_genes={max_genes})")
            continue

        genes = [f"GENE{i}" for i in range(num_genes)]
        fake_metadata = pd.DataFrame({
            "gene_symbol": genes,
            "ensembl_id": [f"ENSG{i:09d}" for i in range(num_genes)]
        })

        # Retrieving fake data
        df, gene_symbols, tissues = generate_fake_expression_data(num_genes, num_tissues)

        corte = CORTE(
            genes_of_interest=genes,
            tissues_of_interest=[],
            threshold=threshold,
            verbose=False,
            metadata=fake_metadata
        )

        corte.gtex_data = df
        corte.metadata = {g: g for g in gene_symbols}
        corte.tissues_of_interest = tissues

        # Run and measure time + memory
        print(f"[TEST] Running for {num_genes} genes...")
        start = time.time()
        mem_usage, temporal_network = memory_usage((corte.construct_temporal_network,), retval=True, max_iterations=1)
        end = time.time()

        mem_peak = max(mem_usage)
        total_edges = sum(g.number_of_edges() for g in temporal_network)

        results.append({
            "num_genes": num_genes,
            "time_sec": round(end - start, 2),
            "memory_mb": round(mem_peak, 2),
            "total_edges": total_edges
        })

    df = pd.DataFrame(results)

    return df

# --- Plotting ---
def plot_scalability(df):
    fig, ax1 = plt.subplots(figsize=(10, 6))

    ax1.set_xlabel("Number of Genes")
    ax1.set_ylabel("Execution Time (s)", color="blue")
    ax1.plot(df["num_genes"], df["time_sec"], color="blue", marker="o", label="Time")
    ax1.tick_params(axis="y", labelcolor="blue")

    ax2 = ax1.twinx()
    ax2.set_ylabel("Memory Usage (MB)", color="red")
    ax2.plot(df["num_genes"], df["memory_mb"], color="red", marker="x", label="Memory")
    ax2.tick_params(axis="y", labelcolor="red")

    plt.grid(True)
    plt.tight_layout()
    plt.show()

# --- Main ---
if __name__ == "__main__":
    gene_sizes = [50,100,150,200,250,300,350,400,450,500,550,600,650,700,750,800,850,900,950,1000]
    output_csv=WORKSPACE+"benchmark_results.csv"
    
    results_df = benchmark_construct_temporal_network(
        gene_counts=gene_sizes,
        num_tissues=1, # the unit has been defined equal to one tissue
    )

    print("[INFO] Saving the results as file...")
    results_df.to_csv(output_csv, index=False)
    print(f"-- data saved to {output_csv}")

    print("[INFO] Plotting...")
    plot_scalability(results_df)

    print("[ DONE ]")
