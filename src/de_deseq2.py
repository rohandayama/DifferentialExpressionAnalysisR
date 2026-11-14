import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

os.makedirs("results", exist_ok=True)

# Load counts (genes x samples), merge normal + tumor
normal = pd.read_csv("data/GSE353851_RNA_raw_readcounts_normal.tsv", sep="\t")
tumor = pd.read_csv("data/GSE353851_RNA_raw_readcounts_tumor.tsv", sep="\t")
symbols = normal.set_index("gene_id")["gene_symbol"]

normal = normal.drop(columns=["gene_symbol", "length"]).set_index("gene_id")
tumor = tumor.drop(columns=["gene_symbol", "length"]).set_index("gene_id")
counts = normal.join(tumor, how="inner")

# Sample labels
ann = pd.read_csv("data/GSE353851_sample_annotation.tsv", sep="\t").set_index("sample_id")
ann = ann[ann["rna_available"].astype(str).str.lower() == "true"]
shared = [s for s in counts.columns if s in ann.index]
counts = counts[shared]
meta = pd.DataFrame({"condition": ann.loc[shared, "sample_type"].values}, index=shared)
print(meta["condition"].value_counts())

# PyDESeq2 wants samples x genes
counts = counts.T.astype(int)

# Filter low-count genes: >= 10 counts in at least 10 samples
counts = counts.loc[:, (counts >= 10).sum(axis=0) >= 10]
print(f"{counts.shape[1]} genes after filtering")

# DESeq2
dds = DeseqDataSet(counts=counts, metadata=meta, design="~condition", n_cpus=4)
dds.deseq2()

stat = DeseqStats(dds, contrast=("condition", "tumor", "normal"), n_cpus=4)
stat.summary()
res = stat.results_df.dropna(subset=["padj"]).copy()
res["gene_symbol"] = symbols.reindex(res.index).values
res = res.sort_values("padj")
res.to_csv("results/deseq2_results.csv")

sig = (res["padj"] < 0.05) & (res["log2FoldChange"].abs() > 1)

# Volcano plot
res["neg_log10_padj"] = -np.log10(res["padj"].clip(lower=1e-300))
plt.figure(figsize=(7, 6))
plt.scatter(res.loc[~sig, "log2FoldChange"], res.loc[~sig, "neg_log10_padj"], s=6, c="gray", alpha=0.5)
plt.scatter(res.loc[sig, "log2FoldChange"], res.loc[sig, "neg_log10_padj"], s=6, c="red", alpha=0.6)
plt.axhline(-np.log10(0.05), ls="--", c="k", lw=0.5)
plt.axvline(-1, ls="--", c="k", lw=0.5)
plt.axvline(1, ls="--", c="k", lw=0.5)
plt.xlabel("log2 fold change (tumor vs normal)")
plt.ylabel("-log10 adj. p-value")
plt.title("Volcano plot")
plt.savefig("results/deseq2_volcano.png", dpi=200, bbox_inches="tight")
plt.close()

# Heatmap of top 50 genes (by padj), on log2 normalized counts
top = res.index[:50]
normed = pd.DataFrame(dds.layers["normed_counts"], index=dds.obs_names, columns=dds.var_names)
heat = np.log2(normed[top] + 1).T
heat.index = [str(symbols.get(g, g)) for g in heat.index]
col_colors = meta["condition"].map({"tumor": "firebrick", "normal": "steelblue"})
g = sns.clustermap(heat, z_score=0, cmap="vlag", figsize=(10, 10),
                   col_colors=col_colors, xticklabels=False, yticklabels=True)
g.savefig("results/deseq2_heatmap.png", dpi=200)

print(f"{len(res)} genes tested, {sig.sum()} significant")
