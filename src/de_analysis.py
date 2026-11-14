import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from statsmodels.stats.multitest import multipletests

os.makedirs("results", exist_ok=True)

# Load normalized expression matrix (genes x samples)
df = pd.read_csv("data/GSE350063_expression_matrix_tsv.tsv", sep="\t", index_col=0)

wt = [c for c in df.columns if c.startswith("WT")]
ko = [c for c in df.columns if c.startswith("KO")]

# Filter: expression > 1 in at least 3 samples
df = df[(df > 1).sum(axis=1) >= 3]

# Log transform
log_df = np.log2(df + 1)

# Welch t-test per gene
t, p = stats.ttest_ind(log_df[ko], log_df[wt], axis=1, equal_var=False)
res = pd.DataFrame({
    "log2FC": log_df[ko].mean(axis=1) - log_df[wt].mean(axis=1),
    "pvalue": p,
}, index=log_df.index).dropna()

# FDR correction
res["padj"] = multipletests(res["pvalue"], method="fdr_bh")[1]
res = res.sort_values("padj")
res.to_csv("results/de_results.csv")

sig = (res["padj"] < 0.05) & (res["log2FC"].abs() > 1)

# Volcano plot
res["neg_log10_padj"] = -np.log10(res["padj"].clip(lower=1e-300))
plt.figure(figsize=(7, 6))
plt.scatter(res.loc[~sig, "log2FC"], res.loc[~sig, "neg_log10_padj"], s=8, c="gray", alpha=0.6)
plt.scatter(res.loc[sig, "log2FC"], res.loc[sig, "neg_log10_padj"], s=8, c="red")
plt.axhline(-np.log10(0.05), ls="--", c="k", lw=0.5)
plt.axvline(-1, ls="--", c="k", lw=0.5)
plt.axvline(1, ls="--", c="k", lw=0.5)
plt.xlabel("log2 fold change (KO vs WT)")
plt.ylabel("-log10 adj. p-value")
plt.title("Volcano plot")
plt.savefig("results/volcano.png", dpi=200, bbox_inches="tight")
plt.close()

# Heatmap of top 50 genes
top = res.index[:50]
g = sns.clustermap(log_df.loc[top], z_score=0, cmap="vlag", figsize=(7, 10),
                   yticklabels=True)
g.savefig("results/heatmap.png", dpi=200)

print(f"{len(res)} genes tested, {sig.sum()} significant")