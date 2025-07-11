import scanpy as sc
import pandas as pd
import gseapy as gp
import matplotlib
import matplotlib.pyplot as plt

adata = sc.read_h5ad("adata_all.h5ad")
# 或
adata_45 = sc.read_h5ad("adata_45.h5ad")
adata_58 = sc.read_h5ad("adata_58.h5ad")
adata_58_filtered = sc.read_h5ad("adata_58_filtered.h5ad")


# 对 Day45 各 cluster 做 marker gene 排名
sc.tl.rank_genes_groups(adata_45, groupby='leiden', method='wilcoxon')
sc.pl.rank_genes_groups(adata_45, n_genes=10, sharey=False)
# 对 Day58 同样操作
sc.tl.rank_genes_groups(adata_58_filtered, groupby='leiden', method='wilcoxon')
sc.pl.rank_genes_groups(adata_58_filtered, n_genes=10, sharey=False)

# Day 45 导出所有 cluster 的差异基因
df_45 = sc.get.rank_genes_groups_df(adata_45, group=None)
df_45.to_csv("markers_day45.csv", index=False)
# Day 58 同理
df_58 = sc.get.rank_genes_groups_df(adata_58_filtered, group=None)
df_58.to_csv("markers_day58.csv", index=False)