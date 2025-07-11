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

# 差异表达分析（在同一cluster中对比时间点）
# 例如对cluster 1比较 D45 vs D58
cluster_id = '1'
subset = adata[adata.obs['leiden'] == cluster_id]
sc.tl.rank_genes_groups(subset, groupby='batch', method='wilcoxon', reference='D45')
sc.pl.rank_genes_groups(subset, n_genes=25, sharey=False)

# 提取上调基因
df = sc.get.rank_genes_groups_df(subset, group='D58')
up_genes = df[(df['pvals_adj'] < 0.05) & (df['logfoldchanges'] > 0.5)]['names'].tolist()

# 通路富集分析（GSEApy）
enr = gp.enrichr(gene_list=up_genes,
                 gene_sets=['GO_Biological_Process_2021', 'KEGG_2021_Human', 'MSigDB_Hallmark_2020'],
                 organism='Human',
                 outdir=None)
gp.barplot(enr.res2d, title='Enriched Pathways (Day58 > Day45)', cutoff=0.05, top_term=10)