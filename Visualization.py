import scanpy as sc
import pandas as pd
import gseapy as gp
import matplotlib
import matplotlib.pyplot as plt

# 设置中文字体（适用于 Mac，自动回退）
matplotlib.rcParams['font.sans-serif'] = ['PingFang SC', 'Arial Unicode MS', 'SimHei', 'Heiti TC', 'STHeiti', 'Microsoft YaHei']
matplotlib.rcParams['axes.unicode_minus'] = False
adata = sc.read_h5ad("adata_all.h5ad")

# 或
adata_45 = sc.read_h5ad("adata_45.h5ad")
adata_58 = sc.read_h5ad("adata_58.h5ad")
adata_58_filtered = sc.read_h5ad("adata_58_filtered.h5ad")


# 可视化
sc.pl.umap(adata, color=['batch', 'celltype'], wspace=0.4)
# Day 45 聚类图
sc.pl.umap(adata_45, color='leiden', title='Day45 Clusters', legend_loc='on data', save='_day45_clusters1.png')
# Day 58 聚类图
sc.pl.umap(adata_58_filtered, color='leiden', title='Day58 Clusters', legend_loc='on data', save='_day58_clusters1.png')
# celltype CHN
sc.pl.umap(adata_45, color='celltype', title='Day45 Cell Types', legend_loc='on data', save='_day45_celltype.png')
sc.pl.umap(adata_58_filtered, color='celltype', title='Day58 Cell Types', legend_loc='on data', save='_day58_celltype.png')
