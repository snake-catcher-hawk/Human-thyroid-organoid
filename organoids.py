import scanpy as sc
import pandas as pd
import gseapy as gp
import matplotlib
import matplotlib.pyplot as plt

# 设置路径
path_45 = "GSE181256_RAW/"
path_58 = "GSE203084_RAW/"

# 读取数据（10X格式）
adata_45 = sc.read_10x_mtx(path_45, var_names='gene_symbols', cache=True)
adata_58 = sc.read_10x_mtx(path_58, var_names='gene_symbols', cache=True)

# 添加时间标签
adata_45.obs['time'] = 'Day45'
adata_58.obs['time'] = 'Day58'

# 合并
adata = adata_45.concatenate(adata_58, batch_key='batch', batch_categories=['D45', 'D58'])


# 质控
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
adata = adata[adata.obs.n_genes_by_counts < 6000, :]
adata = adata[adata.obs.pct_counts_mt < 10, :]

# 归一化 + 高变基因 + PCA
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata = adata[:, adata.var.highly_variable]
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')

# 邻接图 + 聚类 + UMAP
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=30)
sc.tl.leiden(adata, resolution=0.5)
sc.tl.umap(adata)

# 英文字典
cluster_celltype = {
    "0": "Proliferative",
    "1": "Thyroid Progenitor",
    "2": "Thyroid Follicular",
    "3": "Fibroblast",
    "4": "Epithelial",
    "5": "Metabolic-like",
    "6": "Secretory-like",
    "7": "Mucosal progenitor",
    "8": "Cycling",
    "9": "Glandular epithelial",
    "10": "Neuroendocrine-like",
    "11": "Stem-like",
    "12": "Mature Thyroid Follicular",
    "13": "Mesenchymal"
}
# 中文字典
cluster_celltype_CHN = {
    "0": "增殖活跃细胞（Cycling / Ribosomal）",
    "1": "上皮祖细胞 / 甲状腺前体细胞（Thyroid Progenitor）",
    "2": "甲状腺滤泡细胞（Thyroid follicular）",  # 按你的要求已修正
    "3": "成纤维细胞（Fibroblast）",
    "4": "上皮细胞（Epithelial）",
    "5": "代谢样细胞（Metabolic-like）",
    "6": "分泌样细胞（Secretory-like）",
    "7": "粘膜上皮祖细胞（Mucosal progenitor）",
    "8": "强增殖细胞群（Cycling）",
    "9": "腺上皮细胞（Glandular epithelial）",
    "10": "神经内分泌样细胞（Neuroendocrine-like）",
    "11": "未分化/祖细胞样（Progenitor-like）",
    "12": "甲状腺成熟滤泡细胞（Thyroid follicular）",
    "13": "间充质细胞（Mesenchymal-like）"
}
adata.obs["celltype"] = adata.obs["leiden"].map(cluster_celltype_CHN)

# 提取子集
adata_45 = adata[adata.obs['time'] == 'Day45', :].copy()
adata_58 = adata[adata.obs['time'] == 'Day58', :].copy()

# 过滤掉只含1个细胞的cluster
valid_clusters = adata_58.obs['leiden'].value_counts()
valid_clusters = valid_clusters[valid_clusters > 1].index
adata_58_filtered = adata_58[adata_58.obs['leiden'].isin(valid_clusters), :].copy()
sc.tl.rank_genes_groups(adata_58_filtered, groupby='leiden', method='wilcoxon')

# 保存合并后的全部数据
adata.write("adata_all.h5ad")
# 保存 Day45 子集
adata_45.write("adata_45.h5ad")
# 保存 Day58 子集（可选：过滤后）
adata_58.write("adata_58.h5ad")
adata_58_filtered.write("adata_58_filtered.h5ad")
