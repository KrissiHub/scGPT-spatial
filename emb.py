import scanpy as sc
import warnings
import numpy as np
from sklearn.cluster import KMeans
from scipy.sparse import csr_matrix
import scgpt_spatial

warnings.filterwarnings('ignore')

#Multiple samples must be concatenated to take all the data into account
base_path = "/content/scGPT-spatial/data"
sample_dirs = [f"{base_path}/sample_{i}/filtered_feature_bc_matrix" for i in range(1,7)]

adata = []

for i, p in enumerate(sample_dirs, start=1):
    ad = sc.read_visium(path=p, count_file="matrix.h5")
    ad.var_names_make_unique()

    ad.obs["sample"] = f"sample_{i}"

    adata.append(ad)

#Concat the samples
adata = sc.concat(
    adata,
    join="outer",
    label="sample",
    keys=[f"sample_{i}" for i in range(1, 7)],
    index_unique="-",          
    merge="same",              
    fill_value=0
)

adata.X = adata.X.toarray()
sc.pp.filter_genes(adata, min_cells=3) 

coords = adata.obsm['spatial']

model_dir = '/content/scGPT-spatial/scGPT_spatial_v1'
gene_col = 'index'

ref_embed_adata = scgpt_spatial.tasks.embed_data(
    adata,
    model_dir,
    gene_col=gene_col,
    obs_to_save=None,
    batch_size=64,
    return_new_adata=False,
    use_fast_transformer=False
)

# Cluster 
kmeans = KMeans(n_clusters=7).fit(ref_embed_adata.obsm["X_scGPT"])
ref_embed_adata.obs['scGPT_cluster'] = kmeans.labels_.astype(str)

# Visualize in spatial context

sc.settings.figdir = "/content/scGPT-spatial"
sc.pl.spatial(ref_embed_adata, color='scGPT_cluster', spot_size=50, save="/figure_2c.png")
