import scanpy as sc
import warnings
import numpy as np
from sklearn.cluster import KMeans
from scipy.sparse import csr_matrix
import scgpt_spatial

warnings.filterwarnings('ignore')

adata = sc.read_visium(path="/content/scGPT-spatial/data/sample_1/filtered_feature_bc_matrix", count_file="matrix.h5")
adata.X = adata.X.toarray()

adata.var_names_make_unique()
sc.pp.filter_genes(adata, min_cells=3) 

coords = adata.obsm['spatial']
expr = adata.X.toarray() if not isinstance(adata.X, np.ndarray) else adata.X

model_dir = '/content/scGPT-spatial/scGPT_spatial_v1'
gene_col = 'index'

ref_embed_adata = scgpt_spatial.tasks.embed_data(
    adata,
    model_dir,
    gene_col=gene_col,
    obs_to_save=None,
    batch_size=64,
    return_new_adata=True,
    use_fast_transformer=False
)

# Cluster 
kmeans = KMeans(n_clusters=7).fit(ref_embed_adata.obsm["X_scGPT"])
ref_embed_adata.obs['scGPT_cluster'] = kmeans.labels_.astype(str)

# Visualize in spatial context
sc.pl.spatial(ref_embed_adata, color='scGPT_cluster', spot_size=1.2)