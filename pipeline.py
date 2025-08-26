import scanpy as sc
import pandas as pd
import os
import anndata
import loompy
import glob
import warnings

def load_dataset(path, name, file_type=None):
    """Load single-cell data in multi formats with error handling"""
    try:
        if file_type is None:
            # Auto-detect format
            if os.path.isdir(path) and validate_10x_directory(path):
                return load_10x_data(path, name)
            elif path.endswith('.h5ad'):
                return sc.read_h5ad(path)
            elif path.endswith('.loom'):
                return load_loom_data(path, name)
            elif path.endswith(('.csv', '.tsv')):
                return load_csv_data(path, name)
            else:
                raise ValueError(f"Unsupported file format: {path}")
        else:
            if file_type == '10x':
                if not os.path.isdir(path):
                    raise ValueError(f"10X format requires a directory, got: {path}")
                return load_10x_data(path, name)
            elif file_type == 'h5ad':
                if os.path.isdir(path):
                    raise ValueError(f"h5ad format requires a file, got directory: {path}")
                return sc.read_h5ad(path)
            elif file_type == 'loom':
                if os.path.isdir(path):
                    raise ValueError(f"loom format requires a file, got directory: {path}")
                return load_loom_data(path, name)
            elif file_type == 'csv':
                if os.path.isdir(path):
                    raise ValueError(f"CSV format requires a file, got directory: {path}")
                return load_csv_data(path, name)
            else:
                raise ValueError(f"Unsupported file type: {file_type}")
    except Exception as e:
        raise ValueError(f"Failed to load dataset {name} from {path}: {str(e)}")

def load_10x_data(path, name):
    """Load 10X Genomics data with improved error handling"""
    try:
        convert_matrix_extensions(path)
        
        # Check for matrix file variants
        matrix_files = [
            os.path.join(path, "matrix.mtx.gz"),
            os.path.join(path, "matrix.mtx"),
            os.path.join(path, "matrix.mtx.tgz"),
            os.path.join(path, "matrix.mtx.zip")
        ]
        
        matrix_path = None
        for mtx_file in matrix_files:
            if os.path.exists(mtx_file):
                matrix_path = mtx_file
                break
        
        if matrix_path is None:
            available_files = "\n".join(os.listdir(path))
            raise FileNotFoundError(f"No matrix file found in {path}\nAvailable files:\n{available_files}")
        
        # Check for features/genes file
        features_files = [
            os.path.join(path, "features.tsv"),
            os.path.join(path, "genes.tsv"),
            os.path.join(path, "features.tsv.gz"),
            os.path.join(path, "genes.tsv.gz")
        ]
        
        features_path = None
        for f_file in features_files:
            if os.path.exists(f_file):
                features_path = f_file
                break
        
        if features_path is None:
            raise FileNotFoundError(f"No features/genes file found in {path}")
        
        # Load the data
        adata = sc.read_10x_mtx(
            path,
            var_names='gene_symbols',
            cache=True
        )
        adata.var_names_make_unique()
        return process_common_metadata(adata, name)
        
    except Exception as e:
        raise ValueError(f"Failed to load 10X data from {path}: {str(e)}")

def validate_10x_directory(path):
    """Check if directory contains valid 10X files"""
    required = {
        'matrix': ['matrix.mtx', '.matrix.mtx', '*.mtx'],
        'features': ['features.tsv', 'genes.tsv'],
        'barcodes': ['barcodes.tsv']
    }
    
    for file_type, patterns in required.items():
        found = False
        for pattern in patterns:
            if glob.glob(os.path.join(path, pattern)):
                found = True
                break
        if not found:
            return False
    return True

def convert_matrix_extensions(directory):
    """Convert any .matrix.mtx files to standard matrix.mtx"""
    for old_name in glob.glob(os.path.join(directory, '*.matrix.mtx')):
        new_name = os.path.join(directory, 'matrix.mtx')
        if not os.path.exists(new_name):
            os.rename(old_name, new_name)

def load_loom_data(path, name):
    """Load loom format data"""
    adata = sc.read_loom(path)
    return process_common_metadata(adata, name)

def load_csv_data(path, name):
    """Load from CSV/TSV matrix"""
    if path.endswith('.csv'):
        df = pd.read_csv(path, index_col=0)
    else:
        df = pd.read_csv(path, sep='\t', index_col=0)
    adata = anndata.AnnData(X=df.values, var=pd.DataFrame(index=df.columns))
    return process_common_metadata(adata, name)

def process_common_metadata(adata, name):
    """Common processing for all formats"""
    adata.obs['batch'] = name
    if adata.obs_names.empty:
        adata.obs_names = [f"{name}_cell_{i}" for i in range(adata.shape[0])]
    return adata

def integrate_and_cluster(adatas):
    """Concatenate and integrate datasets with error handling"""
    if len(adatas) < 2:
        raise ValueError("Need at least 2 datasets for integration")
    
    try:
        adata = sc.concat(adatas, label='batch', keys=[a.obs['batch'].unique()[0] for a in adatas])
        
        # Select highly variable genes across batches
        sc.pp.highly_variable_genes(adata, flavor='seurat_v3', n_top_genes=2000, batch_key='batch')
        adata = adata[:, adata.var['highly_variable']]

        # Scale and PCA
        sc.pp.scale(adata, max_value=10)
        sc.tl.pca(adata, svd_solver='arpack')

        # Harmony integration
        import harmonypy as hm
        ho = hm.run_harmony(adata.obsm['X_pca'], adata.obs, 'batch')
        adata.obsm['X_pca_harmony'] = ho.Z_corr.T

        # Clustering and UMAP
        sc.pp.neighbors(adata, use_rep='X_pca_harmony')
        sc.tl.umap(adata)
        sc.tl.leiden(adata, resolution=1.0)

        return adata
    except Exception as e:
        raise ValueError(f"Integration failed: {str(e)}")

def run_meta_analysis(adata):
    """Perform meta-analysis with error handling"""
    clusters = sorted(adata.obs['leiden'].unique())
    meta_degs = {}

    for cluster in clusters:
        try:
            subset = adata[adata.obs['leiden'] == cluster].copy()
            
            if subset.obs['batch'].nunique() < 2:
                continue

            sc.tl.rank_genes_groups(
                subset, 
                groupby='batch', 
                method='wilcoxon',
                use_raw=False
            )
            
            result = subset.uns['rank_genes_groups']
            df = pd.DataFrame({
                group: result['names'][group]
                for group in result['names'].dtype.names
            })

            sets = [set(df[col].dropna()) for col in df.columns]
            common = set.intersection(*sets) if sets else set()
            meta_degs[cluster] = list(common)
        except Exception as e:
            warnings.warn(f"Skipped cluster {cluster} due to error: {str(e)}")
            continue

    return meta_degs
