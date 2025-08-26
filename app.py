import streamlit as st
import os
import scanpy as sc
import pandas as pd
import plotly.express as px
import zipfile
import glob
import shutil
from pipeline import load_dataset, integrate_and_cluster, run_meta_analysis

# Configuration
st.set_page_config(layout="wide")
st.title("Single-Cell RNA-seq Meta-Analysis Tool")

# Create directory structure
os.makedirs("data", exist_ok=True)
os.makedirs("results", exist_ok=True)

# File format options
FORMATS = {
    "10X Genomics": "10x",
    "AnnData (.h5ad)": "h5ad",
    "Loom": "loom",
    "CSV/TSV Matrix": "csv"
}

with st.sidebar:
    st.header("Analysis Parameters")
    num_datasets = st.number_input("Number of datasets", min_value=1, max_value=20, value=2)
    resolution = st.slider("Clustering resolution", 0.1, 2.0, 1.0, 0.1)
    min_genes = st.number_input("Minimum genes per cell", min_value=200, value=500)
    min_cells = st.number_input("Minimum cells per gene", min_value=1, value=3)

# File upload section
uploaded_dirs = []
missing_samples = []

for i in range(num_datasets):
    with st.expander(f"Dataset {i+1}", expanded=True):
        col1, col2 = st.columns(2)
        with col1:
            file_type = st.selectbox(
                f"Format for dataset {i+1}",
                options=list(FORMATS.keys()),
                key=f"format_{i}",
                index=0  # Default to 10X
            )
        with col2:
            uploaded = st.file_uploader(
                f"Upload dataset {i+1}",
                type=['zip', 'h5ad', 'loom', 'csv', 'tsv'],
                key=f"file_{i}"
            )
        
        if uploaded:
            dataset_path = f"data/dataset{i+1}"
            os.makedirs(dataset_path, exist_ok=True)
            
            try:
                if uploaded.name.endswith('.zip'):
                    # Handle zip files
                    zip_path = os.path.join(dataset_path, "uploaded.zip")
                    with open(zip_path, "wb") as f:
                        f.write(uploaded.read())
                    
                    # Extract all files preserving directory structure
                    with zipfile.ZipFile(zip_path, 'r') as zip_ref:
                        zip_ref.extractall(dataset_path)
                    
                    # Verify we have the required files
                    has_matrix = len(glob.glob(os.path.join(dataset_path, '*matrix.mtx*'))) > 0
                    has_features = len(glob.glob(os.path.join(dataset_path, '*features.tsv*'))) > 0 or \
                                 len(glob.glob(os.path.join(dataset_path, '*genes.tsv*'))) > 0
                    has_barcodes = len(glob.glob(os.path.join(dataset_path, '*barcodes.tsv*'))) > 0
                    
                    if not all([has_matrix, has_features, has_barcodes]):
                        available = "\n".join(os.listdir(dataset_path))
                        st.error(f"Missing required files in {uploaded.name}\nAvailable files:\n{available}")
                        continue
                    
                    uploaded_dirs.append((dataset_path, FORMATS[file_type]))
                else:
                    # Handle direct file uploads
                    file_path = os.path.join("data", uploaded.name)
                    with open(file_path, "wb") as f:
                        f.write(uploaded.read())
                    uploaded_dirs.append((file_path, FORMATS[file_type]))
                    
            except Exception as e:
                missing_samples.append(f"Dataset {i+1}")
                st.error(f"⚠️ Failed to process dataset {i+1}: {str(e)}")
                continue

# Analysis section
if st.button("Run Analysis", type="primary"):
    if len(uploaded_dirs) < 2:
        st.error("❌ At least 2 valid datasets required for analysis")
        st.stop()
    
    if missing_samples:
        st.warning(f"⚠️ Proceeding with {len(uploaded_dirs)} datasets. Failed to load: {', '.join(missing_samples)}")

    try:
        # Load datasets with progress bar
        st.info("🔍 Preprocessing datasets...")
        adatas = []
        loaded_samples = []
        
        progress_bar = st.progress(0)
        for i, (path, ftype) in enumerate(uploaded_dirs):
            try:
                with st.spinner(f"Loading dataset {i+1}..."):
                    ad = load_dataset(path, f"dataset{i+1}", ftype)
                    
                    # Basic QC filtering
                    sc.pp.filter_cells(ad, min_genes=min_genes)
                    sc.pp.filter_genes(ad, min_cells=min_cells)
                    
                    adatas.append(ad)
                    loaded_samples.append(f"Dataset {i+1}")
                    progress_bar.progress((i + 1) / len(uploaded_dirs))
            except Exception as e:
                st.warning(f"⚠️ Failed to load dataset {i+1}: {str(e)}")
                continue
        
        if len(adatas) < 2:
            st.error("❌ Analysis requires at least 2 successfully loaded datasets")
            st.stop()
        
        st.success(f"✅ Successfully loaded {len(adatas)} datasets: {', '.join(loaded_samples)}")
        
        # Integration and clustering
        st.info("🧬 Integrating and clustering datasets...")
        with st.spinner("Running integration..."):
            adata = integrate_and_cluster(adatas)
            sc.write("results/final_integrated.h5ad", adata)
        st.success("✅ Integration & clustering complete!")
        
        # Visualization
        st.subheader("UMAP Visualizations")
        col1, col2 = st.columns(2)
        
        with col1:
            st.markdown("**By Batch**")
            fig1 = px.scatter(
                x=adata.obsm['X_umap'][:, 0],
                y=adata.obsm['X_umap'][:, 1],
                color=adata.obs['batch'],
                title="Batch Composition",
                labels={"x": "UMAP1", "y": "UMAP2"},
                width=500, height=500
            )
            st.plotly_chart(fig1, use_container_width=True)
        
        with col2:
            st.markdown("**By Cluster**")
            fig2 = px.scatter(
                x=adata.obsm['X_umap'][:, 0],
                y=adata.obsm['X_umap'][:, 1],
                color=adata.obs['leiden'],
                title="Cluster Composition",
                labels={"x": "UMAP1", "y": "UMAP2"},
                width=500, height=500
            )
            st.plotly_chart(fig2, use_container_width=True)
        
        # DEG analysis
        st.info("🔬 Running meta-analysis of DEGs...")
        with st.spinner("Identifying conserved markers..."):
            meta_degs = run_meta_analysis(adata)
        
        # Save and display results
        deg_file = "results/meta_degs.txt"
        with open(deg_file, "w") as f:
            for clus, genes in meta_degs.items():
                f.write(f"Cluster {clus} ({len(genes)} genes):\n")
                f.write(", ".join(genes) + "\n\n")
        
        # Show summary
        st.subheader("Meta-Analysis Results")
        deg_df = pd.DataFrame([
            {"Cluster": clus, "Conserved Genes": len(genes)} 
            for clus, genes in meta_degs.items()
        ])
        st.dataframe(deg_df)
        
        # Download button
        with open(deg_file, "rb") as f:
            st.download_button(
                "💾 Download Full DEG Results",
                f,
                file_name="meta_degs.txt",
                mime="text/plain"
            )
        
        st.success("🎉 Analysis complete!")
        st.balloons()
        
    except Exception as e:
        st.error(f"❌ Analysis failed: {str(e)}")
        st.stop()

# Clean up function
def cleanup():
    if os.path.exists("data"):
        shutil.rmtree("data", ignore_errors=True)
    if os.path.exists("results"):
        shutil.rmtree("results", ignore_errors=True)

# Register cleanup
import atexit
atexit.register(cleanup)
