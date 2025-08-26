
**Single-cell RNA-seq scRNAseq Meta-analysis Tool**
designed to:   

     Integrate heterogeneous datasets from multiple labs/platforms.  
     Identify conserved cell clusters and differentially expressed genes (DEGs) across batches.  
     Provide an interactive Streamlit dashboard for parameter tuning and visualization.

<img width="3755" height="1628" alt="Screenshot from 2025-08-26 11-22-11" src="https://github.com/user-attachments/assets/4e8a2826-52bc-4b76-a3f3-a7169e62aca0" />

     



Built with Scanpy, Harmony, and Streamlit, MetaSC bridges computational complexity and biological insight, enabling reproducible meta-analyses even for non-experts. 
🔍 Key Features 
🧪 Core Analysis Pipeline 

     Batch Correction: Uses Harmony integration to align datasets across technical batches 

.  
 Clustering: Leiden algorithm on UMAP-embedded data ensures robust cell-type identification 
.  
 Conserved DEGs: Identifies genes consistently dysregulated across batches, reducing technical noise 

<img width="3706" height="1944" alt="Screenshot from 2025-08-26 11-26-48" src="https://github.com/user-attachments/assets/6bd30688-8219-418c-9195-66bfec2b23b7" />

<img width="3583" height="1046" alt="Screenshot from 2025-08-26 11-27-02" src="https://github.com/user-attachments/assets/c15a73f8-2732-4732-aabe-7324d2f7e73a" />
    .
     

🌐 User-Friendly Interface 

     Streamlit Dashboard:  
         Upload datasets (supports 10X, AnnData, CSV formats) 

.  
 Visualize UMAP plots colored by batch or cluster 
.  
 Adjust parameters (e.g., clustering resolution, gene filters) in real-time 
.  
 Download DEG lists and integrated datasets for further analysis 

        .
         
     

🛠️ Flexibility & Robustness 

     Input Formats: Process CSV, 10X mtx, or preprocessed AnnData objects.  
     Error Handling: Built-in checks for missing files or incompatible data 

    .  
     Reproducibility: Save full analysis workflows with parameter logs.
     

📦 Quick Start 
     
     Installation 
     bash command-
     pip install metasc  # Coming soon to PyPI!

     git clone https://github.com/your-username/metascrepo.git
cd scRNA_folder/

      pip install -r requirements.txt
 
 
      Run the Dashboard 
       bash Command

     
      streamlit run app.py
 
 For downloading samples- User can use own single cell transcriptomic data to compare for metaanalysis
 or Take NCBI GEO dataset available scRNA dataset
 eg. In this study we used NCBİ GEO dataset to do scRNAs metaanalysis from case study https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE242780
   case study 1 sample 1- GSM7770408.zip 
   case study 2 sample 2- GSM7770409.zip
                          with files barcodes.tsv, genes.tsv, matrıx.mtx
These files are big so cann't be uploaded in github
                         
Open your browser to http://localhost:8501 and start analyzing! 
📊 Example Workflow 

     Upload Datasets: Drag-and-drop multiple 10X or AnnData files.  
     Configure Parameters:  
         Set min_genes=500 to filter low-quality cells 

.  
 Adjust resolution=0.8 for cluster granularity 

        .
         
     Analyze & Visualize:  
         View batch-corrected UMAP plots.  
         Download DEG lists for conserved cell populations.
         
     

📖 Documentation & Tutorials 

     User Guide: Documentation Link   
     Case Study: Analyze pancreatic cancer datasets (Jupyter notebook included).  
     API Reference: Explore functions in pipeline.py for advanced customization.
     

💡 Why Use MetaSC? 

     For Biologists: No coding required. Perform meta-analyses in hours, not weeks.  
     For Developers: Extendable pipeline with modular code.  
     For Researchers: Reproducible, publication-ready results with conserved DEGs.
     

🤝 Contributing 

     Report issues or suggest improvements via GitHub.  
     Contribute code or documentation to expand functionality.
     
     

🌐 Acknowledgments 

     Built with ❤️ using Scanpy, Harmony, and Streamlit.  
     Inspired by the open-science community.
     

License: MIT

Maintained by: Digianalix
This project is licensed under the MIT License.

For any issues or suggestions, please reach out to:

Author: Dr Mukesh nitin

Email: drmukeshnitin@gmail.com

https://github.com/drmukeshnitin/Single-Cell-RNA-seq-Meta-Analysis-Tool


     
