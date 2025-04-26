
# Import necessary libraries
import streamlit as st
import pandas as pd
import anndata as ad
import scanpy as sc
import os
import subprocess
from chatbot import ask_chatbot
from utils import volcano_plot, pathway_analysis

# Add this early in your script
st.markdown("""
    <style>
        * {
            font-family: 'Arial', sans-serif !important;
        }
    </style>
""", unsafe_allow_html=True)

# Load the metadata file containing the information about the compounds used
metadata = pd.read_csv("/home/pranayagarwal/Documents/tahoe_100M/Data/moa_drug_final.csv")
unique_drugs = metadata['Drug_name'].unique().tolist()
unique_drugs = ['Select a drug'] + unique_drugs
cell_line = pd.read_csv("/home/pranayagarwal/Documents/tahoe_100M/Data/cell_line.csv")
unique_indication = cell_line['lineage'].unique().tolist()
unique_indication = ['All'] + unique_indication

st.title("Differential expression analysis on ARC's Single cell dataset...")

# User input: filter
query_1 = st.text_input("Enter a keyword to look for drugs:")
# Filter based on query (case-insensitive search across all string columns)
if query_1:
    mask = metadata.apply(lambda row: row.astype(str).str.contains(query_1, case=False).any(), axis=1)
    filtered_df =metadata[mask]
    st.dataframe(filtered_df[['Drug_name', 'Target', 'Mechanism of action ']])

# Find out if user is interested in any drug
drug = st.selectbox("Choose a drug:", unique_drugs)
if drug != 'Select a drug':
    st.write(ask_chatbot(f"Please find information on {drug}? Please find if it there are known mechanism of action? Please provide sources", "gpt3"))
drug_conc = st.selectbox("Choose the concentration (um):", ['5.0', '0.5', '0.05'])
tumor_indication = st.selectbox("Choose a tumor indication to analyze:", unique_indication)

if st.button("Run DESeq2"):
    try:
        result = subprocess.run(
            ["python", "deseq.py", drug, drug_conc, tumor_indication],
            capture_output=True,
            text=True,
            check=True
        )
        st.success("DESeq2 run completed!")
        # Load result (e.g. from CSV)

        output_path = f"/home/pranayagarwal/Documents/tahoe_100M/Results/{drug}_{drug_conc}_vs_DMSO_TF_in_{tumor_indication}.csv"
        if os.path.exists(output_path):
            res_df = pd.read_csv(output_path)
            res_df = pd.read_csv(f"/home/pranayagarwal/Documents/tahoe_100M/Results/{drug}_{drug_conc}_vs_DMSO_TF_in_{tumor_indication}.csv")
            st.download_button("Download DEG results", res_df.to_csv(index=False), file_name="DESeq2_results.csv")
            st.write('Top 50 differentially expressed genes:')
            st.dataframe(res_df.head(50))
            st.write('Volcano plot:')
            fig = volcano_plot(res_df, treatment = drug, ref = 'DMSO', top_n=3)
            st.pyplot(fig)
            
            # Running pathway analysis on upregulated genes
            up_reg_genes = res_df[(res_df['significant'] == True) & (res_df['log2FoldChange'] > 1)]
            up, fig = pathway_analysis(up_reg_genes.iloc[:,0].to_list())
            st.write('Pathway analysis with upregulated genes:')
            st.dataframe(up.results.head(10))
            # Running pathway analysis on downregulated genes
            down_reg_genes = res_df[(res_df['significant'] == True) & (res_df['log2FoldChange'] < -1)]
            down, fig = pathway_analysis(down_reg_genes.iloc[:,0].to_list())
            st.write('Pathway analysis with downregulated genes:')
            st.dataframe(down.results.head(10))
        else:
            st.error(f"Result file not found at {output_path}")
        
    except subprocess.CalledProcessError as e:
        st.error(f"Error: {e.stderr}")
        
        
        
        
        
        
        
        
        
        
        
        
