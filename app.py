
# Import necessary libraries
import streamlit as st
import pandas as pd
import anndata as ad
import scanpy as sc
import os
import subprocess
from chatbot import ask_chatbot
from utils import volcano_plot, pathway_analysis, extract_data, diff_exp
import requests
import re
import io

# Add this early in your script
st.markdown("""
    <style>
        * {
            font-family: 'Arial', sans-serif !important;
        }
    </style>
""", unsafe_allow_html=True)

# Load the metadata file containing the information about the compounds used
drug_metadata = pd.read_csv("/home/pranayagarwal/Documents/tahoe_100M/Data/moa_drug_final.csv")
unique_drugs = drug_metadata['Drug_name'].unique().tolist()
unique_drugs = ['Select a drug'] + unique_drugs

cell_line = pd.read_csv("/home/pranayagarwal/Documents/tahoe_100M/Data/cell_line.csv")
unique_indication = cell_line['lineage'].unique().tolist()
unique_indication = ['All'] + unique_indication

st.title("Differential expression analysis on ARC's Single cell dataset...")

# User input: filter
query_1 = st.text_input("Enter a keyword to look for drugs:")
# Filter based on query (case-insensitive search across all string columns)
columns_to_search = ['Drug_name', 'Target', 'Mechanism of action']
if query_1:
    mask = drug_metadata[columns_to_search].apply(
        lambda col: col.astype(str).str.contains(query_1, case=False, na=False)
    ).any(axis=1)
    filtered_df = drug_metadata[mask]
    st.dataframe(filtered_df[['Drug_name', 'Target', 'Mechanism of action']])

# Find out if user is interested in any drug
drug = st.selectbox("Choose a drug:", unique_drugs)
if drug != 'Select a drug':
    # Get LLM output value for the drug
    llm_output = drug_metadata.loc[drug_metadata['Drug_name'] == drug, 'LLM_output'].values

    if len(llm_output) == 0 or pd.isna(llm_output[0]) or llm_output[0] == "":
        # Ask GPT if not already present
        answer = ask_chatbot(
            f"Please find information on {drug}? Please find if there are known mechanism of action? Please provide sources", 
            "gpt3"
        )
        st.write(answer)

        # Update LLM_output column
        drug_metadata.loc[drug_metadata['Drug_name'] == drug, 'LLM_output'] = answer
        drug_metadata.to_csv("/home/pranayagarwal/Documents/tahoe_100M/Data/moa_drug_final.csv", index=False)
    else:
        st.write(llm_output[0])
        
drug_conc = st.selectbox("Choose the concentration (um):", ['5.0', '0.5', '0.05'])
tumor_indication = st.selectbox("Choose a tumor indication to analyze:", unique_indication)

# Based on user's input find out which h5ad file to download from google cloud

metadata = pd.read_csv("/home/pranayagarwal/Documents/tahoe_100M/Data/meta_data.csv")
#st.dataframe(metadata.head())
metadata["concentration, uM"] = metadata["concentration, uM"].astype(str).astype('category')
filtered = metadata[
    (metadata['drug_name'] == drug) & (metadata["concentration, uM"] == drug_conc)
]
if not filtered.empty:
    plate_id = filtered['batch'].iloc[0]
else:
     st.warning("No matching plate ID found for the given drug and concentration.")

@st.cache_data
def load_h5ad_from_url(url):
    response = requests.get(url)
    adata = sc.read_h5ad(io.BytesIO(response.content))
    return adata

if st.button("Extract data and Run DESeq2"):
    # fetching data from google cloud
    plateid = re.sub(r"_", "", plate_id)
    url = f'https://storage.googleapis.com/pseudobulk-data/{plateid}_meta.h5ad'

    adata= load_h5ad_from_url(url)

    # Extracting relevant data from the file according to the user input
    adata_subset = extract_data(adata, drug, drug_conc, tumor_indication)
    dlist = adata_subset.obs['drug_name'].unique().to_list()
    #st.write(f"Batch: {str.title(plate_id)}, Shape: {adata_subset.shape}, Drug: {dlist[0]} and {dlist[1]} ")

    res_df = diff_exp(adata_subset, drug)  
    st.success("Analysis completed!")
    st.download_button("Download DEG results", res_df.to_csv(index=False), file_name="DESeq2_results.csv")
    st.write('Top 50 differentially expressed genes:')
    st.dataframe(res_df.head(50))
    st.write('Volcano plot:')
    fig = volcano_plot(res_df, treatment = drug, ref = 'DMSO', top_n=3)
    st.pyplot(fig)
    
    # Running pathway analysis on upregulated genes
    up_reg_genes = res_df[(res_df['significant'] == True) & (res_df['log2FoldChange'] > 1)]
    #st.dataframe(up_reg_genes)
    up, fig = pathway_analysis(up_reg_genes.index.to_list())
    st.write('Pathway analysis with upregulated genes:')
    st.dataframe(up.results.head(10))
    # Running pathway analysis on downregulated genes
    down_reg_genes = res_df[(res_df['significant'] == True) & (res_df['log2FoldChange'] < -1)]
    down, fig = pathway_analysis(down_reg_genes.index.to_list())
    st.write('Pathway analysis with downregulated genes:')
    st.dataframe(down.results.head(10))     
        
        
        
        
        
        
        
        
