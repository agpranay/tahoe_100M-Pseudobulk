# Pseudobulk analysis of Tahoe-100M Single-cell transcriptome dataset

This repository can be used to pseudobulk aggregation and analysis of the Tahoe-100M single-cell RNA-seq dataset from the Arc Institute. The single-cell expression data is aggregated by experimental conditions (cell line and drug combination). 

You can interactively perform differential gene expression analysis at:

https://tahoe-100m-13989002605.us-central1.run.app/

## Features

Streamlit based app for interactive differential gene expression analysis
Drug-related metadata includes drug target and its mechanism of action. It also has additional information for each drug from OpenAI with sources
Cell-related metadata contain tumor-type information
Sample metadata contains information extracted while transforming single-cell to pseudobulk (contains plate-id information, durg name and concentration and cell-line information   

## Data

You can access the pseudobulk-data file from Google Cloud Storage (bucket name: pseudobulk-data).


