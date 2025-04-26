#!/usr/bin/env python3

import os
import sys
import pandas as pd
import scanpy as sc
import numpy as np
import anndata as ad
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr
import seaborn as sns


TOI_1 = sys.argv[1]
conc_1 = sys.argv[2]
tumor_indication = sys.argv[3]

if tumor_indication == 'All':
    LOI = None
else:
    LOI = tumor_indication
 
# Loading the complete dataset       
adata = sc.read_h5ad("/home/pranayagarwal/Documents/tahoe_100M/Data/all_plates.h5ad", backed = 'r')
# Subsetting the data based on user input
# Keep none if you would like to use all cell lines. Otherwise choose tumor indication of your choice
LOI = None # tumor indication of interest
TOI_1 = TOI_1 # 1st treatment of interest
conc_1 = conc_1 # concentration
# If you want to compared to different treatments

TOI_2 = None # 2nd treatment of interest
conc_2 = None

# Start with a universal query (all True)
query = pd.Series(True, index=adata.obs_names)

# Apply filters progressively if each variable is set
if LOI is None and TOI_1 is not None and TOI_2 is not None:
    query &= ((adata.obs["drug_name"] == TOI_1) & (adata.obs["concentration, uM"] == conc_1))
    query |= ((adata.obs["drug_name"] == TOI_2) & (adata.obs["concentration, uM"] == conc_2))
if LOI is not None and TOI_1 is not None and TOI_2 is not None:
    query &= ((adata.obs["lineage"] == LOI) & (adata.obs["drug_name"] == TOI_1) & (adata.obs["concentration, uM"] == conc_1))
    query |= ((adata.obs["lineage"] == LOI) & (adata.obs["drug_name"] == TOI_2) & (adata.obs["concentration, uM"] == conc_2))

if LOI is None and TOI_1 is None and TOI_2 is None:
    query &= (adata.obs["drug_name"] == 'DMSO_TF')
else:
    if LOI is not None and TOI_2 is None:
        query &= adata.obs["lineage"] == LOI

    if TOI_1 is not None and TOI_2 is None:
        query &= (adata.obs["drug_name"] == TOI_1)
            
    if conc_1 is not None and TOI_2 is None:
        query &= adata.obs["concentration, uM"] == conc_1

# Add DMSO treated samples from the same plate as the treatment group
if TOI_2 is None:
    plates = adata.obs[(adata.obs["drug_name"] == TOI_1) & (adata.obs["concentration, uM"] == conc_1)]['batch'].unique()        
    for plate in plates:
        if LOI is None:
            query |= ((adata.obs["drug_name"] == 'DMSO_TF') & (adata.obs["batch"] == plate))
        else:
            query |= ((adata.obs["lineage"] == LOI) & (adata.obs["drug_name"] == 'DMSO_TF') & (adata.obs["batch"] == plate))

adata_subset = adata[query].to_memory().copy()


adata_subset.layers['counts'] = adata_subset.X.copy()

# Compute mean & variance for each gene
adata_subset.var["mean"] = np.array(adata_subset.X.mean(axis=0)).flatten()
adata_subset.var["variance"] = np.array(adata_subset.X.var(axis=0)).flatten()

# Compute dispersion (variance/mean)
adata_subset.var["dispersion"] = adata_subset.var["variance"] / (adata_subset.var["mean"] + 1e-8)  # Avoid division by zero

# Sort by dispersion to find highly variable genes
adata_subset.var["highly_variable"] = adata_subset.var["dispersion"].rank(ascending=False) <= 5000  # Select top 5000 genes


pandas2ri.activate()

# Prepare data for performing differential gene expression analysis in R
r_counts = pandas2ri.py2rpy(pd.DataFrame(adata_subset.X.T, index=adata_subset.var_names, columns=adata_subset.obs_names))
r_coldata = pandas2ri.py2rpy(adata_subset.obs)

# Setting up the contrast
group=list(adata_subset.obs["drug_name"].cat.categories)
if 'DMSO_TF' in group:
    ref = 'DMSO_TF'
    treatment = [g for g in group if g != ref]
    treatment = treatment[0]
else:
    ref = group[0]
    treatment = group[1]   

print(f"Reference: {ref} and Treatment: {treatment}")
# Send to R
robjects.globalenv["countData"] = r_counts
robjects.globalenv["colData"] = r_coldata
robjects.globalenv["ref"] = ref
robjects.globalenv["treatment"] = treatment


# Run DESeq2 in R
robjects.r('''

# Load DESeq2 inside R
library(DESeq2)

# Create DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design = ~ drug_name)

# Apply variance stabilization transformation (VST).

vsd <- vst(dds, blind = TRUE)
vst_mat <- assay(vsd)
                              
# Filter low-count genes, if number of counts per row is less than 10

dds <- dds[rowSums(counts(dds)) > 10, ]

# Run DESeq2 differential expression analysis

dds <- DESeq(dds)

contrast <- c("drug_name", treatment, ref)

res_df <- as.data.frame(results(dds, contrast = contrast))

res_df$significant <- res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1

''')


# Convert results to pandas DataFrame
res_df = pandas2ri.rpy2py(robjects.r('as.data.frame(res_df)')).sort_values('padj')
res_df.to_csv(f"/home/pranayagarwal/Documents/tahoe_100M/Results/{treatment}_{conc_1}_vs_{ref}_in_{tumor_indication}.csv")
#vsd_df = pandas2ri.rpy2py(robjects.r('as.data.frame(vst_mat)'))
#adata_subset.X = vsd_df.T # VST_transformed values are in the main slot

