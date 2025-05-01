#!/usr/bin/env python3


def extract_data(adata, TOI, conc, tumor_indication):
    import os
    import sys
    import pandas as pd
    import scanpy as sc
    import numpy as np
    import anndata as ad
    
    if tumor_indication == 'All':
        LOI = None
    else:
        LOI = tumor_indication
 
    # Keep none if you would like to use all cell lines. Otherwise choose tumor indication of your choice
    # Start with a universal query (all True)
    query = pd.Series(True, index=adata.obs_names)

    # Apply filters progressively if each variable is set
    if LOI is None:
        query &= ((adata.obs["drug_name"] == TOI) & (adata.obs["concentration, uM"] == conc))
    
    if LOI is not None:
        query &= ((adata.obs["lineage"] == LOI) & (adata.obs["drug_name"] == TOI) & (adata.obs["concentration, uM"] == conc))

    
    # Add DMSO treated samples from the same plate as the treatment group
    if LOI is None:
        query |= ((adata.obs["drug_name"] == 'DMSO_TF'))
    else:
        query |= ((adata.obs["lineage"] == LOI) & (adata.obs["drug_name"] == 'DMSO_TF'))

    adata_subset = adata[query].to_memory().copy()


    adata_subset.layers['counts'] = adata_subset.X.copy()

    return adata_subset
    
def diff_exp(adata_subset, TOI):
    import pandas as pd
    import numpy as np
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
    from pydeseq2.default_inference import DefaultInference
    
    # Extract counts matrix
    counts = pd.DataFrame(adata_subset.layers['counts'], columns=adata_subset.var_names, index=adata_subset.obs_names)

    # Extract metadata
    metadata = adata_subset.obs.copy()
    metadata = metadata[metadata['drug_name'].isin(['DMSO_TF', TOI])]
    metadata = metadata.rename(columns={'drug_name': 'condition'})

    genes_to_keep = counts.columns[counts.sum(axis=0) >= 10]
    counts = counts[genes_to_keep]

    inference = DefaultInference(n_cpus=8)
    dds = DeseqDataSet(counts=counts, metadata=metadata, refit_cooks=True, inference = inference)

    dds.deseq2()
    
    ds = DeseqStats(dds, contrast=["condition", TOI, "DMSO_TF"])
    
    ds.summary()
    
    res_df = ds.results_df
    res_df['significant'] = (res_df['padj'] < 0.05) & (np.abs(res_df['log2FoldChange']) > 1)
    return res_df.sort_values('padj')

def volcano_plot(res_df, treatment="treated", ref="control", top_n=3):
    import numpy as np
    import pandas as pd
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    
    fig, ax = plt.subplots(figsize=(5, 5))

    sig = res_df['significant'] == True
    ax.scatter(res_df.loc[sig, 'log2FoldChange'], -np.log10(res_df.loc[sig, 'padj']),
               marker='s', alpha=0.5, color='red', label='Significant')
    ax.scatter(res_df.loc[~sig, 'log2FoldChange'], -np.log10(res_df.loc[~sig, 'padj']),
               marker='s', alpha=0.2, color='grey', label='Not Significant')

    top_genes = res_df[sig].sort_values('padj').head(top_n)
    '''
    for i, row in top_genes.iterrows():
     
        ax.text(row['log2FoldChange'], -np.log10(row['padj']), '  ' +top_genes.iloc[i, 0], color='red',
                ha='left', va='top')
'''
    ax.axhline(-np.log10(0.05), linestyle='--', color='black', linewidth=0.5)
    ax.axvline(0, linestyle='--', color='black', linewidth=0.5)

    ax.set_xlabel("Log2(Fold Change)")
    ax.set_ylabel("-Log10(Padj)")
    ax.set_title(f"{treatment.title()} vs {ref.title()}")
    return fig

def pathway_analysis(gene_list):
    import gseapy as gp
    import pandas as pd
    import matplotlib.pyplot as plt
    enr = gp.enrichr(gene_list=gene_list,
                     gene_sets='Data/ReactomePathways.gmt',
                     organism = 'human',
                     outdir = None
                    )
    
    fig, ax = plt.subplots(figsize=(5,5))
    # Plot
    i= min(enr.results.shape[0], 5)  # Limit to top 5 pathways
    ax.barh(range(i), enr.results['Combined Score'][:i], color='grey', alpha=0.5)
    ax.set_xlabel('Score')
        
    # Optionally, add pathway names as text labels
    for j in range(i):
        ax.text(0.1, i, enr.results['Term'][j], fontsize=12, color='red')

    return enr, fig
