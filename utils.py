#!/usr/bin/env python3

def volcano_plot(res_df, treatment="treated", ref="control", top_n=3):
    import numpy as np
    import pandas as pd
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
                     gene_sets='Reactome_Pathways_2024',
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
