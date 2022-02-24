# Beta-cell_GSMM
Metabolic model of pancreatic β-cells in T2D

Codes_DEG_analysis.m conatins code for obtaining differntially expressed genes, which requires two files: Recon2.v04.mat and meta_gene_exp.mat.

Codes_Constructing_models.m contains code for construting metabolic models, which requires two files: Recon2.v04.mat and meta_gene_exp.mat.

Codes_Flux_distribution.m contains code for predicting metabolic flux rates, which requires the output models obatined from the codes: Codes_Constructing_models.m.

Codes_Comparison_fluxRates.m contains code for the comparison of flux rates which requires four files: Recon2.v04.mat, reaction_description.mat, flux_dia.mat and flux_non_dia.mat.
flux_dia.mat and flux_non_dia.mat files obatine from running the code provided in the file: Codes_Flux_distribution.m.

Codes_task_check.m contains code for checking the metabolic tasks which reqries two files: reaction_description.mat and Tasks.xlsx

Codes_Reporter_pathway.m code for performing reporter pathway analysis which reqries tive files: Recon2.v04.mat, reaction_description.mat, flux_p_val.mat, up_reaction.mat and down_reaction.mat
flux_p_val.mat, up_reaction.mat and down_reaction.mat files can be obatined from the code: Codes_Comparison_fluxRates.m

Code_WGCNA.r code for performing WGCNA which reqires to additional files: Data_coexpress.xlsx and Traits.xlsx.
