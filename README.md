# podocoryna
Podocoryna single cell analysis


There are two snakemake pipelines in this repository: one including all the steps for reference optimization (in folder ref_optimization) and a second for the main analysis (main).

To repeat reference optimization pipeline with your own data:

1. Create a folder resources/rawdata containing folders named by sample including your single cell reads.
2. In `config/download_targets.tsv`, include species name and sources for taxa you would like to include in your gene trees.
3. In the `workflow/Snakefile`, include your sample IDs on line 7
4. Change species name in `config.yaml`.
5. Cellranger needs to be installed separately and made available on your PATH. See intructions for download here: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/tutorial_in
