# conda activate ngs_data_charmer_copy
snakemake --snakefile Snakefile --cores "$1"
# conda deactivate
