
conda create -n test \
    -c conda-forge \
    bioconda::snakemake \
    conda-forge::polars \
    conda-forge::seaborn \
    bioconda::samtools \
    bioconda::bedtools \
    bioconda::nanostat

pip install biopython blue-crab pod5

conda env create -f env_clean.yml -y
