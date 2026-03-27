# Tutorial on an example dataset

In this tutorial, we provide step-by-step instructions to run the benchmarking pipeline. The necessary pod5 and references are provided in the example folder:<br/> 
- Pod5 files - [example/pod5](example/pod5) <br/>
- Reference - [example/references/H.pylori_J99_ATCC700824.fa](example/references/H.pylori_J99_ATCC700824.fa)

You can also run the pipeline using your own dataset by replacing the example dataset folder with your dataset folder and downloading the appropriate genome file from UCSC.

Creating conda environment

## Setting up the environment
1. `sudo apt update && sudo apt install openjdk-17-jdk git wget -y` (if on a debian machine)
2. Install [docker](https://www.digitalocean.com/community/tutorials/how-to-install-and-use-docker-on-ubuntu-22-04)
    - setup docker to run without sudo
        ```bash
        sudo groupadd docker
        sudo usermod -aG docker $USER
        newgrp docker
        ```
    - test if everything is working: 
        ```
            docker run hello-world
        ```
3. Install conda:
    ```
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
    bash miniconda.sh -b -u
    ~/miniconda3/bin/conda init bash 
    source ~/.bashrc 
    conda tos accept
    ```
4. Install [nextflow](https://www.nextflow.io/) 
    ```
    curl -s https://get.nextflow.io | bash
    echo 'alias nextflow='$(realpath nextflow) >> ~/.bashrc
    source ~/.bashrc
    ```
5. Pull github repo:
    ```
    git clone https://github.com/SowpatiLab/ONT_Methylation_Benchmarking.git 
    cd ONT_Methylation_Benchmarking
    ```
4. Index reference files
    All the reference files are expected to be indexed first.
    ```
        samtools faidx example/references/H.pylori_J99_ATCC700824.fa
        cut -f 1,2 \
            example/references/H.pylori_J99_ATCC700824.fa.fai \
            >  example/references/H.pylori_J99_ATCC700824.genome
    ```
    > **Note:** This is for the example reference file in the example folder. All other references need to be indexed in the same fashion.

4. Running the nextflow workflow:
    ```
    cd example

    # run nextflow with conda 
    nextflow run ../benchmark_nextflow \
        -params-file config.yaml \
        -profile conda --conda_prefix $CONDA_PREFIX \
        --threads 10 \
        --memory  32.GB

    # running with docker profile
    # only run once 
    docker pull sowpati/ont-methylation-benchmarking:latest
    [ ! -d tooling/rockfish_models ] && mkdir -p tooling/rockfish_models
    docker run -v $(pwd):$(pwd) -w $(pwd) \
        sowpati/ont-methylation-benchmarking:latest rockfish download \
            -m 5kHz -s tooling/rockfish_models

    # run nextflow with docker
    nextflow run ../benchmark_nextflow \
        -params-file config.yaml \
        -profile docker \
        --threads 10 \
        --memory  32.GB

    # run nextflow with apptainer
    nextflow run ../benchmark_nextflow \
        -params-file config.yaml \
        -profile apptainer \
        --threads 10 \
        --memory  32.GB
    ```
5. Installing [Snakemake](https://snakemake.readthedocs.io/en/stable/)
    ```
        # option 1:
        conda create \
            -c conda-forge \
            -c bioconda -c nodefaults \
            -n benchmark_env snakemake

        # option 2:
        conda env create \
            -n benchmark_env \
            -f benchmark_nextflow/envs/benchmarking_env.yml
        
        # activate conda env
        conda activate benchmark_env
    ```
6. Install apptainer:
    https://apptainer.org/

    on a debian machine:
   ```
    sudo add-apt-repository -y ppa:apptainer/ppa
    sudo apt update
    sudo apt install -y apptainer-suid
   ```

7. Running the snakemake workflow:
    ```
    ## ensure snakemake is installed 
    snakemake --use-apptainer \
        ../benchmark_snakemake/snakefile \
        --cores 10
    ```

## Basic dependencies required to run the workflow
- snakemake
- [samtools](https://www.htslib.org/) >=v1.16
- [bedtools](https://bedtools.readthedocs.io/en/latest/) 
- [pod5-api](https://pod5-file-format.readthedocs.io/en/0.1.21/docs/install.html)

> <b>Note: </b>all of these can be installed directly using the environment setup yml as follows:
```bash
    conda env create -f env.yml -n benchmark_nextflow/envs/benchmarking_env.yml
    conda activate benchmark_env
```

## Final folder structure 
<pre>
.
├── benchmark_nextflow         # nextflow  workflow lives here
├── benchmark_snakemake        # snakemake workflow lives here
├── output
|   │   ├── bam                # bam files go here
|   │   ├── sorted_mod               # mod bams
|   │   ├── sorted_mod_cleansed      # filtered mod bams
|   │   ├── sorted_move              # move table bams
|   │   ├── sorted_move_cleansed     # filtered move table bams
|   │   └── sorted_move_fnord        # move table bams sorted by pod5 id
│   ├── intermediary
|   │   ├── blow5
|   │   └── fastq
│   ├── meta                  # final output beds for each tool
|   │   ├── deepbam
|   │   ├── deepmod2
|   │   ├── deepplant
|   │   ├── dorado
|   │   ├── f5c
|   │   ├── f5c_stranded
|   │   └── rockfish
│   ├── qc
|   │   ├── mosdepth
|   │   ├── nanoplot
|   │   ├── nanoq
|   │   └── nanostat
│   └── tool_out
|       ├── deepbam
|       ├── deepmod2
|       ├── deepplant
|       ├── dorado
|       ├── f5c
|       └── rockfish
├── pod5                       # input raw signal files
├── references                 # references files go here
├── tooling                    # tool install dir
│   ├── DeepBAM_models
│   ├── deepmod2
│   ├── DeepPlant_models
│   ├── dorado
│   ├── dorado_models
│   ├── f5c
│   ├── modkit
│   ├── rockfish
│   └── rockfishmodels
├── config.yaml                # controls for the workflow
├── references.yaml            # reference file map
└── work                       # working folder generated by nextflow
</pre>

> **Note:** For more detailed usage of the workflow check the [nextflow.md](../benchmark_nextflow/readme.md) and [snakemake.md](benchmark_snakemake/readme.md) files 