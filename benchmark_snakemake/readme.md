# Running the snakemake workflow

## Introduction

This folder contains the Snakemake implementation of the ont_methylation_benchmarking. 

## Setup
Since multiple tools are being benchnarked at once, there are multiple dependencies that need to be satisfied before the workflow can be run. These include:

### Base dependencies (have to be installed)
- [snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)
- [conda](https://www.anaconda.com/docs/getting-started/miniconda/install/overview#windows-guides) or [miniconda](https://www.anaconda.com/docs/getting-started/miniconda/main)
- [apptainer](https://apptainer.org/)

### Other tooling dependencies
- [biopython](https://biopython.org/)
- [polars](https://pola.rs/)
- [samtools](http://www.htslib.org/)
- [bedtools](https://bedtools.readthedocs.io/)
- [nanoq](https://github.com/esteinig/nanoq)
- [nanoplot](https://github.com/wdecoster/nanoplot)
- [nanostat](nanostat)

> These can be installed using the [benchmark_snakemake/env/benchmarking_env.yml](benchmark_snakemake/env/benchmarking_env.yml) file.

### Tooling dependencies (optional)
- [Dorado](https://github.com/nanoporetech/dorado) - v0.9.1/v1.1.1 <br/>
- [DeepMod2](https://github.com/WGLab/DeepMod2/tree/main)<br/>
- [DeepBAM](https://github.com/xiaochuanle/DeepBAM) - v0.1.0<br/>
- [DeepPlant](https://github.com/xiaochuanle/DeepPlant/tree/main) - v1.0.0<br/>
- [f5C](https://github.com/hasindu2008/f5c/tree/master) - v1.5.0 <br/>
- [Rockfish](https://github.com/lbcb-sci/rockfish/tree/r10.4.1) - r10.4.1 branch<br/>
- [Mokit](https://github.com/nanoporetech/modkit) - v0.5.1-rc1<br/>

Unlike the nextflow implementstion, the dependedencies will not get automatically installed by this workflow, and hence expects the user to preinstall the dependencies and mention the appropriate install directories and executable names in the [config.yaml](../config.yaml) file.

Alternatively, by using the `--use-apptainer` when running the workflow will auto fetch a pre-packaged [docker-image](https://hub.docker.com/r/sowpati/ont-methylation-benchmarking/tags) that contains all the necessary dependencies for this workflow to run successfully.

or the tools can be installed by running the nextflow workflow using the '--install' flag: 
`nextflow run ../benchmark_nextflow/ --install -params-file config.yaml `, this installs all the tools locally, including their respective conda environments (except for DeepBAM and DeepPlant)

## pre-requisites to workflow setup

## setting up the config file(s)

### Dataset and naming
inside [config.yaml](../config.yaml)
```yaml
pod5dir: pod5
```

#### File name nomenclature rule:
The pod5dir should have a folder for each unique experiment. The pod5 file(s) must me stored under the folder for each corresponding folder.

Folder names must contain the '_5kHz' suffix, this enables nextflow to exclude folder that do not have the suffix. 

<pre>
<b>e.g.:</b>

      ┌-------------------> <b>species name </b>
    |‾‾‾|           # this must be a single word with no 
    |   |              underscores in between
    |   |
    |   | ┌------┬--------> <b>Treatment name</b> # can be underscore delimited
    |   | |      |
    <b>Ecoli_DM_MSssI_5kHz</b> <-------<b>file name</b>
    |            | |  |
    |____________| |__|
        |            |
        |            └----><b>Sample rate</b>
        |
        └-----------------> <b>Experiment name</b>
                             # The experiment name is a combination of 
                                the species name and the treatment
                                eg: 
                                     1. Ecoli_WT - Ecoli - Wild type
                                     2. Ecoli_DM - Ecoli - double mutant
                                     3. Ecoli_DM - Ecoli - double mutant 
                                                          + MSssI treated
</pre>

These names can then be included under the "runControls" key in [config.yaml](../config.yaml).

inside [config.yaml](../config.yaml)
```yaml
experiments: [
    Ecoli_WT,
    Ecoli_DM,
    Ecoli_DM_MSssI,
    HP26695_WT,
    HP26695_WGA
]

## this configuration will cause the workflow to expect their corresponding folders (with the '5kHz') suffix in 'pod5' and exclude the rest
```
#### Selecting which samples to run:
Once the pod5 files are in the designated folder (mentioned in config.yaml), the "exps" key (on line 5) can be used to control which samples to be run.
<pre>
eg: 
    pod5 folder strurture:
        .
        ├── Anabaena_WT_5kHz
        ├── Ecoli_DM_5kHz
        ├── Ecoli_DM_MSssI_5kHz
        ├── Ecoli_WT_5kHz
        ├── HP26695_WGA_5kHz
        ├── HP26695_WT_5kHz
        ├── HPJ99_WT_5kHz
        └── Tdenticola_WT_5kHz
    
    For this setup, to run all the samples, the exps value has to
    be populated like so:
        
        exps: [Anabaena_WT,Ecoli_DM, Ecoli_DM_MSssI, Ecoli_WT,
                HP26695_WGA, HP26695_WT, HPJ99_WT, Tdenticola_WT]

        These are the names of the pod5 folders excluding 
        the '_5kHz' part.

        To limit the number of samples run, you can remove the 
        values that you dont want to be processed from the list.

        exps: [Ecoli_DM, Ecoli_WT]
                here only Ecoli_DM and Ecoli_WT will be processed.
</pre>

#### directory setup:

<pre>
pod5
├── Anabaena_WT_5kHz
│   └── Anabaena_WT_5kHz.pod5
├── Ecoli_DM_5kHz
│   └── Ecoli_DM_5kHz.pod5
├── Ecoli_DM_MSssI_5kHz
│   └── Ecoli_DM_MSssI_5kHz.pod5
├── Ecoli_WT_5kHz
│   └── Ecoli_WT_5kHz.pod5
├── HP26695_WGA_5kHz
│   └── HP26695_WGA_5kHz.pod5
├── HP26695_WT_5kHz
│   └── HP26695_WT_5kHz.pod5
├── HPJ99_WT_5kHz
│   └── HPJ99_WT_5kHz.pod5
└── Tdenticola_WT_5kHz
    └── Tdenticola_WT_5kHz.pod5

## multi pod5's are also accepted as long as they are included under the same folder
</pre>

### Reference files
inside [config.yaml](../config.yaml)
```yaml
references: references.yaml
```
the `references` key points to the  [references.yaml](../references.yaml) file that contains the paths to the reference genomes.

The references for each species used must be included in the [references.yaml](../references.yaml). While the `std_chroms` key can be used to control which chromosomes/contigues are filtered for (leaving it empty will return all available chromosomes/contigues).

```yaml
    references:
        Ecoli:    references/E.coli_MG1655.fa
        HP26695:  references/H.pylori_26695.fa
        HPJ99:    references/H.pylori_J99_ATCC700824.fa
        Nostoc:   references/Nostoc_sp_PCC7120_ATCC27893.fa
        Tdenticola:  references/Treponema_denticola_ATCC35405.fa
        HG002:    references/hg38.fa
        osjaponica: references/GCF_034140825.1_ASM3414082v1_genomic.fa
        arabidopsis: references/arabidopsis_thaliana.fa
        mouse:    references/mm39.fa

    std_chroms:
        HG002: [chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22]
        mouse: [chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19]
        osjaponica: [NC_089035.1,NC_089036.1,NC_089037.1,NC_089038.1,NC_089039.1,NC_089040.1,NC_089041.1,NC_089042.1,NC_089043.1,NC_089044.1,NC_089045.1,NC_089046.1]
        arabidopsis: [NC_003070.9,NC_003071.7,NC_003074.8,NC_003075.7,NC_003076.8]
        Ecoli: ['NC_000913.3']
        HP26695: ['NC_018939.1']
        HPJ99: ['NZ_CP011330.1']
        Nostoc: ['NC_003272.1']
        Tdenticola: ['NC_002967.9']
```

</br>

> **Note:** before the nextflow workflow is run, all reference files used are required to be indexed and their corresponding genome files have to be generated like so:
</br>

```bash
    samtools faidx reference_file.fa
    cut -f 1,2 reference_file.fa.fai >  reference_file.genome
```

#### selecting tools to run and their sorresponding models:
The config file contains the controls to which tools need to be run and where the install directories to each tool resides. This can be modified if you have already installed any of the aforementioned tools.

<pre>
    dorado_mods:
      5mC:  []
      5mCG: [v4r1, v5r1, v5r3, v5.2r2]
      6mA:  []
      4mC:  []
    
    rockfish_run:             []
    deepmod2_transformer_run: [v5r3]
    deepmod2_bilstm_run:      [v5r3]
    f5c_run:                  []
    f5c_stranded_run:         []
    deepbam_run:              []
    deepplnat_run:            []
</pre>

### Tooling setup

The `toolConfig` key contains a list of all the tools that the workflow can auto-configured, including the conda environments names. Alternatively if any of these tools are already installed on your machine you can point Nextflow to those directories, but we recommend allowing Nextflow to install it fresh for the sake of starting from a fresh point.

The following keys need to be setup:
```yaml
tooling_dir: tooling
toolConfig:
    tool_name:
        install_dir: rockfish,      
        executable: rockfish        
        conda: rockfish_test_nf,    
        model_dir: rockfish_models,
        model: rf_5kHz.ckpt,        
        call_flags: -b 8192 -d 0    
    
```

`tooling_dir` This key points to where the tool will be installed. Ideally the workflow expects all the tools to be installed under the same location. If each of your tools are already installed, we recommend setting tooling_dir to "" and giving the absolute path of each tool to its corresponding `install_dir` key.

| key | description | 
|-|-|
| `install_dir` | If this key is present the tool will be installed by the workflow. It can also be changed to use the tool that is pre-installed in your system| 
| `executable` | This attribute is the path to the executable. If using your own install, please ensure the executable is functional| 
| `conda` | This directs the name of the conda environment name, and can be reused by the user since it gets installed globally. If a environment of the same name exists, it will be reinstalled by this workflow. |  
| `model_dir` | This is where the models for the given tool will be downloaded | 
| `model` | Model name |
| `call_flags` | this can be extended to add new arguments to a give tool |

#### Adapting the workflow to your own system architecture
By default this workflow provides the following profiles:
1. conda
2. apptainer

In the current configuration the model names mentioned in the square brackets represents the dorado model that each of the tool will use. Leaving any of the values empty would lead to that specific tool not running.

In this example: 
    only deepmod2_transformer and deepmod2_bilstm would run.
    and for dorado:
        v4r1, v5r1, v5r3, v5.2r2 would be run individually for 5mCG
        whicl 5mC, 4mC and 6mA would not be run.


## Running the Snakemake workflow

```bash
    
    ## setup conda env
    conda env create \
        -n benchmark_env \
        -f ./benchmark_snakemake/env/benchmarking_env.yml

    conda activate benchmark_env
    snakemake \
        -s benchmark_snakemake/snakefile \
        --use-conda \
        --cores 10
    ## cores can be set to the max that your machine supports
```

## Running the Snakemake workflow with apptainer

```bash
    snakemake \
        -s benchmark_snakemake/snakefile \
        --use-apptainer \
        --apptainer-args "--nv " \
        --cores 10
```

## final directory structure

The following is the directory structure generated by the snakemake workflow. Of which, those highlighted cyan are required by the workflow to start running while the others are generated by the workflow as it progresses through the pipeline.

<pre>
.
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
