# Running the Nextflow workflow 

## Introduction

This folder contains the Nextflow implementation of the ont_methylation_benchmarking. 

## Setup
Since multiple tools are being benchnarked at once, there are multiple dependencies that need to be satisfied before the workflow can be run. These include:

- [Dorado](https://github.com/nanoporetech/dorado) - v0.9.1/v1.1.1 <br/>
- [DeepMod2](https://github.com/WGLab/DeepMod2/tree/main)<br/>
- [DeepBAM](https://github.com/xiaochuanle/DeepBAM) - v0.1.0<br/>
- [DeepPlant](https://github.com/xiaochuanle/DeepPlant/tree/main) - v1.0.0<br/>
- [f5C](https://github.com/hasindu2008/f5c/tree/master) - v1.5.0 <br/>
- [Rockfish](https://github.com/lbcb-sci/rockfish/tree/r10.4.1) - r10.4.1 branch<br/>
- [Mokit](https://github.com/nanoporetech/modkit) - v0.5.1-rc1<br/>

The current workflow can be configured to run with existing working installations of said tools, else if not detected, the tools will be installed by the workflow in the directory
mentioned in the [config.yaml](../config.yaml) file.

<b>Note:</b> At the time of writing, DeepBAM and DeepPlant do not provide any precompiled versions of their software and hence they would still have to be installed manually, and the location to their respective executable will have to be provided in the [config.yaml](../config.yaml) file.


### Pre-requisites to run the workflow 
1. [Nextflow](https://www.Nextflow.io/docs/latest/install.html)
2. [Anaconda](https://www.anaconda.com/download)

### Dataset and naming

```
{
    "pod5dir": "datasets/pod5"
}
```
#### File name nomenclature rule:
The pod5dir should have a folder for each unique experiment. The pod5 file(s) must me stored under the folder for each correspoding folder.

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

<pre>
{
    "experiments": [
        "Ecoli_WT",
        "Ecoli_DM",
        "Ecoli_DM_MSssI",
        "HP26695_WT",
        "HP26695_WGA",
    ],
}

<b>## this configuration will cause the workflow to expect their corresponding folders (with the '5kHz') suffix in 'datasets/pod5' and exclude the rest</b>
</pre>


#### directory setup:

```
datasets/pod5
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
```

### Reference files

The references for each species used must be included in the [references.yaml](../references.yaml). The reference locations are included with the 'references' key, while the 'std_chroms' can be used to control which chromosomes/contigues are filtered for (leaving it empty will return all available chromosomes/contigues).

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

### Tooling setup

The 'toolConfig' key contains a list of all the tools that the workflow can autoconfigure, including the conda environments names. Alternatively if any of these tools are already installed on your machine you can point Nextflow to those directories, but we recommend allowing Nextflow to install it fresh for the sake of starting from a fresh point.

The following keys need to be setup:
<pre> 
    {
        "<b>tooling_dir</b>": "tooling"
        "<b>toolConfig</b>": {
            "<b>tool_name</b>": {
                "<b>install_dir</b>": "rockfish",      
                "<b>executable</b>": "rockfish"        
                "<b>conda</b>": "rockfish_test_nf",    
                "<v><b>model_dir</b></v>": "rockfish_models",
                "<b>model</b>": "rf_5kHz.ckpt",        
                "<b>call_flags</b>": "-b 8192 -d 0"    
            }
        }
    }
    
</pre>

<b>tooling_dir: </b> This key points to where the tool will be installed. Ideally the workflow expects all the tools to be installed under the same location. If each of your tools are installed we recommend setting tooling_dir to "" and giving the absolute path of each tool to its corresponding install_dir key.

| key | description | 
|-|-|
| install_dir | If this key is present the tool will be installed by the workflow. It can also be changed to use the tool that is pre-installed in your system| 
| executable | This attribute is the path to the executable. If using your own install, please ensure the executable is functional| 
| conda | This directs the name of the conda environment name, and can be reused by the user since it gets installed globally. If a environment of the same name exists, it will be reinstalled by this workflow. |  
| model_dir | This is where the models for the given tool will be downloaded | 
| model | Model name |
| call_flags | this can be extended to add new arguments to a give tool |
    
    
## Running the workflow 

Asuming the workflow is contained in the <b>benchmark_nextflow</b>:

```
    nextflow run benchmark_nextflow \
        -param-file config.yaml
```

## Running the workflow with apptainer
Since DeepBAM and DeepPlant require sudo permissions to be installed, 
we recommend building a containerised image with apptainer with the 
[build config](../apptainer_setup.def) included in this repo for reproducibility sake.

### Setting up the apptainer image

```bash
    # from this directory:
    [ ! -d ../apptainer_build ] && mkdir ../apptainer_build;
    sudo apptainer build apptainer_build/ontMethylationBenchmarking.sif ../apptainer_setup.def
```
### running with apptainer

```
    nextflow run benchmark_nextflow \
        -param-file config.yaml \
        -with-apptainer
```

## Final folder structure 
```
├── bam                        # bam files go here
│   ├── sorted_mod               # mod bams
│   ├── sorted_mod_cleansed      # filtered mod bams
│   ├── sorted_move              # move table bams
│   ├── sorted_move_cleansed     # filtered move table bams
│   └── sorted_move_fnord        # move table bams sorted by pod5 id
├── intermediary               # intermediate files stored here
│   ├── blow5
│   └── fastq
├── meta                       # final output beds
│   ├── deepmod2
│   ├── deepbam
│   ├── deepplant
│   ├── f5c
│   ├── rockfish
│   ├── deepbam
│   └── dorado
├── tool_out
│   ├── deepmod2
│   ├── deepbam
│   ├── deepplant
│   ├── f5c
│   ├── rockfish
│   └── dorado
├── tooling                    # tool install dir
│   ├── DeepMod2
│   ├── dorado
│   ├── dorado_models
│   ├── f5c
│   ├── modkit
│   ├── rockfish
│   └── rockfishmodels
├── benchmark_nextflow         # nextflow workflow is contained here
├── config.yaml
└── references.yaml
```
