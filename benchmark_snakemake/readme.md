# Running the snakemake workflow

## Dependencies

- [Dorado](https://github.com/nanoporetech/dorado)<br/>
- [DeepMod2](https://github.com/WGLab/DeepMod2/tree/main)<br/>
- [DeepBAM](https://github.com/xiaochuanle/DeepBAM)<br/>
- [DeepPlant](https://github.com/xiaochuanle/DeepPlant/tree/main)<br/>
- [f5C](https://github.com/hasindu2008/f5c/tree/master)<br/>
- [Rockfish](https://github.com/lbcb-sci/rockfish/tree/r10.4.1)<br/>

other dependencies can be installed using the [env.yml](/env.yml) included in this repo

```bash
    conda env create -f env.yml -n ont_basemod_benchmarking_env
    conda activate ont_basemod_benchmarking_env
```

## pre-requisites to workflow setup

### File name nomenclature rule:
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

These names can then be included under the "runControls" key in [config.json](../config.json).

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


### directory setup:

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

### setting up the config file

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

In the current configuration the model names mentioned in the square brackets represents the dorado model that each of the tool will use. Leaving any of the values empty would lead to that specific tool not running.

In this example: 
    only deepmod2_transformer and deepmod2_bilstm would run.
    and for dorado:
        v4r1, v5r1, v5r3, v5.2r2 would be run individually for 5mCG
        whicl 5mC, 4mC and 6mA would not be run.

#### Selecting which samples to run:
Once the pod5 files are in the designated folder (mentioned in config.json), the "exps" key (on line 5) can be used to control which samples to be run.
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

## Running the Snakemake workflow

```bash
    conda activate ont_basemod_benchmarking_env
    snakemake --use-conda -s snakemake_pipelines/snakefile --cores 40 
    ## cores can be set to the max that your machine supports
```
## final directory structure

The following is the directory structure generated by the snakemake workflow. Of which, those highlighted cyan are required by the workflow to start running while the others are generated by the workflow as it progresses through the pipeline.

<pre>
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
├── benchmark_snakemake        # snakemake workflow is contained here
├── config.json
└── references.yaml
</pre>