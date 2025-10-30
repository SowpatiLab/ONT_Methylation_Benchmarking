## Issue Report
We wanted to use the new v5.2.0_5mC@v1 (v5.2r1 hereon) modbase models with dorado v1.0.0 to look at 5mC methylation and ran into some weird behavior. We ran both the 5mC and the 5mCG models, and we see expected CpG calls with  both, but the 5mC model isn't reporting expected non CpG calls. The v5.0.0_5mC_5hmC@v3 (v5r3 from hereon) models on the same datasets seem to be doing fine. Here are more details on how we tested this:

### Experimental design:
| dataset | description                              |
|---------|------------------------------------------|
| Ecoli_WT | Wild Type _*E.coli*_ expected to have 5mC in CmCWGG context                        |
| Ecoli_DM | A mutant strain of _*E.coli*_ which has no native 5mC |
| Ecoli_DM_MSssI | Genomic DNA of above mutant strain treated with M.SssI, expected to have 5mC at CpG sites | 

From prior literature as well as bisulfite data, we confirmed close to 100% methylation at expected sites (CCWGG in WT and CG in M.SssI treated DNA)

We ran the methylation calling models on these datasets as follows, and compiled the output using modkit v0.4.4:

```
# For v5.0.0_5mC@v3 (v5r3):
dorado-v0.9.1/bin/dorado basecaller \
   doradoModels/dna_r10.4.1_e8.2_400bps_sup@v5.0.0 \
   --min-qscore 10 \
   --reference <reference.fa> \
   --modified-bases-models doradoModels/dna_r10.4.1_e8.2_400bps_sup_sup@v5.0.0_5mC_5hmC@v3 \
   | samtools sort -@ 20 -O BAM -o $modbam --write-index;

# For v5.2.0_5mC@v1 (v5.2r1):
dorado-v1.0.0/bin/dorado basecaller \
   doradoModels/dna_r10.4.1_e8.2_400bps_sup@v5.2.0 \
   --min-qscore 10 \
   --reference <reference.fa> \
   --modified-bases-models doradoModels/dna_r10.4.1_e8.2_400bps_sup_sup@v5.2.0_5mC_5hmC@v1 \
   | samtools sort -@ 20 -O BAM -o $modbam --write-index;

modkit pileup $modbam $modbed -t 20 --ref <reference.fa> --ignore h 2;
```

The corresponding 5mCG_5hmCG models were also run on the same data in the same manner. We then counted the number of sites with >80% methylation. We saw that while both v5.2r1 and v5r3 5mC models report similar number of CpG sites in DM_MSssI data (as expected), only the v5r3 model reported methylated sites in the WT sample at CCWGG sites. We noticed that the default threshold used by modkit for v5.2r1 was too high (1), and that modkit v0.5.0 uses a fallback mechanism in such cases. But, even when we repeated the pileup with modkit v0.5.0, the results were similar. The methylated CCWGG sites reported in WT sample, while not 0 like modkit v0.4.4, was still significantly lower for v5.2r1 model.

The results, and the default thresholds that various models called are summarized below in the table, and also in the plot that follows.

### Default thresholds set by modkit, and number of reported sites

```
+--------------------------------------------------------------------------------------------------------------------------+
|                 |          |          |       |           modkit threshold    |    n sites > 80% methylated   |          |
|    dataset      | accuracy |  version | model |-------------------------------|-------------------------------|  motif   |
|                 |          |          |       | modkit_v0.4.4 | modkit_v0.5.0 | modkit_v0.4.4 | modkit_v0.5.0 |          |
|-----------------|----------|----------|-------|---------------|---------------|---------------|---------------|----------|
| Ecoli_WT        | hac      |  v5.2r1  | 5mC   | 1             | 0.9941406     | 0             | 392           |          |
| Ecoli_WT        | sup      |  v5.2r1  | 5mC   | 1             | 0.9941406     | 0             | 240           |          |
| Ecoli_WT        | hac      |  v5.2r1  | 5mCG  | 0.86816406    | 0.86816406    | 0             | 0             |          |
| Ecoli_WT        | sup      |  v5.2r1  | 5mCG  | 0.91308594    | 0.91308594    | 0             | 0             |  CCWGG   |
| Ecoli_WT        | hac      |  v5r3    | 5mC   | 0.72558594    | 0.72558594    | 23,619        | 23,619        |          |
| Ecoli_WT        | sup      |  v5r3    | 5mC   | 0.72753906    | 0.72753906    | 23,745        | 23,745        |          |
| Ecoli_WT        | hac      |  v5r3    | 5mCG  | 0.82910156    | 0.82910156    | 0             | 0             |          |
| Ecoli_WT        | sup      |  v5r3    | 5mCG  | 0.82910156    | 0.82910156    | 0             | 0             |          |
|-----------------|----------|----------|-------|---------------|---------------|---------------|---------------|----------|
| Ecoli_DM        | hac      |  v5.2r1  | 5mC   | 1             | 0.9941406     | 0             | 0             |          |
| Ecoli_DM        | sup      |  v5.2r1  | 5mC   | 1             | 0.9941406     | 0             | 0             |          |
| Ecoli_DM        | hac      |  v5.2r1  | 5mCG  | 0.8798828     | 0.8798828     | 0             | 0             |          |
| Ecoli_DM        | sup      |  v5.2r1  | 5mCG  | 0.91503906    | 0.91503906    | 0             | 0             |          |
| Ecoli_DM        | hac      |  v5r3    | 5mC   | 0.70410156    | 0.70410156    | 2             | 2             |  CCWGG   |
| Ecoli_DM        | sup      |  v5r3    | 5mC   | 0.7001953     | 0.7001953     | 1             | 1             |          |
| Ecoli_DM        | hac      |  v5r3    | 5mCG  | 0.8251953     | 0.8251953     | 0             | 0             |          |
| Ecoli_DM        | sup      |  v5r3    | 5mCG  | 0.82128906    | 0.82128906    | 0             | 0             |          |
|-----------------|----------|----------|-------|---------------|---------------|---------------|---------------|----------|
| Ecoli_DM_MSssI  | hac      |  v5.2r1  | 5mC   | 0.82128906    | 0.82128906    | 600,811       | 600,811       |          |
| Ecoli_DM_MSssI  | sup      |  v5.2r1  | 5mC   | 0.89160156    | 0.89160156    | 611,741       | 611,741       |          |
| Ecoli_DM_MSssI  | hac      |  v5.2r1  | 5mCG  | 0.6826172     | 0.6826172     | 652,291       | 652,291       |          |
| Ecoli_DM_MSssI  | sup      |  v5.2r1  | 5mCG  | 0.72558594    | 0.72558594    | 656,007       | 656,007       |          |
| Ecoli_DM_MSssI  | hac      |  v5r3    | 5mC   | 0.7373047     | 0.7373047     | 686,724       | 686,724       |    CG    |
| Ecoli_DM_MSssI  | sup      |  v5r3    | 5mC   | 0.7626953     | 0.7626953     | 683,520       | 683,520       |          |
| Ecoli_DM_MSssI  | hac      |  v5r3    | 5mCG  | 0.82910156    | 0.82910156    | 680,601       | 680,601       |          |
| Ecoli_DM_MSssI  | sup      |  v5r3    | 5mCG  | 0.7861328     | 0.7861328     | 666,124       | 666,124       |          |
+--------------------------------------------------------------------------------------------------------------------------+
```
<img width="50%" alt="Image" src="https://github.com/user-attachments/assets/55de3529-06ad-472d-803f-970ddaf45fc4" />

As you can see from above, the issue is consistent with both hac and sup models, as well as both modkit v0.4.4 and 0.5.0. Further, we rebasecalled these data with v5r3 models using Dorado v1.0, and we get correct results, ruling out issue with the Dorado version.

We then looked at the read-level probability scores using modkit extract. As you can see below, at CpG sites in DM_MSssI, both v5r3 and v5.2r1 models report most reads to be methylated with high confidence (though the overall read count is lower with v5.2r1). But, at CCWGG sites in WT, only ~30% of reads get classified as methylated and rest are called as unmethylated with high confidence.

Please let us know what we might be doing wrong. Other details of our run environment are below:

## Run environment
- Dorado version: v1.0.0
- Dorado command: dorado basecaller
- Operating system: linux x86 64
- Hardware specs: 128 cores, 256gb RAM, Nvidia L40s
- Source data type: pod5
- Source data location: HDD
- Flow cell version : R10.4.1
- Library prep kit: LSK114

## Reproducibility

A snakemakew workflow has been provided along with this repository. This can be used to replicate the data processing pipeline. [snakemake workflow](./snakemake_pipelines/snakemake.md)

The example workflow can be run by refering to the [tutorial.md](tutorial.md) file.