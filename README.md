# ONT Methylation Benchmarking
This repository contains the scripts and code that was used for benchmarking various tools and models for Oxford Nanopore (ONT) sequencing based identification of DNA methylation. The corresponding preprint is [here](https://doi.org/10.1101/2024.11.09.622763)

## Contents
The code is organized into two main folders:
* Processing scripts
* Plotting scripts

### Processing scripts
This folder contains all the code that was required to convert variout modkit/bismark outputs into dataframes that can be used for comparison. In addition, there are also scripts that were used to subsample the data to various coverages, and to filter out reads under a specific q score.

### Plotting scripts
These are the code snippets that show how the plots used in the preprint were generated.

## Tools benchmarked:
|sr| Tool     | Model            | mods profiled | 
|--|----------|------------------|---------------|
| 1| Dorado   | 4kHz_v4          | 5mC, 6mA      |
| 2| Dorado   | 5kHz_v4          | 5mC, 6mA, 4mC |
| 3| Dorado   | 5kHz_v5          | 5mC, 6mA, 4mC |
| 4| DeepMod2 | 4kHz_Transformer | 5mC           |
| 5| DeepMod2 | 4kHz_BiLSTM      | 5mC           |
| 6| DeepMod2 | 5kHz_Transformer | 5mC           |
| 7| DeepMod2 | 5kHz_BiLSTM      | 5mC           |

## Datasets Benchmarked
|sr|Organism| Treatment |
|--|--------|-----------|
|1 | Human  HG002 | Native |
|2 | <i>Escherichia coli</i> str. K-12 substr. MG1655 | Native (WT) <br/> Double Mutant (DM)  <br/>  Double Mutant M.SssI Treated (DM_M.SssI) <br/>  Double Mutant M.CviP Treated (DM_M.CvIP)
|3 | <i>Helicobacter pylori</i> str. 26695 | Native (WT) <br/> Whole Genome Amplified (WGA)


## Contact
In case of any queries/suggestions, contact

Onkar Kulkarni - onkar {at} ccmb {dot} res {dot} in <br>
Divya Tej Sowpati - tej {at} ccmb {dot} res {dot} in
