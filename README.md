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
|sr| Tool     | Model            | 
|--|----------|------------------|
| 1| Dorado   | 4kHz_v4          |
| 2| Dorado   | 5kHz_v4          |
| 3| Dorado   | 5kHz_v5          |
| 4| DeepMod2 | 4kHz_Transformer |
| 5| DeepMod2 | 4kHz_BiLSTM      |
| 6| DeepMod2 | 5kHz_Transformer |
| 7| DeepMod2 | 5kHz_BiLSTM      |

## Datasets Benchmarked
|sr|Organism|
|--|--------|
|1 | Human  HG002 |
|2 | <i>Escherichia coli</i> str. K-12 substr. MG1655 |
|3 | <i>Helicobacter pylori</i> str. 26695 |


## Contact
In case of any queries/suggestions, contact

Onkar Kulkarni - onkar@ccmb.res.in <br>
Divya Tej Sowpati - tej@ccmb.res.in
