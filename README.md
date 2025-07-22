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

<table>
    <tr>
        <th>Sr</th>
        <th>Tool</th>
        <th>SampleRate</th>
        <th>Model</th>
        <th>Mods</th>
    </tr>
    <tr>
        <td rowspan="6">1</td>
        <td rowspan="6">Dorado</td>
        <td>4kHz</td>
        <td>res_dna_r10.4.1_e8.2_400bps_sup@v4.0.1_5mC@v2</td>
        <td>5mC</td>
    </tr>
    <tr>
        <td>5kHz</td>
        <td>dna_r10.4.1_e8.2_400bps_sup@v4.3.0_5mC_5hmC@v1 <br/>dna_r10.4.1_e8.2_400bps_sup@v4.3.0_5mCG_5hmCG@v1 <br/>res_dna_r10.4.1_e8.2_400bps_sup@v4.3.0_4mC_5mC@v1 <br/>dna_r10.4.1_e8.2_400bps_sup@v4.3.0_6mA@v2</td>
        <td> 5mC<br/> 5mCG<br/> 6mA<br/> 4mC</td>
    </tr>
    <tr>
        <td>5kHz</td>
        <td>dna_r10.4.1_e8.2_400bps_sup@v4.3.0_5mC_5hmC@v1 <br/>dna_r10.4.1_e8.2_400bps_sup@v4.3.0_5mCG_5hmCG@v1 <br/>res_dna_r10.4.1_e8.2_400bps_sup@v4.3.0_4mC_5mC@v1 <br/>dna_r10.4.1_e8.2_400bps_sup@v4.3.0_6mA@v2</td>
        <td> 5mC<br/> 5mCG<br/> 6mA<br/> 4mC</td>
    </tr>
    <tr>
        <td>5kHz</td>
        <td>dna_r10.4.1_e8.2_400bps_sup@v5.0.0_5mC_5hmC@v1<br/>dna_r10.4.1_e8.2_400bps_sup@v5.0.0_5mCG_5hmCG@v1<br/>dna_r10.4.1_e8.2_400bps_sup@v5.0.0_6mA@v1<br/>dna_r10.4.1_e8.2_400bps_sup@v5.0.0_4mC_5mC@v1</td>
        <td> 5mC<br/> 5mCG<br/> 6mA<br/> 4mC</td>
    </tr>
    <tr>
        <td>5kHz</td>
        <td>dna_r10.4.1_e8.2_400bps_sup@v5.0.0_5mC_5hmC@v3<br/>dna_r10.4.1_e8.2_400bps_sup@v5.0.0_5mCG_5hmCG@v3<br/>dna_r10.4.1_e8.2_400bps_sup@v5.0.0_6mA@v3<br/>dna_r10.4.1_e8.2_400bps_sup@v5.0.0_4mC_5mC@v3</td>
        <td> 5mC<br/> 5mCG<br/> 6mA<br/> 4mC</td>
    </tr>
    <tr>
        <td>5kHz</td>
        <td>dna_r10.4.1_e8.2_400bps_sup@v5.2.0_5mC_5hmC@v1<br/>dna_r10.4.1_e8.2_400bps_sup@v5.2.0_5mCG_5hmCG@v1<br/>dna_r10.4.1_e8.2_400bps_sup@v5.2.0_6mA@v1<br/>dna_r10.4.1_e8.2_400bps_sup@v5.2.0_4mC_5mC@v1</td>
        <td> 5mC<br/> 5mCG<br/> 6mA<br/> 4mC</td>
    </tr>
    <tr>
        <td>2</td>
        <td>DeepMod2</td>
        <td>5kHz</td>
        <td>5kHz_Transformer<br/>5kHz_BiLSTM</td>
        <td>5mCG</td>
    </tr>
    <tr>
        <td>3</td>
        <td>F5C</td>
        <td>5kHz</td>
        <td>-</td>
        <td>5mCG</td>
    </tr>
    <tr>
        <td>4</td>
        <td>Rockfish</td>
        <td>5kHz</td>
        <td>rf_5kHz.ckpt</td>
        <td>5mCG</td>
    </tr>
    <tr>
        <td>5</td>
        <td>DeepBAM</td>
        <td>5kHz</td>
        <td>LSTM_20240524_newfeature_script_b9_s15_epoch25_accuracy0.9742.pt</td>
        <td>5mCG</td>
    </tr>
    <tr>
        <td>5</td>
        <td>DeepPlant</td>
        <td>5kHz</td>
        <td>both_bilstm.b51_s15_epoch8.cpg<br/>  both_bilstm.b51_s15_epoch9.chg <br/> both_bilstm.b13_s15_epoch8.chh </td>
        <td>5mCG<br/>5mCHG<br/>5mCHH</td>
    </tr>



</table>

## Datasets Benchmarked

<table>
    <tr>
        <th></th>
        <th></th>
        <th>organism</th>
        <th>Sample</th>
    </tr>
    <tr>
        <td rowspan="5">Bacteria</td>
        <td rowspan="1">1<br/>2<br/>3</td>
        <td rowspan="1"><i>Escherichia coli</i> str. K-12 substr. MG1655</td>
        <td rowspan="1">Native (WT) <br/>Double Mutant (DM) <br/> Double Mutant M.SssI Treated (DM_M.SssI)</td>
    </tr>
    <tr>
        <td rowspan="1">4<br/>5</td>
        <td rowspan="1"><i>Helicobacter pylori</i> str. 26695</td>
        <td rowspan="1">Native (WT) <br/> Whole Genome Amplified (WGA)</td>
    </tr>
   <tr>
        <td rowspan="1">6</td>
        <td rowspan="1"><i>Helicobacter pylori</i> str. J99</td>
        <td rowspan="1">Native (WT)</td>
    </tr>
    <tr>
        <td rowspan="1">7</td>
        <td rowspan="1"><i>Anabaena variabilis</i> ATCC 27983</td>
        <td rowspan="1">Native (WT)</td>
    </tr>
    <tr>
        <td rowspan="1">8</td>
        <td rowspan="1"><i>Treponema denticola</i> ATCC 35405</td>
        <td rowspan="1">Native (WT)</td>
    </tr>
    <tr>
        <td rowspan="2">Mammalian</td>
        <td rowspan="1">9</td>
        <td rowspan="1">Human</td>
        <td rowspan="1">HG002</td>
    </tr>
    <tr>
        <td rowspan="1">10</td>
        <td rowspan="1">Mouse</td>
        <td rowspan="1">mouse_Brain <br/> mouse_ESC</td>
    </tr>
    <tr>
        <td rowspan="2">Plant</td>
        <td rowspan="1">11</td>
        <td rowspan="1"><i>Arabidopsis thaliana</i></td>
        <td rowspan="1">Native (WT)</td>
    </tr>
    <tr>
        <td rowspan="1">12</td>
        <td rowspan="1"><i>Oryza sativa japonica</i></td>
        <td rowspan="1">Native (WT)</td>
    </tr>
</table>



## Contact
In case of any queries/suggestions, contact

Onkar Kulkarni - onkar {at} ccmb {dot} res {dot} in <br>
Divya Tej Sowpati - tej {at} ccmb {dot} res {dot} in
