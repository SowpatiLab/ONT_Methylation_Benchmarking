# This is a scaffold script that contains the steps used
# to perform the benchmark study. This script is not 
# intended to be used as it, the variables need to we 
# adapted to the specific sample name and usecase


########################################################
# basecalling with dorado
########################################################

# for 4kHz data
dorado-0.5.3-linux-x64/bin/dorado basecaller \
    ${rerio_model_directory}/res_dna_r10.4.1_e8.2_400bps_sup@v4.0.1 \
    ${input_pod5} \
    --modified-bases-models \
    ${rerio_model_directory}/res_dna_r10.4.1_e8.2_400bps_sup@v4.0.1_5mC@v2,${rerio_model_directory}/res_dna_r10.4.1_e8.2_400bps_sup@v4.0.1_6mA@v2 \
    > ${unaligned_modbam}_4kHz_v4.bam

# for 5kHz data v4 models
dorado-0.5.3-linux-x64/bin/dorado basecaller \
    ${rerio_model_directory}/res_dna_r10.4.1_e8.2_400bps_sup@v4.0.1 \
    ${input_pod5} \
    --modified-bases-models \
        ${dorado_model}/dna_r10.4.1_e8.2_400bps_sup@v4.3.0_6mA@v1,${rerio_model_directory}/res_dna_r10.4.1_e8.2_400bps_sup@v4.3.0_4mC_5mC@v1 \
    > ${unaligned_modbam}_5kHz_v4.bam

# for 5kHz data v5 models
dorado-0.7.0-linux-x64/bin/dorado basecaller ${input_pod5} --modified-bases-models > ${unaligned_modbam}.bam \
    ${dorado_model}/dna_r10.4.1_e8.2_400bps_sup@v5.0.0 \
    ${input_pod5} \
    --modified-bases-models \
        ${dorado_model}/dna_r10.4.1_e8.2_400bps_sup@v5.0.0_4mC_5mC@v1,${dorado_model}/dna_r10.4.1_e8.2_400bps_sup@v5.0.0_6mA@v1 \
    > ${unaligned_modbam}_5kHz_v4.bam

########################################################
# Converting unaligned modbam file to fastq
########################################################


samtools fastq -T '*' ${unaligned_modbam}.bam > ${mod_fastq}.fastq


########################################################
# aligning fastq reads to reference genome:
########################################################

minimap2 -ax map-ont -y $reference_genome ${mod_fastq}.fastq | samtools sort -O BAM -o ${aligned_modbam}.bam

########################################################
# nanostat
########################################################

NanoStat --bam ${aligned_modbam}.bam > ${nanostat_modbam}.bam


########################################################
# modkit | converting modbam files to per site bed files:
########################################################

modkit pileup ${aligned_modbam}.bam ${modbed}.bed --ref $reference_genome
awk '$4=="5mC"' ${modbed}.bed > ${modbed}_5mC.bed
awk '$4=="6mA"' ${modbed}.bed > ${modbed}_5mC.bed
awk '$4=="4mC"' ${modbed}.bed > ${modbed}_5mC.bed # 4mC is denoted as 21839

########################################################
# intersecting 5mC 
########################################################

bedtools intersect -a ${modbed}_5mC.bed -b ${bisulfite}_CpG.bed -wa -wb | cut -f 1-6,10-13,22-24 > ${modbed}_5mC_CpG_intersected.bed
bedtools intersect -a ${modbed}_5mC.bed -b ${bisulfite}_CHG.bed -wa -wb | cut -f 1-6,10-13,22-24 > ${modbed}_5mC_CHG_intersected.bed
bedtools intersect -a ${modbed}_5mC.bed -b ${bisulfite}_CHH.bed -wa -wb | cut -f 1-6,10-13,22-24 > ${modbed}_5mC_CHH_intersected.bed

########################################################
# adding fasta sequence to bed files
########################################################

samtools faidx ${reference}.fasta -o ${reference}.fasta.fai
cut -f 1,2 ${reference}.fasta.fai > ${reference}.genome


left=5   # number of nucleotides to expand to left
right=11 # number of nucleotides to expand to right

bedtools slop -i ${modbed}_5mC_${context}_intersected.bed -l ${left} -r ${right} -g ${reference}.genome -s |\
    bedtools getfasta -fi ${reference}.fasta -s -bed - -tab | \
    cut -f2 | paste ${modbed}_5mC_${context}_intersected.bed - > ${modbed}_5mC_${context}_intersected_fasta.bed


python reshape_modbed.py \
    -i ${modbed}_5mC_${context}_intersected_fasta.bed \
    -o ${modbed}_5mC_${context}_intersected_fasta_cleansed.bed \
    -c ${CONTEXT} \
    --species ${species} \
    --condition ${condition} \
    --model ${model} \
    --tool ${tool} \
    --sample-rate ${sample_rate};

