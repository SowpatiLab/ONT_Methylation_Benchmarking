# This is a scaffold script that contains the steps used
# to perform the benchmark study. This script is not 
# intended to be used as it, the variables need to we 
# adapted to the specific sample name and usecase


########################################################
# basecalling with dorado
########################################################

# for 4kHz data
dorado-0.5.3-linux-x64/bin/dorado basecaller \
    ${rero_model_directory}/res_dna_r10.4.1_e8.2_400bps_sup@v4.0.1 \
    ${input_pod5} \
    --modified-bases-models \
    --emit-moves \ 
    > ${unaligned_bam}_4kHz_v4.bam

dorado-0.5.3-linux-x64/bin/dorado basecaller \
    ${rero_model_directory}/res_dna_r10.4.1_e8.2_400bps_sup@v4.0.1 \
    ${input_pod5} \
    --modified-bases-models \
    --emit-moves \ 
    > ${unaligned_bam}_5kHz_v4.bam

dorado-0.7.0-linux-x64/bin/dorado basecaller ${input_pod5} --modified-bases-models > ${unaligned_bam}.bam \
    ${dorado_model}/dna_r10.4.1_e8.2_400bps_sup@v5.0.0 \
    ${input_pod5} \
    --modified-bases-models \
    --emit-moves \ 
    > ${unaligned_bam}_5kHz_v5.bam

########################################################
# rebase basecalling with deepmod2
########################################################

samtools fastq ${unaligned_bam}.bam -T MM,ML,mv,ts | \
    minimap2 -ax map-ont  $reference_genome - -y | \
    samtools sort -o ${aligned_modbam}.bam --write-index

# BiLSTM 4kHz 
python ${PATH_TO_DEEPMOD2_REPOSITORY}/deepmod2 detect \
    --bam ${aligned_modbam}.bam \
    --input INPUT_DIR \
    --model models/bilstm/R10.4.1_4kHz_v4.1 \
    --file_type FILE_TYPE \
    --ref $reference_genome \
    --output ${deepmod2_out}_bilstm

# BiLSTM 5kHz 
python ${PATH_TO_DEEPMOD2_REPOSITORY}/deepmod2 detect \
    --bam ${aligned_modbam}.bam \
    --input INPUT_DIR \
    --model models/bilstm/R10.4.1_5kHz_v4.3 \
    --file_type FILE_TYPE \
    --ref $reference_genome \
    --output ${deepmod2_out}_bilstm


# Transformer 4kHz 
python ${PATH_TO_DEEPMOD2_REPOSITORY}/deepmod2 detect \
    --bam ${aligned_modbam}.bam \
    --input INPUT_DIR \
    --model models/transformer/R10.4.1_4kHz_v4.1 \
    --file_type FILE_TYPE \
    --ref $reference_genome \
    --output ${deepmod2_out}_transformer

# Transformer 5kHz 
python ${PATH_TO_DEEPMOD2_REPOSITORY}/deepmod2 detect \
    --bam ${aligned_modbam}.bam \
    --input INPUT_DIR \
    --model models/transformer/R10.4.1_5kHz_v4.3 \
    --file_type FILE_TYPE \
    --ref $reference_genome \
    --output ${deepmod2_out}_transformer

########################################################
# reordering deepmod2 output
########################################################

bedtools intersect \
    -a ${deepmod2_out}/output.per_site.aggregated \
    -b ${bisulfite}_CpG.bed -wa -wb | awk -v OFS="\t" '{print  $1,$2,$3,"m",$5,$4,$5,$8,$6,$7,$20,$21,$22}' \
    > ${deepmod2_out}_modbef_intersected.bed

left=5   # number of nucleotides to expand to left
right=11 # number of nucleotides to expand to right

bedtools slop -i ${deepmod2_out}_modbef_intersected.bed -l ${left} -r ${right} -g ${reference}.genome -s |\
    bedtools getfasta -fi ${reference}.fasta -s -bed - -tab | \
    cut -f2 | paste ${deepmod2_out}_modbef_intersected.bed - > ${deepmod2_out}_5mC_${context}_intersected_fasta.bed

python reshape_modbed.py \
    -i ${deepmod2_out}_5mC_${context}_intersected_fasta.bed \
    -o ${deepmod2_out}_5mC_${context}_intersected_fasta_cleansed.bed \
    -c ${CONTEXT} \
    --species ${species} \
    --condition ${condition} \
    --model ${model} \
    --tool ${tool} \
    --sample-rate ${sample_rate};
