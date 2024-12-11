# This is a scaffold script that contains the steps used to perform the subsampling of dorado output based on coverage.
# This script is not intended to be used as is.
# The variables need to be adapted to the specific sample name and usecase.


######################################################
# Generating the subsampled bam using the tool rasusa
######################################################

rasusa aln \
    -s ${seed} \
    -c ${cov} ${aligned_modbam}.bam | \
    samtools sort \
        -@ 10 -O BAM \
        --write-index \
        -o ${sorted_subsampled}.bam #The 100x aligned modbam was provided as input for subsampling

##################################################################
# Generating the modbed from the subsampled bam file using modkit
##################################################################

modkit \
    pileup -t 30 \
    --ref ${reference}.fasta \
    --only-tabs \
    --ignore h ${sorted_subsampled}.bam \
    ${modbed_subsampled}.bed


###############################
# Reordering the modkit output
###############################

# header.tsv has the columns: chrom, p1, p2, mod, coverage, strand, M, UM, other_mod, per
# where M refers to methylated reads, UM refers to unmethylated reads and per refers to the percent methylation

cat ${header}.tsv \
    <(awk 'OFS="\t" {print $1,$2,$3,$4,$10,$6,$12,$13,$14,$11}' ${modbed_subsampled}.bed) \
    > ${modified_modbed}.bed

# Done




