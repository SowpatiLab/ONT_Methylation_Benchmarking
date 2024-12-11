# This is a scaffold script that contains the steps used to perform the subsampling of dorado 5mC output based on minimum read quality and calculate performance metrics.
# This script is not intended to be used as is.
# The variables need to be adapted to the specific sample name and usecase.


###############################################################
# Filtering the input fastq file based on minimum read quality
###############################################################

nanoq \
    -i ${mod_fastq}.fastq \
    -q ${minQual} \
    -vvv -t 5 \
    -o ${mod_fastq}_${minQual}.fastq \
    -O g \
    -r ${fastq_stats}.txt


#################################
# Aligning the subset fastq file
#################################

minimap2 \
    -x map-ont -a -t 25 -y \
    ${reference}.fasta ${mod_fastq}_${minQual}.fastq | \
    samtools sort \
        --write-index \
        -@ 10 -O BAM -o ${aligned_bam}_${minQual}.bam

#############################
# Generating the modbed file
#############################

modkit pileup \
    -t 20 --only-tabs --ref ${reference}.fasta \
    --ignore h ${aligned_bam}_${minQual}.bam \
    ${modbed_minQual}_5mC.bed


############################
# Reordering the modbed file
############################

awk 'OFS="\t" {if ($4=="m") print $1,$2,$3,$4,$10,$6,$12,$13,$14,$11}' ${modbed_minQual}_5mC.bed > ${modified_minQual}_5mC.bed


#####################################################################
# Intersecting the 5mC modbed with the modbed of the negative control
#####################################################################

bedtools intersect \
    -wb -s \
    -a ${modified_minQual}_5mC.bed \
    -b ${DM_minQual}_5mC.bed \
    > ${intersect_DM}_5mC.bed		# Intersecting the double mutant (DM) modbed subsampled based on the minimum read quality


###############################################
# Intersecting with the bisulfite ground truth
###############################################

bedtools intersect \
    -wb -s \
    -a ${intersect_DM}_5mC.bed \
    -b ${bisulfite}.bed \
    > ${intersect_bis}_5mC.bed		# The bisulfite bed has both 5mC WT/MSssI and DM negative control data

cat ${header}.tsv <(awk 'OFS="\t" {print $1,$2,$3,$38,$56,$6,$7,$8,$4,$17,$18}' ${intersect_bis}_5mC.bed) > df_${intersect_bis}_5mC.bed

# The example header.tsv for WT data has the following columns (M is methylated reads and UM is unmethylated reads. per is percent methylation):
# chrom, p1, p2, bis_per_DM, bis_per_WT, strand, M_WT, UM_WT, mod, M_DM, UM_DM


##############################################################
# Calculating the F1 score, precision, recall and specificity
##############################################################

# $4 refers to the percent bisulfite methylation of the double mutant (DM) negative control. $5 refers to the percent bisulfite methylation of the wildtype/MSssI-treated data.
# $11 refers to the number of unmethylated reads and $10 is the number of methylated reads in DM. $7 is the number of methylated reads and $8 is the number of unmethylated reads of wildtype/MSssI-treated data.

awk '{if ($4==0) {TN+=$11; FP+=$10}; if ($5==100) {TP+=$7; FN+=$8}} END{recall=(TP/(TP+FN)); precision=(TP/(TP+FP)); F1=(2*((precision*recall)/(precision+recall))); specificity=(TN/(TN+FP)); print "TN="TN,"\t","FP="FP,"\t","TP="TP,"\t","FN="FN,"\t","precision="precision,"\t","recall="recall,"\t","F1="F1,"\t","specificity="specificity}' df_${intersect_bis}_5mC.bed > metric_${intersect_bis}_5mC.txt

# Done
