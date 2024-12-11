# This is a scaffold script that contains the steps used to perform the subsampling of dorado 6mA output based on minimum read quality and calculate performance metrics.
# This script is not intended to be used as is.
# The variables need to be adapted to the specific sample name and usecase.


###############################################################
# Filtering the input fastq file based on minimum read quality
###############################################################

nanoq \
    -i ${mod_fastq}.fastq \
    -q ${minQual} -vvv -t 5 \
    -o ${mod_fastq}_${minQual}.fastq \
    -O g -r ${fastq_stats}.txt


#################################
# Aligning the subset fastq file
#################################

minimap2 \
    -x map-ont -a -t 25 -y 
    ${reference}.fasta ${mod_fastq}_${minQual}.fastq \
    | samtools sort -@ 10 -O BAM \
    --write-index \
    -o ${aligned_bam}_${minQual}.bam 


#############################
# Generating the modbed file
#############################

modkit pileup \
    -t 20 \
    --ref ${reference}.fasta \
    --only-tabs \
    --ignore h \
    ${aligned_bam}_${minQual}.bam \
    ${modbed_minQual}_6mA.bed

############################
# Reordering the modbed file
############################

awk 'OFS="\t" {if ($4=="a") print $1,$2,$3,$4,$10,$6,$12,$13,$14,$11}' ${modbed_minQual}_6mA.bed > ${modified_minQual}_6mA.bed


######################################################################
# Intersecting the 6mA modbed with the modbed of the negative control
######################################################################

bedtools intersect \
    -wb -s \
    -a ${modified_minQual}_6mA.bed \
    -b ${DM_minQual}_6mA.bed \
    > ${intersect_DM}_6mA.bed		# Intersecting the double mutant (DM) modbed subsampled based on the minimum read quality


#####################################
# Reordering the intersected bed file
#####################################

awk 'OFS="\t" {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$15,$17,$18,$19,$20}' ${intersect_DM}_6mA.bed > df_${intersect_DM}_6mA.bed


##################################
# Adding the sequence information
##################################

bedtools slop \
    -i df_${intersect_DM}_6mA.bed \
    -s -l 5 -r 11 \
    -g ${reference}.genome \
    > ${slop_coordinates}.bed		# The bed coordinates were extended by 5 nucleotides on the left and 11 nucleotides on the right

bedtools getfasta \
    -fi ${reference}.fasta \
    -bed ${slop_coordinates}.bed -tab -s \
    | cut -f 2 | paste df_${intersect_DM}_6mA.bed - \
    > df_${intersect_DM}_6mA_fasta.bed

##############################################################
# Calculating the F1 score, precision, recall and specificity
##############################################################


# Profiling all GmATC loci. $15 refers to the wildtype percent methylation and $10 refers to the percent methylation of the double mutant negative control.

awk '{if ($16~/^....GATC/) print $0}' df_${intersect_DM}_6mA_fasta.bed | awk '{if ($15<=10) {TN+=1} if ($15>10) {FP+=1} if ($10>=90) {TP+=1} if ($10<90) {FN+=1}} END{recall=(TP/(TP+FN)); precision=(TP/(TP+FP)); F1=(2*((precision*recall)/(precision+recall))); specificity=(TN/(TN+FP)); print "TN="TN,"\t","FP="FP,"\t","TP="TP,"\t","FN="FN,"\t","precision="precision,"\t","recall="recall,"\t","F1="F1,"\t","specificity="specificity}'> metric_${intersect_DM}_6mA.txt

rm ${slop_coordinates}.bed

# Done
