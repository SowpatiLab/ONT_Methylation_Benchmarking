# This is a scaffold script that contains the steps used to perform the subsampling of deepmod2 output based on minimum read quality and calculating the performance metrics.
# This script is not intended to be used as is.
# The variables need to be adapted to the specific sample name and usecase.


##################################################################################################
# Running deepmod2 merge to obtain per_site information for the given minimum read quality cutoff
##################################################################################################

prefix="deepmod2_${minQual}"
outdir="deepmod2_merge_minQual"

conda activate deepmod2

python ./DeepMod2/deepmod2 merge \
    --prefix ${prefix} \
    --output ${outdir} \
    --qscore_cutoff ${minQual} \
    --input ${output_per_read} \
    --cpg_out		


###############################################################
# Intersecting the MSssI-treated data and the negative control
###############################################################

bedtools intersect -wb -s \
    -a <(awk 'OFS="\t" {if ($5=="True") print $1,$2,$3,$7,$8,$4,$6,$9}' ${MSssI_minQual_per_site}) \
    -b <(awk 'OFS="\t" {if ($5=="True") print $1,$2,$3,$7,$8,$4,$6,$9}' ${DM_minQual_per_site}) \
    > ${intersect_MSssI_DM}.bed


#################################################################
# Intersecting the output bed with the reference bisulfite file 
#################################################################

bedtools intersect -wb -s \
    -a ${intersect_MSssI_DM}.bed \
    -b ${bisulfite}.bed \
    > ${intersect_bis}.bed	# The bisulfite reference has methylation information for both MSssI and DM data


##########################################################################
# Reordering the intersected bed file for performance metrics calculations
##########################################################################

awk 'OFS="\t" {print $1,$2,$3,$4,$5,$6,$7,$8*100,$12,$13,$14,$15,$16*100,$39,$41,$42,$52,$21,$23,$24,$34,$50}' \
    ${intersect_bis}.bed \
    > df_${intersect_bis}.bed


##############################################################
# Calculating the F1 score, precision, recall and specificity
##############################################################

# $21 refers to the percent bisulfite methylation of the double mutant (DM) negative control. $17 refers to the percent bisulfite methylation of the MSssI-treated data.
# $10 refers to the number of unmethylated reads and $9 is the number of methylated reads in DM. $4 is the number of methylated reads and $5 is the number of unmethylated reads of MSssI data. 
awk '{if ($21==0) {TN+=$10; FP+=$9}; if ($17==100) {TP+=$4; FN+=$5}} END{recall=(TP/(TP+FN)); precision=(TP/(TP+FP)); F1=(2*((precision*recall)/(precision+recall))); specificity=(TN/(TN+FP)); print "TN="TN,"\t","FP="FP,"\t","TP="TP,"\t","FN="FN,"\t","precision="precision,"\t","recall="recall,"\t","F1="F1,"\t","specificity="specificity}' df_${intersect_bis}.bed > metrics_${intersect_bis}.txt

# Done
