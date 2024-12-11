# This is a scaffold script that contains the steps used to extract human TSS sites.
# This script is not intended to be used as is.
# The variables need to be adapted to the specific sample name and usecase.


###################################################
# Extracting TSS information from the gencode file
###################################################

# Extracting the TSS sites and the coordinates of regions flanking 5kb upstream and downstream of the TSS site

awk -F ";|\t" 'OFS="\t" {if ($3=="gene" && $7=="+") print $1,$4-5000,$4+5000,"TSS_flank","0",$7,$9; if ($3=="gene" && $7=="-") print $1,$5-5000,$5+5000,"TSS_flank","0",$7,$9; if ($3=="gene" && $7=="+") print $1,$4,$4+1,"TSS","0",$7,$9; if ($3=="gene" && $7=="-") print $1,$5-1,$5,"TSS","0",$7,$9}' ${gencode_annotation}.gtf |sed 's/gene_id "//g'|sed 's/"//g' > ${TSS_flank_5kb}.bed

awk 'OFS="\t" {if ($4=="TSS") print $0}' ${TSS_flank_5kb}.bed > ${TSS_coordinates}.bed


#############################################################
# Intersecting the CpG combined dataframe with the TSS sites
############################################################

bedtools intersect -wb -s -a ${TSS_flank_5kb}.bed -b ${CpG_modbed}.bed > ${CpG_modbed}_intersect_TSS.bed


################################
# Reordering the intersected bed
################################

awk 'OFS="\t" {print $1,$2,$3,$4,$5,$6,$7,$12,$14,$15,$16,$17,$18,$19,$20}' ${CpG_modbed}_intersect_TSS.bed > df_${CpG_modbed}_intersect_TSS.bed


############################################################
# Indexing the flanked regions based on the TSS coordinates
############################################################

awk 'NR==FNR {col7[$7]=$7; col2[$7]=$2; next}; ($7==col7[$7]){print $0=$0,"\t",$2-col2[$7]}' ${TSS_coordinates}.bed df_${CpG_modbed}_intersect_TSS.bed > df_indexed_${CpG_modbed}_intersect_TSS.bed

# These are the columns that were retained: 'chr', 'p1', 'p2', 'feature', '.', 'strand', 'ID', 'mod', 'per_bis', '4kHz_sup_v4', '5kHz_hac_v4', '5kHz_sup_v4', '5kHz_hac_v5', '5kHz_sup_v5', 'context', 'index'
# This output is used to plot the TSS plot

# Done


