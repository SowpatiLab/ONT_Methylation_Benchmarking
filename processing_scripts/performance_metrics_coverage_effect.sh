# This is a scaffold script that contains the steps used to calculate the performance metrics of subsampled modbeds based on coverage.
# This script is not intended to be used as is.
# The variables need to be adapted to the specific sample name and usecase.


# Setting the array with the different coverage values
array=(5 10 20 30 40 50 60 70)


#####################################################################################
# Calculating the performance metrics for the 5mC E.coli data of different coverages
#####################################################################################

for cov in "${array[@]}"; do
	
	#####################################################################
	# Intersecting the WT/MSssI and DM data with the bisulfite reference
	#####################################################################
		
	# WT refers to wildtype and DM refers to double mutant (negative control)
	# The bisulfite data has all locations where WT/MSssI data is >=90% methylated and DM is <=10% methylated
		
	bedtools intersect -wb -s -a ${modbed}_5mC_${cov}.bed -b ${DM}_5mC_${cov}.bed > ${intersect_DM}_5mC_${cov}.bed
	bedtools intersect -wb -s -a ${intersect_DM}_5mC_${cov}.bed -b  ${bisulfite}.bed > ${intersect_bis}_5mC_${cov}.bed

	awk 'OFS="\t" {if ($4=="m" && $38>=90 && $56<=10) print $1,$2,$3,$4,$10,$6,$20,$36,$38,$56}' ${intersect_bis}_5mC_${cov}.bed) > df_${intersect_bis}_5mC_${cov}.bed
		
	# This file contains methylation information for wildtype/MSssI, negative control (DM) and the corresponding bisulfite data for each location.
	# The columns are chr, p1, p2, mod, per_5mC, strand, per_DM, full_context, bis_5mC, bis_DM

		
	##############################################################
	# Calculating the F1 score, precision, recall and specificity
	##############################################################
		
	# Here, $5 is the column for percent methylation of the WT/MSssI data and $7 is the column for percent methylation of the double mutant (DM) negative control. 
	# The number of actual positives ($AP1) refer to the number of locations where percent bisulfite methylation is >=90%
	# The number of actual negatives ($AN1) refer to the number of locations where percent bisulfite methylation of the negative control is <=10%
	awk -v var1=$AP1 -v var2=$AN1 '{if ($7<=10) {TN+=1} if ($7>10) {FP+=1} if ($5>=90) {TP+=1} if ($5<90) {FN+=1}} END{recall=(TP/var1); precision=(TP/(TP+FP)); F1=(2*((precision*recall)/(precision+recall))); specificity=(TN/var2); print "TN="TN,"\t","FP="FP,"\t","TP="TP,"\t","FN="FN,"\t","precision="precision,"\t","recall="recall,"\t","F1="F1,"\t","specificity="specificity}' df_${intersect_bis}_5mC_${cov}.bed > metric_${intersect_bis}_5mC_${cov}.txt
		
done
	
	
	
#################################################################################	
# Calculating the performance metrics for 6mA E.coli data of different coverages
#################################################################################
		
for cov in "${array[@]}"; do
		
	################################################
	# Intersecting the WT modbed with the DM modbed
	################################################
			
	bedtools intersect -wb -s -a ${DM}_6mA_${cov}.bed -b ${modbed_WT}_6mA_${cov}.bed | awk 'OFS="\t" {if ($4=="a") print $1,$2,$3,$4,$5,$6,$7,$8,$10,$16,$18,$19,$21,$22}' > ${intersect_DM}_6mA_${cov}.bed
	bedtools intersect -s -a ${intersect_DM}_6mA_${cov}.bed -b ${Ecoli_reference_6mA}.bed > ${intersect_ref}_6mA_${cov}.bed
			
	# Here, the reference file ($Ecoli_reference_6mA.bed) is the 100x modbed file containing all locations where wildtype data is >=90% methylated and DM data is <=10% methylated 
			
			
	##############################################################
	# Calculating the F1 score, precision, recall and specificity
	##############################################################
		
	# Here, $15 is the column for percent methylation of the DM negative control and $13 is the column for percent methylation of the 6mA wildtype data in the coverage bed being profiled.
	# The number of actual positives ($AP2) refers to the number of locations in the 100x dataset where percent 6mA methylation is >=90%.
	# The number of actual negatives ($AN2) refers to the number of locations in the 100x negative control dataset where percent 6mA methylation is <=10%.
	awk -v var1=$AP2 -v var2=$AN2 '{if ($9<=15) {TN+=1} if ($9>15) {FP+=1} if ($13>=85) {TP+=1} if ($13<85) {FN+=1}} END{recall=(TP/var1); precision=(TP/(TP+FP)); F1=(2*((precision*recall)/(precision+recall))); specificity=(TN/var2); print "TN="TN,"\t","FP="FP,"\t","TP="TP,"\t","FN="FN,"\t","precision="precision,"\t","recall="recall,"\t","F1="F1,"\t","specificity="specificity}' ${intersect_ref}_6mA_${cov}.bed > metric_${intersect_ref}_6mA_${cov}.txt
	
done


###################################################################################	
# Calculating the performance metrics for 6mA H.pylori data of different coverages
###################################################################################
		
for cov in "${array[@]}"; do
		
	#################################################
	# Intersecting the NAT modbed with the WGA modbed
	#################################################
			
	# NAT refers to native and WGA refers to whole genome amplified (negative control)
			
	bedtools intersect -wb -s -a ${modbed_NAT}_6mA_${cov}.bed -b ${modbed_WGA}_6mA_${cov}.bed | awk 'OFS="\t" {if (($4=="a") && ($23=="a")) print $1,$2,$3,$4,".",$6,$10,$11,$12,$13,$14,$19,$29,$30,$31,$32,$33}' > ${intersect_WGA}_6mA_${cov}.bed
	# The columns here are: chr, p1, p2, mod, ., strand, NAT_cov, NAT_per, NAT_M, NAT_UM, NAT_other, full_context, WGA_cov, WGA_per, WGA_M, WGA_UM, WGA_other

	bedtools intersect -s -a ${intersect_WGA}_6mA_${cov}.bed -b ${Hpylori_reference_6mA} > ${intersect_ref}_6mA_${cov}.bed
			
	# The reference file here is the 100x modbed file containing all locations where wildtype data is >=90% methylated and WGA data is <=10% methylated 
			
			
			
	##############################################################
	# Calculating the F1 score, precision, recall and specificity
	##############################################################
			
	# Here, $14 is the column for percent methylation of the DM negative control and $8 is the column for percent methylation of the 6mA wildtype data in the coverage bed being profiled.
	# The number of actual positives refers to the number of locations in the 100x dataset where percent 6mA methylation is >=90%.
	# The number of actual negatives refers to the number of locations in the 100x negative control dataset where percent 6mA methylation is <=10%
			
	awk -v var1=$AP3 -v var2=$AN3 '{if ($14<=15) {TN+=1} if ($14>15) {FP+=1} if ($8>=85) {TP+=1} if ($8<85) {FN+=1}} END{recall=(TP/var1); precision=(TP/(TP+FP)); F1=(2*((precision*recall)/(precision+recall))); specificity=(TN/var2); print "TN="TN,"\t","FP="FP,"\t","TP="TP,"\t","FN="FN,"\t","precision="precision,"\t","recall="recall,"\t","F1="F1,"\t","specificity="specificity}' ${intersect_WGA}_6mA_${cov}.bed > "$metric3" 
done


# Done
