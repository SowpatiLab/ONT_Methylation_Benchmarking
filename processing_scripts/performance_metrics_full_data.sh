# This is a scaffold script that contains the steps used to calculate performance metrics for 5mC and 6mA data.
# This script is not intended to be used as is.
# The variables need to be adapted to the specific sample name and usecase.


################################################################
# Calculating the performance metrics for 5mC base modification
################################################################

# In the following calculations, $14 is the column containing sequence information. $20 is the column for bisulfite percent methylation of the double mutant (DM) negative control. $13 is the number of unmenthylated reads of DM data. $12 is the methylated reads of DM data. $16 is the bisulfite percent methylation of the 5mC data. $8 is the methylated reads and $9 is the unmethylated reads of the E.coli 5mC data.

############################
# Calculating for mCG motif
############################

# The input bed file is generated by intersecting the E.coli 5mCG modbed with the negative control (DM), and then intersecting the resultant bed with the bisulfite reference.

awk '{if ($14~/^.....CG/) print $0}' ${modbed_mCG}.bed|  awk '{if ($20==0) {TN+=$13; FP+=$12}; if ($16==100) {TP+=$8; FN+=$9}} END{recall=(TP/(TP+FN)); precision=(TP/(TP+FP)); F1=(2*((precision*recall)/(precision+recall))); specificity=(TN/(TN+FP)); print "TN="TN,"\t","FP="FP,"\t","TP="TP,"\t","FN="FN,"\t","precision="precision,"\t","recall="recall,"\t","F1="F1,"\t","specificity="specificity}' > metric_${modbed_mCG}.txt


###############################
# Calculating for CmCWGG motif
###############################

# The input bed file is generated by intersecting the E.coli CmCWGG modbed with the negative control (DM), and then intersecting the resultant with the bisulfite reference.

awk '{if ($14~/^....CCAGG/ || $14~/^....CCTGG/) print $0}' ${modbed_CmCWGG}.bed| awk '{if ($20==0) {TN+=$13; FP+=$12}; if ($16>=95) {TP+=$8; FN+=$9}} END{recall=(TP/(TP+FN)); precision=(TP/(TP+FP)); F1=(2*((precision*recall)/(precision+recall))); specificity=(TN/(TN+FP)); print "TN="TN,"\t","FP="FP,"\t","TP="TP,"\t","FN="FN,"\t","precision="precision,"\t","recall="recall,"\t","F1="F1,"\t","specificity="specificity}' > metric_${modbed_CmCWGG}.txt


############################
# Calculating for GmC motif
############################

# The input bed file is generated by intersecting the E.coli GmC modbed with the negative control (DM), and then intersecting the resultant with the bisulfite reference.

awk '{if ($14~/^....GC/) print $0}' ${modbed_GmC}.bed |awk '{if ($20==0) {TN+=$13; FP+=$12}; if ($16>=95) {TP+=$8; FN+=$9}} END{recall=(TP/(TP+FN)); precision=(TP/(TP+FP)); F1=(2*((precision*recall)/(precision+recall))); specificity=(TN/(TN+FP)); print "TN="TN,"\t","FP="FP,"\t","TP="TP,"\t","FN="FN,"\t","precision="precision,"\t","recall="recall,"\t","F1="F1,"\t","specificity="specificity}' > metric_${modbed_GmC}.txt



################################################################
# Calculating the performance metrics for 6mA base modification
################################################################

# The input bed files are generated by intersecting the 6mA bed file with the negative control. They are then filtered for the motif of interest based on the sequence information.
# In the absence of a reference bisulfite file, the assumption is that all locations with the specific motif must be fully methylated in wildtype data.

# In case of E.coli data, $5 and $10 are the columns for covered depth per location of DM negative control and wildtype data, respectively. $14 is the column containing sequence information. $11 is the column for the percent methylation of the double mutant (DM) negative control. $7 is the column for the percent methylation for E.coli 6mA data.

# In case of H.pylori data, $7 and $13 are the columns for covered depth per location of native and whole genome amplified (WGA) negative control, resepectively. $12 is the column containing the sequence information. $14 is the column for percent methylation of WGA data and the column $8 has percent methylation information of H.pylori 6mA data.

#####################################
# Calculating for E.coli GmATC motif
#####################################

awk '{if ($5>=20 && $10>=20 && $14~/^....GATC/) print $0}' ${modbed_Ecoli_GmATC}.bed |awk '{if ($11<=15) {TN+=1} if ($11>15) {FP+=1} if ($7>=85) {TP+=1} if ($7<85) {FN+=1}} END{recall=(TP/(TP+FN)); precision=(TP/(TP+FP)); F1=(2*((precision*recall)/(precision+recall))); specificity=(TN/(TN+FP)); print "TN="TN,"\t","FP="FP,"\t","TP="TP,"\t","FN="FN,"\t","precision="precision,"\t","recall="recall,"\t","F1="F1,"\t","specificity="specificity}' > metric_${modbed_Ecoli_GmATC}_Ecoli.txt


##############################################
# Calculating for E.coli AmACNNNNNNGTGC motif
##############################################

awk '{if ($5>=20 && $10>=20 && $14~/^....AAC......GTGC/) print $0}' ${modbed_AmACNNNNNNGTGC}.bed |awk '{if ($11<=15) {TN+=1} if ($11>15) {FP+=1} if ($7>=85) {TP+=1} if ($7<85) {FN+=1}} END{recall=(TP/(TP+FN)); precision=(TP/(TP+FP)); F1=(2*((precision*recall)/(precision+recall))); specificity=(TN/(TN+FP)); print "TN="TN,"\t","FP="FP,"\t","TP="TP,"\t","FN="FN,"\t","precision="precision,"\t","recall="recall,"\t","F1="F1,"\t","specificity="specificity}' > metric_${modbed_AmACNNNNNNGTGC}.txt


#######################################
# Calculating for H.pylori GmATC motif
#######################################

awk '{if ($7>=20 && $13>=20 && $12~/^....GATC/) print $0}' ${modbed_Hpylori_GmATC}.bed | awk '{if ($14<=15) {TN+=1} if ($14>15) {FP+=1} if ($8>=85) {TP+=1} if ($8<85) {FN+=1}} END{recall=(TP/(TP+FN)); precision=(TP/(TP+FP)); F1=(2*((precision*recall)/(precision+recall))); specificity=(TN/(TN+FP)); print "TN="TN,"\t","FP="FP,"\t","TP="TP,"\t","FN="FN,"\t","precision="precision,"\t","recall="recall,"\t","F1="F1,"\t","specificity="specificity}'> metric_${modbed_Hpylori_GmATC}.txt


########################################
# Calculating for H.pylori GmANTC motif
########################################

awk '{if ($7>=20 && $13>=20 && $12~/^....GA.TC/) print $0}' ${modbed_GmANTC}.bed | awk '{if ($14<=15) {TN+=1} if ($14>15) {FP+=1} if ($8>=85) {TP+=1} if ($8<85) {FN+=1}} END{recall=(TP/(TP+FN)); precision=(TP/(TP+FP)); F1=(2*((precision*recall)/(precision+recall))); specificity=(TN/(TN+FP)); print "TN="TN,"\t","FP="FP,"\t","TP="TP,"\t","FN="FN,"\t","precision="precision,"\t","recall="recall,"\t","F1="F1,"\t","specificity="specificity}'> metric_${modbed_GmANTC}.txt


#########################################
# Calculating for H.pylori ATTAmAT motif
#########################################

awk '{if ($7>=20 && $13>=20 && $12~/^.ATTAAT/) print $0}' ${modbed_ATTAmAT}.bed | awk '{if ($14<=15) {TN+=1} if ($14>15) {FP+=1} if ($8>=85) {TP+=1} if ($8<85) {FN+=1}} END{recall=(TP/(TP+FN)); precision=(TP/(TP+FP)); F1=(2*((precision*recall)/(precision+recall))); specificity=(TN/(TN+FP)); print "TN="TN,"\t","FP="FP,"\t","TP="TP,"\t","FN="FN,"\t","precision="precision,"\t","recall="recall,"\t","F1="F1,"\t","specificity="specificity}'> metric_${modbed_ATTAmAT}.txt


#######################################
# Calculating for H.pylori CmATG motif
#######################################

awk '{if ($7>=20 && $13>=20 && $12~/^....CATG/) print $0}' ${modbed_CmATG}.bed | awk '{if ($14<=15) {TN+=1} if ($14>15) {FP+=1} if ($8>=85) {TP+=1} if ($8<85) {FN+=1}} END{recall=(TP/(TP+FN)); precision=(TP/(TP+FP)); F1=(2*((precision*recall)/(precision+recall))); specificity=(TN/(TN+FP)); print "TN="TN,"\t","FP="FP,"\t","TP="TP,"\t","FN="FN,"\t","precision="precision,"\t","recall="recall,"\t","F1="F1,"\t","specificity="specificity}'> metric_${modbed_CmATG}.txt


#######################################
# Calculating for H.pylori TCGmA motif
#######################################

awk '{if ($7>=20 && $13>=20 && $12~/^..TCGA/) print $0}' ${modbed_TCGmA}.bed | awk '{if ($14<=15) {TN+=1} if ($14>15) {FP+=1} if ($8>=85) {TP+=1} if ($8<85) {FN+=1}} END{recall=(TP/(TP+FN)); precision=(TP/(TP+FP)); F1=(2*((precision*recall)/(precision+recall))); specificity=(TN/(TN+FP)); print "TN="TN,"\t","FP="FP,"\t","TP="TP,"\t","FN="FN,"\t","precision="precision,"\t","recall="recall,"\t","F1="F1,"\t","specificity="specificity}'> metric_${modbed_TCGmA}.txt


#################################################
# Calculating for H.pylori CTmANNNNNNNNTGT motif
#################################################

awk '{if ($7>=20 && $13>=20 && $12~/^...CTA........TGT/) print $0}' ${modbed_CTmANNNNNNNNTGT}.bed | awk '{if ($14<=15) {TN+=1} if ($14>15) {FP+=1} if ($8>=85) {TP+=1} if ($8<85) {FN+=1}} END{recall=(TP/(TP+FN)); precision=(TP/(TP+FP)); F1=(2*((precision*recall)/(precision+recall))); specificity=(TN/(TN+FP)); print "TN="TN,"\t","FP="FP,"\t","TP="TP,"\t","FN="FN,"\t","precision="precision,"\t","recall="recall,"\t","F1="F1,"\t","specificity="specificity}'> metric_${modbed_CTmANNNNNNNNTGT}.txt


########################################
# Calculating for H.pylori GAAGmA motif
########################################

awk '{if ($7>=20 && $13>=20 && $12~/^.GAAGA/) print $0}' ${modbed_GAAGmA}.bed | awk '{if ($14<=15) {TN+=1} if ($14>15) {FP+=1} if ($8>=85) {TP+=1} if ($8<85) {FN+=1}} END{recall=(TP/(TP+FN)); precision=(TP/(TP+FP)); F1=(2*((precision*recall)/(precision+recall))); specificity=(TN/(TN+FP)); print "TN="TN,"\t","FP="FP,"\t","TP="TP,"\t","FN="FN,"\t","precision="precision,"\t","recall="recall,"\t","F1="F1,"\t","specificity="specificity}'> metric_${modbed_GAAGmA}.txt


#######################################
# Calculating for H.pylori GCmAG motif
#######################################

awk '{if ($7>=20 && $13>=20 && $12~/^...GCAG/) print $0}' ${modbed_GCmAG}.bed | awk '{if ($14<=15) {TN+=1} if ($14>15) {FP+=1} if ($8>=85) {TP+=1} if ($8<85) {FN+=1}} END{recall=(TP/(TP+FN)); precision=(TP/(TP+FP)); F1=(2*((precision*recall)/(precision+recall))); specificity=(TN/(TN+FP)); print "TN="TN,"\t","FP="FP,"\t","TP="TP,"\t","FN="FN,"\t","precision="precision,"\t","recall="recall,"\t","F1="F1,"\t","specificity="specificity}'> metric_${modbed_GCmAG}.txt


# Done


