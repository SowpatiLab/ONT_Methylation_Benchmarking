# This is a scaffold script that contains the steps used to generate the methylation dataframe of cytosines neighbouring human CpG sites.
# This script is not intended to be used as is.
# The variables need to be adapted to the specific sample name and usecase.

# The cytosines flanking 5 nucleotides upstream and downstream of these CpGs will be profiled. 
# In our study, cytosines flanking only fully methylated (100% bisulfite methylation) and fully unmethylated (0% bisulfite methylation) were looked at.


#########################################################
# Generating bed files with C methylation upstream of CpG
#########################################################

# Bedtools intersect intersects locations with >=20 depth in the CHG/CHH bed file with the locations where a C is upstream within five nucleotides of a CG. The methylation percent is determined by $per. The position of C from CG is relatively changed in each iteration. The output is a bed file with all the Cs upstream of the CG and their respective methylation information.

out1_us="${CpG_modbed}_per"$per"_us1.bed"
bedtools intersect \
	-s -a <(grep -v "strand" ${CHG_modbed}.bed | awk '{if ($15>=20 && $5>=20) print $0}') -b <(awk 'IGNORECASE=1 {if ($14 ~/^....CCG/) print $0}' ${CpG_modbed}.bed | awk -v var=$per '{if ($15>=20 && $5>=20 && $16==var) print $0}' | awk 'OFS="\t" {if ($6=="+") {$2=$2-1;$3=$3-1;} else if ($6=="-") {$2=$2+1;$3=$3+1;} print}'| grep -v "strand")> $out1_us

out2_us="${CpG_modbed}_per"$per"_us2.bed"
bedtools intersect \
	-s -a <(grep -v "strand" ${CHG_modbed}.bed | awk '{if ($15>=20 && $5>=20) print $0}') -b <(awk 'IGNORECASE=1 {if ($14 ~/^...C.CG/) print $0}' ${CpG_modbed}.bed | awk -v var=$per '{if ($15>=20 && $5>=20 && $16==var) print $0}' | awk 'OFS="\t" {if ($6=="+") {$2=$2-2;$3=$3-2;} else if ($6=="-") {$2=$2+2;$3=$3+2;} print}'| grep -v "strand")> $out2_us

out3_us="${CpG_modbed}_per"$per"_us3.bed"
bedtools intersect \
	-s -a <(grep -v "strand" ${CHG_modbed}.bed | awk '{if ($15>=20 && $5>=20) print $0}') -b <(awk 'IGNORECASE=1 {if ($14 ~/^..C..CG/) print $0}' ${CpG_modbed}.bed | awk -v var=$per '{if ($15>=20 && $5>=20 && $16==var) print $0}' | awk 'OFS="\t" {if ($6=="+") {$2=$2-3;$3=$3-3;} else if ($6=="-") {$2=$2+3;$3=$3+3;} print}'| grep -v "strand")> $out3_us

out4_us="${CpG_modbed}_per"$per"_us4.bed"
bedtools intersect \
	-s -a <(grep -v "strand" ${CHG_modbed}.bed | awk '{if ($15>=20 && $5>=20) print $0}') -b <(awk 'IGNORECASE=1 {if ($14 ~/^.C...CG/) print $0}' ${CpG_modbed}.bed | awk -v var=$per '{if ($15>=20 && $5>=20 && $16==var) print $0}' | awk 'OFS="\t" {if ($6=="+") {$2=$2-4;$3=$3-4;} else if ($6=="-") {$2=$2+4;$3=$3+4;} print}'| grep -v "strand") | awk 'OFS="\t" {if ($14!~/^.....C.CGCG/) print $0}' > $out4_us

out5_us="${CpG_modbed}_per"$per"_us5.bed"
bedtools intersect \
	-s -a <(grep -v "strand" ${CHG_modbed}.bed | awk '{if ($15>=20 && $5>=20) print $0}') -b <(awk 'IGNORECASE=1 {if ($14 ~/^C....CG/) print $0}' ${CpG_modbed}.bed | awk -v var=$per '{if ($15>=20 && $5>=20 && $16==var) print $0}' | awk 'OFS="\t" {if ($6=="+") {$2=$2-5;$3=$3-5;} else if ($6=="-") {$2=$2+5;$3=$3+5;} print}'| grep -v "strand") | awk 'OFS="\t" {if ($14!~/^.....C..CGCG/) print $0}' > $out5_us


###########################################################
# Generating bed files with C methylation downstream of CpG
###########################################################

# Bedtools intersect intersects locations with >=20 depth in the CHG/CHH bed file with the locations where a C is downstream within five nucleotides of a CG. The methylation percent is determined by $per. The position of C from CG is relatively changed in each iteration. The output is a bed file with all the Cs downstream of the CG and their respective methylation information.

out1_ds="${CpG_modbed}_per"$per"_ds1.bed"
bedtools intersect -s -a <(grep -v "strand" ${CHG_modbed}.bed| awk '{if ($15>=20 && $5>=20) print $0}') -b <(awk 'IGNORECASE=1 {if ($14 ~/^.....CGC/) print $0}' ${CpG_modbed}.bed | awk -v var=$per '{if ($15>=20 && $5>=20 && $16==var) print $0}' | awk 'OFS="\t" {if ($6=="+") {$2=$2+2;$3=$3+2;} else if ($6=="-") {$2=$2-2;$3=$3-2;} print}'| grep -v "strand") > $out1_ds

out2_ds="${CpG_modbed}_per"$per"_ds2.bed"
bedtools intersect -s -a <(grep -v "strand" ${CHG_modbed}.bed| awk '{if ($15>=20 && $5>=20) print $0}') -b <(awk 'IGNORECASE=1 {if ($14 ~/^.....CG.C/) print $0}' ${CpG_modbed}.bed | awk -v var=$per '{if ($15>=20 && $5>=20 && $16==var) print $0}' | awk 'OFS="\t" {if ($6=="+") {$2=$2+3;$3=$3+3;} else if ($6=="-") {$2=$2-3;$3=$3-3;} print}'| grep -v "strand") > $out2_ds

out3_ds="${CpG_modbed}_per"$per"_ds3.bed"
bedtools intersect -s -a <(grep -v "strand" ${CHG_modbed}.bed| awk '{if ($15>=20 && $5>=20) print $0}' ) -b <(awk 'IGNORECASE=1 {if ($14 ~/^.....CG..C/) print $0}' ${CpG_modbed}.bed | awk -v var=$per '{if ($15>=20 && $5>=20 && $16==var) print $0}' | awk 'OFS="\t" {if ($6=="+") {$2=$2+4;$3=$3+4;} else if ($6=="-") {$2=$2-4;$3=$3-4;} print}'| grep -v "strand") | awk 'OFS="\t" {if ($14!~/^.CGCGC/) print $0}' > $out3_ds

out4_ds="${CpG_modbed}_per"$per"_ds4.bed"
bedtools intersect -s -a <(grep -v "strand" ${CHG_modbed}.bed| awk '{if ($15>=20 && $5>=20) print $0}' ) -b <(awk 'IGNORECASE=1 {if ($14 ~/^.....CG...C/) print $0}' ${CpG_modbed}.bed | awk -v var=$per '{if ($15>=20 && $5>=20 && $16==var) print $0}' | awk 'OFS="\t" {if ($6=="+") {$2=$2+5;$3=$3+5;} else if ($6=="-") {$2=$2-5;$3=$3-5;} print}'| grep -v "strand") | awk 'OFS="\t" {if (($14!~/^CGCG.C/) && ($14!~/^CG.CGC/)) print $0}' > $out4_ds

out5_ds="${CpG_modbed}_per"$per"_ds5.bed"
bedtools intersect -s -a <(grep -v "strand" ${CHG_modbed}.bed| awk '{if ($15>=20 && $5>=20) print $0}' ) -b <(awk 'IGNORECASE=1 {if ($14 ~/^.....CG....C/) print $0}' ${CpG_modbed}.bed | awk -v var=$per '{if ($15>=20 && $5>=20 && $16==var) print $0}' | awk 'OFS="\t" {if ($6=="+") {$2=$2+6;$3=$3+6;} else if ($6=="-") {$2=$2-6;$3=$3-6;} print}'| grep -v "strand") | awk 'OFS="\t" {if (($14!~/^GCG..C/) && ($14!~/^G.CG.C/) && ($14!~/^G..CGC/)) print $0}' > $out5_ds


####################################################################################################################
# Finding the number of cytosines upstream and downstream of the CpG site (at each location upto 5 nucleotides away)
####################################################################################################################
 
wc -l *_per${per}_us*.bed

wc -l *_per${per}_ds*.bed


##########################################################################
# Setting the number of common sites to be profiled from each output file.
##########################################################################

num=lowest_num	# This number is the lowest number of loci as per the previous step.


##############################################################
# Choosing a random subset of locations from each output file
##############################################################

for file in *_per${per}_us*.bed; do 
	shuf ${file} | head -n ${num} > shuf_${file}
done

for file in *_per${per}_ds*.bed; do 
	shuf ${file} | head -n ${num} > shuf_${file}
done


#######################
# Making the dataframe
#######################

# These columns were retained: "chr", "p1", "p2", "per_ONT", "per_bis", "strand", "index", "model"
cat <(awk -v name=$model 'OFS="\t" {print $1,$2,$3,$7,$16,$6,"-5",name}' shuf_${out5_us}) <(awk -v name=$model 'OFS="\t" {print $1,$2,$3,$7,$16,$6,"-4",name}' shuf_${out4_us}) <(awk -v name=$model 'OFS="\t" {print $1,$2,$3,$7,$16,$6,"-3",name}' shuf_${out3_us}) <(awk -v name=$model 'OFS="\t" {print $1,$2,$3,$7,$16,$6,"-2",name}' shuf_${out2_us}) <(awk -v name=$model 'OFS="\t" {print $1,$2,$3,$7,$16,$6,"-1",name}' shuf_${out1_us}) > df_indexed_${model}_upstream_CpG_per${per}_shuf.bed

cat <(awk -v name=$model 'OFS="\t" {print $1,$2,$3,$7,$16,$6,"1",name}' shuf_${out1_ds}) <(awk -v name=$model 'OFS="\t" {print $1,$2,$3,$7,$16,$6,"2",name}' shuf_${out2_ds}) <(awk -v name=$model 'OFS="\t" {print $1,$2,$3,$7,$16,$6,"3",name}' shuf_${out3_ds}) <(awk -v name=$model 'OFS="\t" {print $1,$2,$3,$7,$16,$6,"4",name}' shuf_${out4_ds}) <(awk -v name=$model 'OFS="\t" {print $1,$2,$3,$7,$16,$6,"5",name}' shuf_${out5_ds}) > df_indexed_${model}_downstream_CpG_per${per}_shuf.bed

# Done
