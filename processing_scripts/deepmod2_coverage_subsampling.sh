# This is a scaffold script that contains the steps used to perform the subsampling of deepmod2 output based on coverage.
# This script is not intended to be used as is.
# The variables need to be adapted to the specific sample name and usecase.


##########################################################
# Extracting unique read names from the out.per_read file
##########################################################

cut -f1 ${out_per_read} | sort -u > ${unique_reads}.txt	# Extracting unique read names from the out.per_read file (100x coverage)

read_count=$(grep -v "read" ${unique_reads} | wc -l)	# Counting the number of unique reads


##########################################################
# Calculating the number of reads needed per coverage
##########################################################

cov=50		# Setting the coverage level (e.g., '10', '20') to be subsampled to

cov_frac=$(echo "scale=0; (${read_count} * ${cov}) / 100" | bc)	
cov_frac=$(printf "%.0f" ${cov_frac})  			


#####################################################################
# Generating a random subset of reads based on the coverage fraction
#####################################################################

grep -v "read_name" ${unique_reads} | shuf | head -n ${cov_frac} > ${cov}x_uniq.per_read	

####################################################################
# Generating the subset output.per_read file based on the coverage
####################################################################

awk 'FNR==NR{ arr[$1]; next }$1 in arr' ${cov}x_uniq.per_read ${out_per_read} > ${cov}x_output.per_read

####################################################
# Using deepmod2 merge to obtain the per_site file
####################################################

prefix="deepmod2_${cov}x"
outdir="deepmod2_merge_out"

# Running deepmod2 merge
conda activate deepmod2

python ./DeepMod2/deepmod2 merge \
    --prefix ${prefix} \
    --output ${outdir} \
    --input ${cov}x_output.per_read \
    --cpg_out

# Done
