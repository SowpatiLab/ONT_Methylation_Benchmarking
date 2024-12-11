[ ! -d combined_filt ] && mkdir combined_filt
for exp in { HG002_sup_4kHz_v4_5mC,HG002_sup_5kHz_v4_5mC,HG002_sup_5kHz_v5_5mC,HG002_deepmod2_4kHz_5mC,HG002_deepmod2_5kHz_5mC};
do
    # preserve only the percentage methylation values along with the context name
    head -n 1 HG002_sup_4kHz_v4_5mC.bed | cut -f 11,12,15 > combined_filt/${exp}.bed
    sed '1d' ${exp}.bed | awk -v OFS='\t' '{if($13+$14>=20 && $5>=20) print $11,$12,$15}' >> combined_filt/${exp}.bed
    echo "done $exp"
done

[ ! -d combined_filt ] && mkdir combined_filt;


# This script does the following:
#    1. simplifty datasets for quicker loading
#    2. Aggregates the data points to better suit heatmap requirements
#    3. Bins the count values into bins of 50
#    4. pre-calculates correlation values for each experimental set
python consolidate_datasets.py

[ ! -d plots ] && mkdir plots 
# plot nonCpG context heatmaps
Rscript ../human_plot_non_CpG_contexts.R
