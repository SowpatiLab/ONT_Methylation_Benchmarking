import polars as pl


def load_dataframe(file, expnam):
    return pl.read_csv(file, separator='\t', has_header=True).with_columns(pl.lit(expnam).alias('exp'))#.filter(pl.col('bis_per')>=1)

def bin_heatmap(df, xlab, ylab, bins=50):
    bins = int(bins)
    binsize = 100/bins

    fdf = (
        df
            .with_columns([
                ((pl.col(xlab)/binsize).floor()*binsize).alias('x_bin'),
                ((pl.col(ylab)/binsize).floor()*binsize).alias('y_bin'),
                pl.lit(1).alias('count')
            ])
            .with_columns([
                pl.col("x_bin").cast(pl.Int64),
                pl.col("y_bin").cast(pl.Int64),
            ])
            .select(['x_bin', 'y_bin', 'exp', 'count'])
            .group_by(['x_bin', 'y_bin', 'exp']).agg(pl.col('count').sum())
    )
    return fdf

d1 = load_dataframe('HG002_sup_4kHz_v4_5mC.bed',   '4kHz_Dorado_v4')
d2 = load_dataframe('HG002_sup_5kHz_v4_5mC.bed',   '5kHz_Dorado_v4')
d3 = load_dataframe('HG002_sup_5kHz_v5_5mC.bed',   '5kHz_Dorado_v5')
d4 = load_dataframe('HG002_deepmod2_4kHz_5mC.bed', '4kHz_DeepMod2').with_columns(pl.col('per')*100)
d5 = load_dataframe('HG002_deepmod2_5kHz_5mC.bed', '5kHz_DeepMod2').with_columns(pl.col('per')*100)

full_dataset = pl.concat([d1,d2,d3,d4,d5])

######################################################################
# non CpG
######################################################################
d = full_dataset.filter(pl.col('bis_context')!='CpG').select(['per', 'bis_per', 'exp'])
binned = bin_heatmap(d, 'bis_per', 'per')
corr = d.group_by(['exp']).agg([
    pl.corr('per', 'bis_per').alias('r'),
    pl.count('per').alias('n')
])

binned.write_csv('full_dataset/consolidated_binned/nonCpG_binned.tsv', separator='\t', include_header=True)
corr.write_csv('full_dataset/consolidated_binned/nonCpG_correl.tsv', separator='\t', include_header=True)


######################################################################
# CHH
######################################################################
d = full_dataset.filter(pl.col('bis_context')=='CHH').select(['per', 'bis_per', 'exp'])
binned = bin_heatmap(d, 'bis_per', 'per')
corr = d.group_by(['exp']).agg([
    pl.corr('per', 'bis_per').alias('r'),
    pl.count('per').alias('n')
])

binned.write_csv('full_dataset/consolidated_binned/CHH_binned.tsv', separator='\t', include_header=True)
corr.write_csv('full_dataset/consolidated_binned/CHH_correl.tsv', separator='\t', include_header=True)

######################################################################
# CHG
######################################################################
d = full_dataset.filter(pl.col('bis_context')=='CHH').select(['per', 'bis_per', 'exp'])
binned = bin_heatmap(d, 'bis_per', 'per')
corr = d.group_by(['exp']).agg([
    pl.corr('per', 'bis_per').alias('r'),
    pl.count('per').alias('n')
])

binned.write_csv('full_dataset/consolidated_binned/CHH_binned.tsv', separator='\t', include_header=True)
corr.write_csv('full_dataset/consolidated_binned/CHH_correl.tsv', separator='\t', include_header=True)
