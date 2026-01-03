import polars as pl
import argparse, sys

def get_args():
    parser = argparse.ArgumentParser(description=(f"""This script can be used to aggregate rockfish readwise output
    to generate the site-wise methylation output
    """
),formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-i', '--input-readwise', type=str, help='Rockfish readwise data. (predictions.tsv)', required=True)
    parser.add_argument('-u', '--upper-thresh', type=float, help='upper threshold | values above this is considered methylated', required=False, default=0.5)
    parser.add_argument('-l', '--lower-thresh', type=float, help='lower threshold | values below this is considered unmethylated', required=False, default=0.5)
    parser.add_argument('-o', '--output', type=str, help='output file to save aggregatred data to ', required=True)
    
    return parser.parse_args()

def parse_rockfish(rockfish, get_predicted=False, type_suffix=''):
    if get_predicted:
        value_col = pl.when(pl.col('prob') > 0.5).then(1).otherwise(0).alias(
            f'Rockfish{type_suffix}').cast(pl.UInt8)
    else:
        value_col = pl.col('prob').alias(f'prob{type_suffix}').cast(pl.Float32)

    rockfish = rockfish.select([
        pl.col('read_id'),
        pl.col('ctg').cast(pl.Categorical),
        pl.col('pos').cast(pl.UInt32), value_col
    ])

    return rockfish

def consolidate(df, ut=0.5, lt=0.5):
    N_NOMOD_COL = (pl.col('prob') <= lt).sum().alias('UM')
    N_MOD_COL = (pl.col('prob') > ut).sum().alias('M')
    FREQ_COL = (pl.col('M') /
                (pl.col('UM') + pl.col('M'))).alias('freq')
    END_COL = (pl.col('pos') + 1).cast(pl.UInt32).alias('end')

    read_level = parse_rockfish(df)
    read_level = read_level.with_columns((pl.col('prob') - 0.5).abs().alias('abs')) \
                        .sort('abs', descending=True) \
                        .group_by(['ctg', 'pos']) \
                        .head(n=1000) \
                        .group_by(['ctg', 'pos']) \
                        .agg([N_MOD_COL, N_NOMOD_COL]) \
                        .sort(['ctg', 'pos'])#.collect()

    read_level = read_level.select([
        'ctg',
        pl.col('pos').alias('start'), END_COL, FREQ_COL, 'M', 'UM'
    ])
    return read_level

if __name__=="__main__":
    args = get_args()
    pred_file, out_file = args.input_readwise, args.output
    up_thresh, lw_thresh =  args.upper_thresh,  args.lower_thresh

    d = pl.read_csv(pred_file, separator='\t', has_header=True)
    e = consolidate(d, ut=up_thresh, lt=lw_thresh)
    e.write_csv(out_file, separator='\t', include_header=False)
