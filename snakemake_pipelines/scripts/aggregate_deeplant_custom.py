import polars as pl
import argparse

def load_df(fi):
    return pl.read_csv(fi, 
        separator='\t', columns=[3,4,6],
        new_columns=['chrom', 'p1','methylation_rate'],
        schema_overrides={ 'chrom': pl.Categorical, 'p1': pl.UInt32, 'methylation_rate': pl.Float32}
    ).lazy()

def aggregate(df):
    return (
        df.group_by(['chrom', 'p1'])
        .agg([
            pl.col('methylation_rate').len().alias('coverage'),
            pl.col('methylation_rate').filter(pl.col('methylation_rate')>0.1).len().alias('M1'),
            pl.col('methylation_rate').filter(pl.col('methylation_rate')>0.2).len().alias('M2'),
            pl.col('methylation_rate').filter(pl.col('methylation_rate')>0.3).len().alias('M3'),
            pl.col('methylation_rate').filter(pl.col('methylation_rate')>0.4).len().alias('M4'),
            pl.col('methylation_rate').filter(pl.col('methylation_rate')>0.5).len().alias('M5'),
            pl.col('methylation_rate').filter(pl.col('methylation_rate')>0.6).len().alias('M6'),
            pl.col('methylation_rate').filter(pl.col('methylation_rate')>0.7).len().alias('M7'),
            pl.col('methylation_rate').filter(pl.col('methylation_rate')>0.8).len().alias('M8'),
            pl.col('methylation_rate').filter(pl.col('methylation_rate')>0.9).len().alias('M9'),
        ])
    ).collect(engine='streaming')

def get_args():
    parser = argparse.ArgumentParser(description=(f"""deepplant aggregation for aggregating at various thresholds
                                                  
    output fmt:                                                
    chrom	p1	  coverage	M1	M2	M3	M4	M5	M6	M7	M8	M9
    chr1	48527163	23	8	0	0	0	0	0	0	0	0
    chr1	39347824	28	19	2	1	1	1	1	1	1	0
    chr1	93552966	26	3	1	0	0	0	0	0	0	0
    chr1	89991010	25	8	4	3	2	2	1	1	1	1
    chr16	58520235	1	1	1	0	0	0	0	0	0	0
                                                  
    columns explaines: 
    M1 => M for threshold = 0.1
    M2 => M for threshold = 0.2
    M3 => M for threshold = 0.3
    M4 => M for threshold = 0.4
    M5 => M for threshold = 0.5
    M6 => M for threshold = 0.6
    M7 => M for threshold = 0.7
    M8 => M for threshold = 0.8
    M9 => M for threshold = 0.9
    """
),formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-i', '--input-readwise', type=str, help='deepplant readwise file', required=True)
    parser.add_argument('-o', '--output', type=str, help='output file to write to', required=True)
    
    return parser.parse_args()

if __name__=="__main__":
    args = get_args()
    input_file, out_file = args.input_readwise, args.output

    e = load_df(input_file) 
    e_aggregated = aggregate(e)
    e_aggregated.write_csv(out_file, separator='\t', include_header=True)