import polars as pl
import argparse

def get_args():
    parser = argparse.ArgumentParser(description=(f"""This script is used to filter and reorganise bed files to a manageable format
    """
),formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-i',  '--input',       type=str, help='Input modbed file with fatsa', required=True)
    parser.add_argument('-o',  '--output',      type=str, help='Output File',                  required=True)
    parser.add_argument('-c',  '--context',     type=str, help='context',                      required=True)
    parser.add_argument('-sp', '--species',     type=str, help='context',                      required=True)
    parser.add_argument('-cd', '--condition',   type=str, help='condition',                    required=True)
    parser.add_argument('-m',  '--model',       type=str, help='model',                        required=True)
    parser.add_argument('-t',  '--tool',        type=str, help='tool',                         required=True)
    parser.add_argument('-se', '--sample-rate', type=str, help='sample rate',                  required=True)
    
    return parser.parse_args()

args = get_args()
CONTEXT = args.context
INPUT   = args.output
OUTPUT  = args.output



def load_intersect(fi):
    mod_list = { 'a':"6mA", 'm':"5mC", '21839':"4mC" }

    ## for human specific contexts
    # strMatcher = {
    #     'CpG': r'^CG[ATGC]$',
    #     'CHG': r'^C[CTA]G$',
    #     'CHH': r'^C[CTA][CTA]$',
    # }
    strMatcher = {
        '5mC': r'^\w{5}C\w{11}$',
        '4mC': r'^\w{5}C\w{11}$',
        '6mA': r'^\w{5}A\w{11}$',
    }

    cols = [ 'chrom', 'p1', 'p2', 'mod', 's', 'strand', 'coverage', 'per', 'M', 'UM', 'per_bis', 'M_bis', 'UM_bis', 'cotext' ]
    cols_vals = 'chrom,p1,p2,mod,coverage,strand,M,UM,tool,model,sample,sample_rate,full_context,species,per'.split(',')

    dataframe = pl.read_csv(fi, separator='\t', has_header=False, new_columns=cols,
        schema={
            'chrom': pl.String,
            'p1': pl.Int64,
            'p2': pl.Int64,
            'mod': pl.String,
            's': pl.String,
            'strand': pl.String,
            'coverage': pl.Int64,
            'per': pl.Float64,
            'M':pl.Int64,
            'UM': pl.Int64,
            'per_bis': pl.Float64,
            'M_bis': pl.Int64,
            'UM_bis': pl.Int64,
            'cotext': pl.String
        }
    )

    species     = args.species
    condition   = args.condition
    model       = args.model
    tool        = args.tool
    sample_rate = args.sample_rate

    return (
        dataframe
        .filter(pl.col('full_context').str.to_uppercase().str.contains(strMatcher[CONTEXT]))
        .with_columns([
            pl.col('mod').replace(mod_list),
            pl.lit(sample_rate).alias('sample_rate'),
            pl.lit(species).alias('species'),
            pl.lit(condition).alias('sample'),
            pl.lit(tool).alias('tool'),
            pl.lit(model).alias('model')
        ])
        .select(cols_vals)
    )

df = load_intersect(INPUT)
df.write_csv(OUTPUT, separator='\t', include_header=False)
