from pathlib import Path
import polars as pl
import re, sys

select  = ['chrom', 'p1', 'p2', 'mod', 'coverage', 'strand', 'M', 'UM', 'per',
                   'flowcell', 'tool', 'model', 'sample', 'acc', 'sample_rate', 'species', 'full_context']

def load_modkit(fi):
    exp, sr, acc, model, mod = re.match(r'(.*)_([4|5])kHz_(sup|hac)_(v[45.2]+r[123])_([564]m[CAG]+)', Path(fi).stem).group(1, 2, 3, 4, 5)
    species = exp.split('_')[0]
    modlookup = {'a': "6mA", 'm': "5mC", 'h': "5hmC", '21839': '4mC'}
    filter_mod = {'5mC':'m','5mCG':'m', '5hmC':'h','5hmCG':'h','4mC':'21839','6mA':'a'}
    strMatcher = {
        '5mC':   r'^\w{5}C\w{11}$',
        '5mCG':  r'^\w{5}C\w{11}$',
        '5hmC':  r'^\w{5}C\w{11}$',
        '5hmCG': r'^\w{5}C\w{11}$',
        '4mC':   r'^\w{5}C\w{11}$',
        '6mA':   r'^\w{5}A\w{11}$',
    }
    return (
        pl.read_csv(fi, separator='\t', has_header=False,
            columns = [0,1,2,3,4,5,10,11,12,18],
            new_columns=['chrom', 'p1','p2','mod','coverage','strand','per','M','UM','full_context']
        )
        .filter(pl.col('mod')==filter_mod[mod])
        .filter(
            pl.col('full_context').str.to_uppercase().str.contains(strMatcher[mod])
        ).with_columns([
            pl.lit(exp).alias('sample'),
            pl.col('mod').replace(modlookup).alias('mod'),
            pl.lit('r10.4.1').alias('flowcell'),
            pl.lit(f"dorado_{model}_{mod}").alias('tool'),
            pl.lit(f"{model}_{mod}").alias('model'),
            pl.lit(acc).alias('acc'),
            pl.lit(f"{sr}kHz").alias('sample_rate'),
            pl.lit(species).alias('species')
        ]).select(select)
    )

df = load_modkit(sys.argv[1])
df.write_csv(sys.argv[2], separator='\t', include_header=True)