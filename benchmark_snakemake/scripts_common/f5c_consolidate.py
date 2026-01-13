import polars as pl
import sys, re
from pathlib import Path

input_file  = sys.argv[1]
output_file = sys.argv[2]

t = pl.read_csv(input_file, separator='\t', has_header=False, new_columns=['chrom', 'p1', 'p2', 'coverage', 's', 'strand', 'M', 'UM', 'per', 'full_context'])
select  = ['chrom', 'p1', 'p2', 'mod', 'coverage', 'strand', 'M', 'UM', 'per', 'flowcell', 'tool', 'model', 'sample', 'acc', 'sample_rate', 'species', 'full_context']

experiment, sr, acc, ver = re.match(r'(.*)_([45]kHz)_(hac|sup)_(v[45.2]+r[123])', Path(input_file).stem).group(1,2,3,4)
species = experiment.split('_')[0]
strandstate = '_stranded' if Path(input_file).stem.count('.stranded')==1 else ''

t = (
    t.with_columns([
        pl.lit('5mC').alias('mod'),
        pl.lit('r10.4.1').alias('flowcell'),
        (pl.col('M')+pl.col('UM')).alias('coverage'),
        pl.lit(experiment).alias('sample'),
        pl.lit(f'f5c{strandstate}').alias('tool'),
        pl.lit(acc).alias('acc'),
        pl.lit('.').alias('model'),
        pl.lit(f"{sr}").alias('sample_rate'),
        pl.lit(species).alias('species'),
        (pl.col('per')*100).alias('per')
    ])
)
t.select(select).write_csv(output_file, separator='\t', include_header=True)