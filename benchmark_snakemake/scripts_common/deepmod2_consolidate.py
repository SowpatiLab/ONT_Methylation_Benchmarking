import polars as pl
from pathlib import Path
import re, sys

input_file  = sys.argv[1]
output_file = sys.argv[2]
deepmod2_model = sys.argv[3]

experiment, sr, acc, ver = re.match(r'(.*)_([45]kHz)_(hac|sup)_(v[45.2]+r[123])', Path(input_file).stem).group(1,2,3,4)
species = experiment.split('_')[0]

t = pl.read_csv(input_file, separator='\t', has_header=True, new_columns=['chrom', 'p1', 'p2', 'coverage', 'is_cpg', 'strand', 'M', 'UM', 'per', 'full_context'])
select  = ['chrom', 'p1', 'p2', 'mod', 'coverage', 'strand', 'M', 'UM', 'per', 'flowcell', 'tool', 'model', 'sample', 'acc', 'sample_rate', 'species', 'full_context']


(
    t.with_columns([
        pl.lit('5mC').alias('mod'),
        pl.lit('r10.4.1').alias('flowcell'),
        pl.lit(experiment).alias('sample'),
        pl.lit(f'deepmod2_{deepmod2_model}').alias('tool'),
        pl.lit(acc).alias('acc'),
        pl.lit(deepmod2_model).alias('model'),
        pl.lit(f"{sr}").alias('sample_rate'),
        pl.lit(species).alias('species'),
        (pl.col('per')*100).alias('per')
    ])
    .select(select)
).write_csv(output_file, separator='\t', include_header=True)
