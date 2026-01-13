
import polars as pl
import sys

df = (
    pl.read_csv(input[0], separator='\t', has_header=True)
    .with_columns(
        pl.when(pl.col('full_context').str.to_uppercase().str.contains(r'\w{5}CG\w{10}')).then(1)
        .when(pl.col('full_context').str.to_uppercase().str.contains(r'\w{5}C[ATC]G\w{9}')).then(2)
        .when(pl.col('full_context').str.to_uppercase().str.contains(r'\w{5}C[ATC][ATC]\w{9}')).then(3)
        .otherwise(0)
        .cast(pl.String)
        .replace({ '0': 'other', '1': 'CpG', '2': 'CHG', '3': 'CHH' })
        .alias('context')
    )
)
df.write_csv(output[0], separator='\t', include_header=True)