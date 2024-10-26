import pandas as pd
import numpy as np

weights=pd.read_csv(snakemake.input[0],sep=' ',header=None)
fam=pd.read_csv(snakemake.input[1],sep=' ',header=None)

fam.columns = ['col1']
fam['col2'] = fam['col1']
weights.columns = ['col1', 'col2', 'weight_value']
merged_df = pd.merge(fam, weights, on=['col1', 'col2'], how='left')
merged_df['weight_value'] = merged_df['weight_value'].fillna(0)


sumval = merged_df['weight_value'].sum()
print(f'Sum is: {sumval}')



merged_df.to_csv(snakemake.output[0], sep=" ", index=False, header=False)