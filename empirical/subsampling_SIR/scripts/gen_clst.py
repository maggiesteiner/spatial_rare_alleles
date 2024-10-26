import pandas as pd

fam=pd.read_csv(snakemake.input[0],sep=' ',header=None)

fam.columns = ['col1']
fam['col2'] = fam['col1']
fam['col3'] = 0
fam.to_csv(snakemake.output[0], sep=" ", index=False, header=False)