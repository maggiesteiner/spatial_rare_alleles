import pandas as pd
import glob
import numpy as np
import re

label = snakemake.params.label
n = snakemake.params.nsamp

aclab = "AC_"+label
file_list = glob.glob("results/sfs_files/*merged*"+label+"*")

sfs_list = []

for file in file_list:
    iter = int(re.search(r'_iter(\d+)_', file).group(1)) if re.search(r'_iter(\d+)_', file) else []
    data = pd.read_csv(file,sep='\t',compression='gzip')
    data[aclab] = pd.to_numeric(data[aclab], errors='coerce')
    ac, ac_counts = np.unique(data[aclab], return_counts=True)
    counts_truc = ac_counts[:1001]
    row_temp = np.insert(counts_truc,0,iter)
    sfs_list.append(row_temp)

column_names=["iter"]+[str(x) for x in ac[:1000]]#1]]
df = pd.DataFrame(sfs_list,columns=column_names)
df.to_csv(snakemake.output[0],sep='\t',index=False)

