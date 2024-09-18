import pandas as pd
import glob
import numpy as np
import re

label = snakemake.params.label
n = snakemake.params.nsamp
type = snakemake.params.type

aclab = "AC_"+label
file_list = glob.glob("results/sfs_files/*merged*"+label+"*")

sfs_list = []
max_count = 1000
all_acs = np.arange(max_count+1)

for file in file_list:
    seed = int(re.search(r'_seed(\d+)_', file).group(1)) if re.search(r'_seed(\d+)_', file) else []
    data = pd.read_csv(file,sep='\t',compression='gzip')
    data[aclab] = pd.to_numeric(data[aclab], errors='coerce')
    if type == "all":
        ac, ac_counts = np.unique(data[aclab], return_counts=True)
    elif type == "LOF":
        ac, ac_counts = np.unique(data.loc[data['Annot']=="LoF",aclab], return_counts=True)
    counts_trunc = np.zeros(max_count+1,dtype=int)
    ac = ac.astype(int)
    
    valid_ac_indices = (ac >= 0) & (ac <= max_count)
    valid_ac = ac[valid_ac_indices]
    valid_ac_counts = ac_counts[valid_ac_indices]
    
    counts_trunc[valid_ac] = valid_ac_counts

    #counts_trunc[ac] = ac_counts[:max_count+1]
    row_temp = np.insert(counts_trunc,0,seed)
    sfs_list.append(row_temp)

column_names=["seed"]+[str(x) for x in all_acs]
df = pd.DataFrame(sfs_list,columns=column_names)
df.to_csv(snakemake.output[0],sep='\t',index=False)
