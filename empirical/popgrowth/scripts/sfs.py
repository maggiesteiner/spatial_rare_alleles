import pandas as pd
import glob
import numpy as np
import re

group = snakemake.params.group
n = snakemake.params.nsamp
type = snakemake.params.type

aclab = "AC_"+group
file_list = glob.glob("results/sfs_files/*merged*"+group+"*nsamp"+str(n)+'_*.gz')

sfs_list = []

for file in file_list:
    iter = int(re.search(r'_numiter(\d+).sfs', file).group(1)) if re.search(r'_numiter(\d+).', file) else []
    #print(file)
    data = pd.read_csv(file,sep='\t',compression='gzip')
    data[aclab] = pd.to_numeric(data[aclab], errors='coerce')
    point_mutations = data[(data['REF'].apply(len) == 1) & (data['ALT'].apply(len) == 1)]
    if type == "all":
        ac, ac_counts = np.unique(point_mutations[aclab], return_counts=True)
    elif type == "lof":
        ac, ac_counts = np.unique(point_mutations.loc[point_mutations['Annot']=="LoF",aclab], return_counts=True)
    #ac, ac_counts = np.unique(data[aclab], return_counts=True)
    counts_truc = ac_counts[:1001]
    row_temp = np.insert(counts_truc,0,iter)
    sfs_list.append(row_temp)

column_names=["iter"]+[str(x) for x in ac[:1001]]
#column_names = ["iter"]+[str(x) for x in ac[:min(1001, len(ac))]]
df = pd.DataFrame(sfs_list,columns=column_names)
df.to_csv(snakemake.output[0],sep='\t',index=False)

