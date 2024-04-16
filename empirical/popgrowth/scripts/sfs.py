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

def process_ac(data, aclab, type):
    if type == "all":
        ac, ac_counts = np.unique(data[aclab], return_counts=True)
    elif type == "lof":
        ac, ac_counts = np.unique(data.loc[data['Annot'] == "LoF", aclab], return_counts=True)
    return ac, ac_counts

for file in file_list:
    iter = int(re.search(r'_numiter(\d+).sfs', file).group(1)) if re.search(r'_numiter(\d+).', file) else []
    # print(file)
    data = pd.read_csv(file, sep='\t', compression='gzip')
    data[aclab] = pd.to_numeric(data[aclab], errors='coerce')
    point_mutations = data[(data['REF'].apply(len) == 1) & (data['ALT'].apply(len) == 1)]
    
    ac, ac_counts = process_ac(point_mutations, aclab, type)
    
    row_temp = np.insert(ac_counts, 0, iter)
    sfs_list.append(row_temp)

max_ac = max(len(row) - 1 for row in sfs_list)
column_names = ["iter"] + [str(x) for x in range(max_ac)]
df = pd.DataFrame(sfs_list, columns=column_names)
df.to_csv(snakemake.output[0], sep='\t', index=False)

#for file in file_list:
#    iter = int(re.search(r'_numiter(\d+).sfs', file).group(1)) if re.search(r'_numiter(\d+).', file) else []
#    #print(file)
#    data = pd.read_csv(file,sep='\t',compression='gzip')
#    data[aclab] = pd.to_numeric(data[aclab], errors='coerce')
#    point_mutations = data[(data['REF'].apply(len) == 1) & (data['ALT'].apply(len) == 1)]
#    if type == "all":
#        ac, ac_counts = np.unique(point_mutations[aclab], return_counts=True)
#    elif type == "lof":
#        ac, ac_counts = np.unique(point_mutations.loc[point_mutations['Annot']=="LoF",aclab], return_counts=True)
    #ac, ac_counts = np.unique(data[aclab], return_counts=True)#
#    counts_truc = ac_counts[:1001]
#    row_temp = np.insert(counts_truc,0,iter)
#    sfs_list.append(row_temp)

#column_names=["iter"]+[str(x) for x in ac[:1001]]
#column_names = ["iter"]+[str(x) for x in ac[:min(1001, len(ac))]]
#df = pd.DataFrame(sfs_list,columns=column_names)
#df.to_csv(snakemake.output[0],sep='\t',index=False)

