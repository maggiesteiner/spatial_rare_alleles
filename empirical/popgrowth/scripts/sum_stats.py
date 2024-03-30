import pandas as pd
import glob
import numpy as np
import re

label = snakemake.params.label
n = snakemake.params.nsamp
aclab = "AC_"+label
file_list = glob.glob("results/sfs_files/*merged*"+label+"*")
segsites_list = []
mono_list = []
ac_mean_list = []
ac_sd_list = []
seg_ac_mean_list = []
seg_ac_sd_list = []
singletons_list = []
iter_list = []
for file in file_list: 
    data = pd.read_csv(file,sep='\t',compression='gzip')
    data[aclab] = pd.to_numeric(data[aclab], errors='coerce')
    iter_list.append(int(re.search(r'_iter(\d+)_', file).group(1)) if re.search(r'_iter(\d+)_', file) else [])
    segsites_list.append(data[(data[aclab]>0) & (data[aclab]<n)].groupby(['CHROM','POS']).size().shape[0])
    mono_list.append(data[data[aclab]==0].groupby(['CHROM','POS']).size().shape[0])
    ac_mean_list.append(np.mean(data[aclab]))
    ac_sd_list.append(np.std(data[aclab]))    
    seg_ac_mean_list.append(np.mean(data[aclab][data[aclab]>0]))
    seg_ac_sd_list.append(np.std(data[aclab][data[aclab]>0]))    
    singletons_list.append(data[data[aclab]==1].groupby(['CHROM','POS']).size().shape[0])

df = pd.DataFrame({
    'iter':iter_list,
    'segsites':segsites_list,
    'monosites':mono_list,
    'ac_mean':ac_mean_list,
    'ac_sd':ac_sd_list,
    'seg_ac_mean':seg_ac_mean_list,
    'seg_ac_sd':seg_ac_sd_list,
    'singletons':singletons_list})

df.to_csv(snakemake.output[0],sep='\t',index=False)

