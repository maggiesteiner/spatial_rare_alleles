import pandas as pd
import glob
import numpy as np
import re

group = snakemake.params.group
n = snakemake.params.nsamp
type = snakemake.params.type

aclab = "AC_"+group
file_list = glob.glob("results/sfs_files/*merged*"+group+"*nsamp"+str(n)+'_*')#glob.glob("results/sfs_files/*merged*"+label+"*.gz")
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

    point_mutations = data[(data['REF'].apply(len) == 1) & (data['ALT'].apply(len) == 1)]
    
    iter_list.append(int(re.search(r'_iter(\d+)_', file).group(1)) if re.search(r'_iter(\d+)_', file) else [])
    
    if type == "all":
        segsites_list.append(point_mutations[(point_mutations[aclab] > 0) & (point_mutations[aclab] < n)].groupby(['CHROM', 'POS']).size().shape[0])
        mono_list.append(point_mutations[point_mutations[aclab] == 0].groupby(['CHROM', 'POS']).size().shape[0])
        ac_mean_list.append(np.mean(point_mutations[aclab]))
        ac_sd_list.append(np.std(point_mutations[aclab]))
        seg_ac_mean_list.append(np.mean(point_mutations[aclab][point_mutations[aclab] > 0]))
        seg_ac_sd_list.append(np.std(point_mutations[aclab][point_mutations[aclab] > 0]))
        singletons_list.append(point_mutations[point_mutations[aclab] == 1].groupby(['CHROM', 'POS']).size().shape[0])
    elif type == "lof":
        filtered_data = point_mutations[point_mutations['Annot'] == 'LoF']
        segsites_list.append(filtered_data[(filtered_data[aclab] > 0) & (filtered_data[aclab] < n)].groupby(['CHROM', 'POS']).size().shape[0])
        mono_list.append(filtered_data[filtered_data[aclab] == 0].groupby(['CHROM', 'POS']).size().shape[0])
        ac_mean_list.append(np.mean(filtered_data[aclab]))
        ac_sd_list.append(np.std(filtered_data[aclab]))
        seg_ac_mean_list.append(np.mean(filtered_data[aclab][filtered_data[aclab] > 0]))
        seg_ac_sd_list.append(np.std(filtered_data[aclab][filtered_data[aclab] > 0]))
        singletons_list.append(filtered_data[filtered_data[aclab] == 1].groupby(['CHROM', 'POS']).size().shape[0])


#    iter_list.append(int(re.search(r'_iter(\d+)_', file).group(1)) if re.search(r'_iter(\d+)_', file) else [])
#    if type == "all":
#    	segsites_list.append(data[(data[aclab]>0) & (data[aclab]<n)].groupby(['CHROM','POS']).size().shape[0])
#    	mono_list.append(data[data[aclab]==0].groupby(['CHROM','POS']).size().shape[0])
#    	ac_mean_list.append(np.mean(data[aclab]))
#    	ac_sd_list.append(np.std(data[aclab]))    
#    	seg_ac_mean_list.append(np.mean(data[aclab][data[aclab]>0]))
#    	seg_ac_sd_list.append(np.std(data[aclab][data[aclab]>0]))    
#    	singletons_list.append(data[data[aclab]==1].groupby(['CHROM','POS']).size().shape[0])
#    elif type == "lof":
#        filtered_data = data[data['Annot'] == 'LoF']
#        segsites_list.append(filtered_data[(filtered_data[aclab]>0) & (filtered_data[aclab]<n)].groupby(['CHROM','POS']).size().shape[0])
#        mono_list.append(filtered_data[filtered_data[aclab]==0].groupby(['CHROM','POS']).size().shape[0])
#        ac_mean_list.append(np.mean(filtered_data[aclab]))
#        ac_sd_list.append(np.std(filtered_data[aclab]))
#        seg_ac_mean_list.append(np.mean(filtered_data[aclab][filtered_data[aclab]>0]))
#        seg_ac_sd_list.append(np.std(filtered_data[aclab][filtered_data[aclab]>0]))
#        singletons_list.append(filtered_data[filtered_data[aclab]==1].groupby(['CHROM','POS']).size().shape[0])



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

