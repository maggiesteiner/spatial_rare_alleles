import pandas as pd
import glob
import numpy as np
import re

label = snakemake.params.label
n = snakemake.params.nsamp
type = snakemake.params.type
aclab = "AC_"+label
anlab = "AN_"+label
file_list = glob.glob("results/sfs_files/*merged*"+label+"_*")
segsites_list = []
mono_list = []
#ac_mean_list = []
#ac_sd_list = []
#seg_ac_mean_list = []
#seg_ac_sd_list = []
singletons_list = []
freq_mean_list = []
freq_sd_list = []
seg_freq_mean_list = []
seg_freq_sd_list = []
het_list = []
seed_list = []
for file in file_list: 
    data = pd.read_csv(file,sep='\t',compression='gzip')
    if type == 'LOF':
        data = data.loc[data['Annot']=="LoF",:]
    data[aclab] = pd.to_numeric(data[aclab], errors='coerce')
    data[anlab] = pd.to_numeric(data[anlab], errors='coerce')
    data = data[data[anlab] > 0]
    seed_list.append(int(re.search(r'_seed(\d+)_', file).group(1)) if re.search(r'_seed(\d+)_', file) else [])
    segsites_list.append(data[(data[aclab]>0) & (data[aclab]<n)].groupby(['CHROM','POS']).size().shape[0])
    mono_list.append(data[data[aclab]==0].groupby(['CHROM','POS']).size().shape[0])
    #ac_mean_list.append(np.mean(data[aclab]))
    #ac_sd_list.append(np.std(data[aclab]))    
    #seg_ac_mean_list.append(np.mean(data[aclab][data[aclab]>0]))
    #seg_ac_sd_list.append(np.std(data[aclab][data[aclab]>0]))    
    singletons_list.append(data[data[aclab]==1].groupby(['CHROM','POS']).size().shape[0])
    data['het'] = 2 * (data[aclab] / data[anlab]) * (1 - (data[aclab] / data[anlab]))
    het_list.append(data['het'].mean())
    data['allele_freq'] = data[aclab] / data[anlab]
    freq_mean_list.append(np.mean(data['allele_freq']))
    freq_sd_list.append(np.std(data['allele_freq']))
    seg_freq_mean_list.append(np.mean(data[data[aclab] > 0]['allele_freq']))
    seg_freq_sd_list.append(np.std(data[data[aclab] > 0]['allele_freq']))


df = pd.DataFrame({
    'seed':seed_list,
    'segsites':segsites_list,
    'monosites':mono_list,
    #'ac_mean':ac_mean_list,
    #'ac_sd':ac_sd_list,
    #'seg_ac_mean':seg_ac_mean_list,
    #'seg_ac_sd':seg_ac_sd_list,
    'singletons':singletons_list,
    'heterozygosity': het_list,
    'af_mean':freq_mean_list,
    'af_sd':freq_sd_list,
    'seg_af_mean':seg_freq_mean_list,
    'seg_af_sd':seg_freq_sd_list})

df.to_csv(snakemake.output[0],sep='\t',index=False)
