import pandas as pd
import numpy as np

meta = pd.read_csv(snakemake.input[0])
wes = pd.read_csv(snakemake.input[1],header=None)
evecs = pd.read_csv(snakemake.input[2],delim_whitespace=True)

# format metadata
meta = meta.rename(columns={'participant.eid':'id'})

# undo codings
meta['participant.p54_i0'] = meta['participant.p54_i0'].astype(str)
center_coding = pd.read_csv(snakemake.input[3],delimiter="\t",names=["participant.p54_i0","assessment_center"],header=0)
center_dict = dict(zip(center_coding['participant.p54_i0'].astype(str),center_coding['assessment_center']))
meta['assessment_center'] = meta['participant.p54_i0'].map(center_dict)
meta = meta.drop('participant.p54_i0',axis=1)

meta['participant.p1647_i0'] = meta['participant.p1647_i0'].astype(str)
meta['participant.p1647_i0'] = meta['participant.p1647_i0'].str.replace(r'\.0$', '',regex=True)
birth_UKelsewhere_coding = pd.read_csv(snakemake.input[4],delimiter="\t",usecols=[0,1],names=["participant.p1647_i0","birth_UKelsewhere"],header=0)
birth_UKelsewhere_dict = dict(zip(birth_UKelsewhere_coding['participant.p1647_i0'].astype(str),birth_UKelsewhere_coding['birth_UKelsewhere']))
meta['birth_UKelsewhere'] = meta['participant.p1647_i0'].map(birth_UKelsewhere_dict)
meta = meta.drop('participant.p1647_i0',axis=1)

# in coordinates, replace -1 with NA
meta = meta.rename(columns={'participant.p130_i0':'birth_east_coord','participant.p129_i0':'birth_north_coord','participant.p22006':'whitebritish','participant.p22020':'used_in_pca'})
meta.loc[meta['birth_east_coord'] == -1, 'birth_east_coord'] = pd.NA
meta.loc[meta['birth_north_coord'] == -1, 'birth_north_coord'] = pd.NA

# fill NAs
meta['whitebritish'] = meta['whitebritish'].fillna(0)
meta['used_in_pca'] = meta['used_in_pca'].fillna(0)

# merge data
wes = wes.rename(columns={wes.columns[0]:'id'})
filtered_meta = pd.merge(meta,wes,on='id')

evecs = evecs.rename(columns={'IID':'id'})
meta_pcs = pd.merge(filtered_meta,evecs.iloc[:,1:],on='id',how='left')

# calculate median centroid
meta_pcs = meta_pcs[(meta_pcs['used_in_pca']==True)]
centroid_pca = meta_pcs.iloc[:,7:9].median()

# calculate distances from centroid
distances = np.sqrt(((meta_pcs.iloc[:,7:9] - centroid_pca) ** 2).sum(axis=1))
meta_pcs['distance_to_centroid_pca'] = distances

epsilon_pca = snakemake.params.eps_pca

meta_pcs['within_1epsilon_pca'] = meta_pcs['distance_to_centroid_pca'] < 1 * epsilon_pca

# save file
meta_pcs.to_csv("metadata/metadata_cleaned_wes.csv",index=False)
