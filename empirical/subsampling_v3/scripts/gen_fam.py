import pandas as pd
import numpy as np


df = pd.read_csv(snakemake.input[0],skiprows=1)
seed = int(snakemake.params.seed)
np.random.seed(seed)

if snakemake.params.label=="uniform":
    filtered_df = df[df.iloc[:,2]==1] # filter to used in PCA YES
    sampled_df = filtered_df.sample(n=snakemake.params.numsamp)
    new_df = pd.DataFrame({"X1": sampled_df.iloc[:, 0], "X2": sampled_df.iloc[:, 0], "Clust": snakemake.params.label})
    new_df.to_csv(snakemake.output[0], sep="\t", index=False)

elif snakemake.params.label=="whitebritish":
    filtered_df = df[(df.iloc[:,1]==1) & (df.iloc[:,2]==1)] # filter to used in PCA YES and white british YES
    sampled_df = filtered_df.sample(n=snakemake.params.numsamp)
    new_df = pd.DataFrame({"X1": sampled_df.iloc[:, 0], "X2": sampled_df.iloc[:, 0], "Clust": snakemake.params.label})    
    new_df.to_csv(snakemake.output[0], sep="\t", index=False)

elif snakemake.params.label=="newcastle":
    filtered_df = df[(df.iloc[:,1]==1) & (df.iloc[:,2]==1) & (df.iloc[:,3]==11009)] # filter to used in PCA YES and white british YES and NEWCASTLE
    sampled_df = filtered_df.sample(n=snakemake.params.numsamp)
    new_df = pd.DataFrame({"X1": sampled_df.iloc[:, 0], "X2": sampled_df.iloc[:, 0], "Clust": snakemake.params.label})
    new_df.to_csv(snakemake.output[0], sep="\t", index=False)    

elif snakemake.params.label=="wales":
    filtered_df = df[(df.iloc[:,1]==1) & (df.iloc[:,2]==1) & (df.iloc[:,6]==3)] # filter to used in PCA YES and white british YES and WALES
    sampled_df = filtered_df.sample(n=snakemake.params.numsamp)
    new_df = pd.DataFrame({"X1": sampled_df.iloc[:, 0], "X2": sampled_df.iloc[:, 0], "Clust": snakemake.params.label})
    new_df.to_csv(snakemake.output[0], sep="\t", index=False)

elif snakemake.params.label=="grid":
    filtered_df = df[(df.iloc[:,1]==1) & (df.iloc[:,2]==1) & (df.iloc[:,4]>=420000) & (df.iloc[:,4]<=450000) & (df.iloc[:,5]>=420000) & (df.iloc[:,5]<=450000)] # filter to used in PCA YES and white british YES and GRID
    sampled_df = filtered_df.sample(n=snakemake.params.numsamp)
    new_df = pd.DataFrame({"X1": sampled_df.iloc[:, 0], "X2": sampled_df.iloc[:, 0], "Clust": snakemake.params.label})
    new_df.to_csv(snakemake.output[0], sep="\t", index=False)
