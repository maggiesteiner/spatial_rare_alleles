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

elif snakemake.params.label=="2epsilon":
    filtered_df = df[(df.iloc[:,1]==1) & (df.iloc[:,-4]==True)] # filter to used in PCA YES and within_2epsilon TRUE
    sampled_df = filtered_df.sample(n=snakemake.params.numsamp)
    new_df = pd.DataFrame({"X1": sampled_df.iloc[:, 0], "X2": sampled_df.iloc[:, 0], "Clust": snakemake.params.label})    
    new_df.to_csv(snakemake.output[0], sep="\t", index=False)

elif snakemake.params.label=="3epsilon":
    filtered_df = df[(df.iloc[:,1]==1) & (df.iloc[:,-3]==True)] # filter to used in PCA YES and within_3epsilon TRUE
    sampled_df = filtered_df.sample(n=snakemake.params.numsamp)
    new_df = pd.DataFrame({"X1": sampled_df.iloc[:, 0], "X2": sampled_df.iloc[:, 0], "Clust": snakemake.params.label})    
    new_df.to_csv(snakemake.output[0], sep="\t", index=False)
    
elif snakemake.params.label=="4epsilon":
    filtered_df = df[(df.iloc[:,1]==1) & (df.iloc[:,-2]==True)] # filter to used in PCA YES and within_4epsilon TRUE
    sampled_df = filtered_df.sample(n=snakemake.params.numsamp)
    new_df = pd.DataFrame({"X1": sampled_df.iloc[:, 0], "X2": sampled_df.iloc[:, 0], "Clust": snakemake.params.label})    
    new_df.to_csv(snakemake.output[0], sep="\t", index=False)

elif snakemake.params.label=="10epsilon":
    filtered_df = df[(df.iloc[:,1]==1) & (df.iloc[:,-1]==True)] # filter to used in PCA YES and within_10epsilon TRUE
    sampled_df = filtered_df.sample(n=snakemake.params.numsamp)
    new_df = pd.DataFrame({"X1": sampled_df.iloc[:, 0], "X2": sampled_df.iloc[:, 0], "Clust": snakemake.params.label})    
    new_df.to_csv(snakemake.output[0], sep="\t", index=False)

