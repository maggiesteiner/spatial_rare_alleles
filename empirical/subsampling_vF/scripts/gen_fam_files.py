import pandas as pd
import numpy as np

df = pd.read_csv(snakemake.input[0])
seed = int(snakemake.params.seed)
np.random.seed(seed)

pca_column = "used_in_pca"
birth_column = "birth_UKelsewhere"
east_coord_column = "birth_east_coord"

pca_columns = {
    "closest10kpca": "closest_10k_pca",
    "1epsilonpca": "within_1epsilon_pca",
    "5epsilonpca": "within_5epsilon_pca",
    "10epsilonpca": "within_10epsilon_pca",
    "50epsilonpca": "within_50epsilon_pca",
    "100epsilonpca": "within_100epsilon_pca"
}
geo_columns = {
    "closest10kgeo": "closest_10k_geo",
    "5epsilongeo": "within_5epsilon_geo",
    "10epsilongeo": "within_10epsilon_geo",
    "50epsilongeo": "within_50epsilon_geo",   
}

# Filtering and sampling based on the label
label = snakemake.params.label
filtered_df = pd.DataFrame()  

if label == "uniform":
    filtered_df = df[df[pca_column] == 1]  # Filter to 'used_in_pca == 1'
elif label in pca_columns:
    filtered_df = df[(df[pca_column] == 1) & (df[pca_columns[label]] == True)] # Filter to `used_in_pca ==1` and label
elif label in geo_columns:
    filtered_df = df[(df[pca_column] == 1) & (df[geo_columns[label]] == True) & (df['within_1epsilon_pca']==True) & (df[east_coord_column].notna()) & (df[birth_column]!='Elsewhere')] 
elif label == "uniformgeo":
    filtered_df = df[(df[pca_column] == 1) & (df['within_1epsilon_pca']==True) & (df[east_coord_column].notna()) & (df[birth_column]!='Elsewhere')]
else:
    raise ValueError(f"Unexpected label: {label}")

# Check if filtering worked
if label in pca_columns and not filtered_df[pca_columns[label]].all():
    raise ValueError(f"Filtering error: Not all entries in {pca_columns[label]} are True.")

if label in geo_columns and not filtered_df[geo_columns[label]].all():
    raise ValueError(f"Filtering error: Not all entries in {geo_columns[label]} are True.")


sampled_df = filtered_df.sample(n=snakemake.params.numsamp)
new_df = pd.DataFrame({"X1": sampled_df["id"], "X2": sampled_df["id"], "Clust": label})
new_df.to_csv(snakemake.output[0], sep="\t", index=False)

