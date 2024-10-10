import pandas as pd
import numpy as np

def gaussian_2d(r, sigma):
    return (1 / (2 * np.pi * sigma ** 2)) * np.exp(-r ** 2 / (2 * sigma ** 2))

df = pd.read_csv(snakemake.input[0])
df_IS_geo = pd.read_csv(snakemake.input[1])
df_IS_pca = pd.read_csv(snakemake.input[2])
seed = int(snakemake.params.seed)
np.random.seed(seed)

pca_column = "used_in_pca"
birth_column = "birth_UKelsewhere"
east_coord_column = "birth_east_coord"

IS_labels = [
    '5epsilonpcaIS',
    '10epsilonpcaIS',
    '50epsilonpcaIS',
    '100epsilonpcaIS',
    'uniformpcaIS',
    '5epsilongeoIS',
    '10epsilongeoIS',
    '50epsilongeoIS',
    'uniformgeoIS'
]

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

epsilon_geo = snakemake.params.eps_geo
epsilon_pca = snakemake.params.eps_pca

if label == "uniform":
    filtered_df = df[df[pca_column] == 1]  # Filter to 'used_in_pca == 1'
elif label in pca_columns:
    filtered_df = df[(df[pca_column] == 1) & (df[pca_columns[label]] == True)] # Filter to `used_in_pca ==1` and label
elif label in geo_columns:
    filtered_df = df[(df[pca_column] == 1) & (df[geo_columns[label]] == True) & (df['within_1epsilon_pca']==True) & (df[east_coord_column].notna()) & (df[birth_column]!='Elsewhere')] 
elif label == "uniformgeo":
    filtered_df = df[(df[pca_column] == 1) & (df['within_1epsilon_pca']==True) & (df[east_coord_column].notna()) & (df[birth_column]!='Elsewhere')]
elif label in IS_labels:
    if label == "uniformpcaIS":
        df_IS_pca['weights'] = 1/df_IS_pca['pc_density']
    if label == "5epsilonpcaIS":
        df_IS_pca['weights'] = df_IS_pca.apply(lambda row: gaussian_2d(row['distance_to_centroid_pca'],5*epsilon_pca), axis=1)/df_IS_pca['pc_density']
    if label == "10epsilonpcaIS":
        df_IS_pca['weights'] = df_IS_pca.apply(lambda row: gaussian_2d(row['distance_to_centroid_pca'],10*epsilon_pca), axis=1)/df_IS_pca['pc_density']
    if label == "50epsilonpcaIS":
        df_IS_pca['weights'] = df_IS_pca.apply(lambda row: gaussian_2d(row['distance_to_centroid_pca'],50*epsilon_pca), axis=1)/df_IS_pca['pc_density']
    if label == "100epsilonpcaIS":
        df_IS_pca['weights'] = df_IS_pca.apply(lambda row: gaussian_2d(row['distance_to_centroid_pca'],100*epsilon_pca), axis=1)/df_IS_pca['pc_density']
    if label == "uniformgeoIS":
        df_IS_geo = df_IS_geo[df_IS_geo['within_1epsilon_pca']==True]
        df_IS_geo['weights'] = 1/df_IS_geo['geo_density']
    if label == "5epsilongeoIS":
        df_IS_geo = df_IS_geo[df_IS_geo['within_1epsilon_pca'] == True]
        df_IS_geo['weights'] = df_IS_geo.apply(
            lambda row: gaussian_2d(row['distance_to_centroid_geo'], 5 * epsilon_geo), axis=1) / df_IS_geo['geo_density']
    if label == "10epsilongeoIS":
        df_IS_geo = df_IS_geo[df_IS_geo['within_1epsilon_pca'] == True]
        df_IS_geo['weights'] = df_IS_geo.apply(
            lambda row: gaussian_2d(row['distance_to_centroid_geo'], 10 * epsilon_geo), axis=1) / df_IS_geo['geo_density']
    if label == "50epsilongeoIS":
        df_IS_geo = df_IS_geo[df_IS_geo['within_1epsilon_pca'] == True]
        df_IS_geo['weights'] = df_IS_geo.apply(
            lambda row: gaussian_2d(row['distance_to_centroid_geo'], 50 * epsilon_geo), axis=1) / df_IS_geo['geo_density']
else:
    raise ValueError(f"Unexpected label: {label}")

# # Check if filtering worked
# if label in pca_columns and not filtered_df[pca_columns[label]].all():
#     raise ValueError(f"Filtering error: Not all entries in {pca_columns[label]} are True.")
#
# if label in geo_columns and not filtered_df[geo_columns[label]].all():
#     raise ValueError(f"Filtering error: Not all entries in {geo_columns[label]} are True.")

if label not in IS_labels:
    sampled_df = filtered_df.sample(n=snakemake.params.numsamp)

elif label in IS_labels:
    if 'geo' in label:
        sampled_df = df_IS_geo.sample(n=snakemake.params.numsamp, weights='weights')
    elif 'pca' in label:
        sampled_df = df_IS_pca.sample(n=snakemake.params.numsamp, weights='weights')

new_df = pd.DataFrame({"X1": sampled_df["id"], "X2": sampled_df["id"], "Clust": label})
new_df.to_csv(snakemake.output[0], sep="\t", index=False)

