import pandas as pd
import numpy as np

df = pd.read_csv('metadata/metadata_cleaned_wes.csv')

pca_column = "used_in_pca"


filtered_df = pd.DataFrame()  
filtered_df = df[(df[pca_column] == 1) & (df['within_1epsilon_pca'] == True)]  # Filter to 'used_in_pca == 1' and within 1eps PCA
new_df = pd.DataFrame({"X1": filtered_df["id"], "X2": filtered_df["id"], "Clust": "within_1epsilon_pca_all"})
new_df.to_csv('data/id_lists/within_1epsilon_pca_all.tsv', sep="\t", index=False)

filtered_df = pd.DataFrame()
filtered_df = df[(df[pca_column] == 1) & (df['within_5epsilon_geo'] == True)]  # Filter to 'used_in_pca == 1' and within 5eps Geo
new_df = pd.DataFrame({"X1": filtered_df["id"], "X2": filtered_df["id"], "Clust": "within_5epsilon_geo_all"})
new_df.to_csv('data/id_lists/within_5epsilon_geo_all.tsv', sep="\t", index=False)

filtered_df = pd.DataFrame()
filtered_df = df[(df[pca_column] == 1) & (df['within_1epsilon_pca'] == False)]  # Filter to 'used_in_pca == 1' and !within 1eps PCA
new_df = pd.DataFrame({"X1": filtered_df["id"], "X2": filtered_df["id"], "Clust": "notwithin_1epsilon_pca_all"})
new_df.to_csv('data/id_lists/notwithin_1epsilon_pca_all.tsv', sep="\t", index=False)

filtered_df = pd.DataFrame()
filtered_df = df[(df[pca_column] == 1) & (df['within_5epsilon_geo'] == False)]  # Filter to 'used_in_pca == 1' and !within 5eps Geo
new_df = pd.DataFrame({"X1": filtered_df["id"], "X2": filtered_df["id"], "Clust": "notwithin_5epsilon_geo_all"})
new_df.to_csv('data/id_lists/notwithin_5epsilon_geo_all.tsv', sep="\t", index=False)
