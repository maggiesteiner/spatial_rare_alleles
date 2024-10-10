import pandas as pd
import numpy as np
from sklearn.neighbors import KernelDensity

df = pd.read_csv(snakemake.input[0])

filtered_df_pca = df[(df['used_in_pca'] == 1) & (df['PC1'].notna()) & (df['PC2'].notna())]
coords_pca = filtered_df_pca.iloc[:,7:9].values

kde_pca = KernelDensity()
kde_pca.fit(coords_pca)
log_dens_pca = kde_pca.score_samples(coords_pca)
filtered_df_pca['pc_density'] = np.exp(log_dens_pca)

filtered_df_pca.to_csv(snakemake.output[0],index=False)

print("KDE for PCA fitted")

filtered_df_geo = df[(df['used_in_pca'] == 1) & (df['birth_east_coord'].notna()) & (df['birth_north_coord'].notna()) & (df['birth_UKelsewhere']!='Elsewhere')]
coords_geo = filtered_df_geo.iloc[:,3:5].values

kde_geo = KernelDensity()
kde_geo.fit(coords_geo)
log_dens_geo = kde_geo.score_samples(coords_geo)
filtered_df_geo['geo_density'] = np.exp(log_dens_geo)

print("KDE for Geo fitted")

filtered_df_geo.to_csv(snakemake.output[1],index=False)
