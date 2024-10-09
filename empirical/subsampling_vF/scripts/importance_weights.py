import pandas as pd
from sklearn.neighbors import KernelDensity

df = pd.read_csv(snakemake.input[0])
coords_pca = df.iloc[:,7:9].values
coords_geo = df.iloc[:,3:5].values

kde_pca = KernelDensity()
kde_pca.fit(coords_pca)
log_dens_pca = kde_pca.score_samples(coords_pca)
df['pc_density'] = np.exp(log_dens_pca)

print("KDE for PCA fitted")

kde_geo = KernelDensity()
kde_geo.fit(coords_geo)
log_dens_geo = kde_geo.score_samples(coords_geo)
data['geo_density'] = np.exp(log_dens_geo)

print("KDE for Geo fitted")

df.to_csv(snakemake.output[0],index=False)