import pandas as pd
import numpy as np

# function for Gaussian densities
def gaussian_2d(r, sigma):
    return (1 / (2 * np.pi * sigma ** 2)) * np.exp(-r ** 2 / (2 * sigma ** 2))

# function to map to bin
def map_to_bin(x,y,xedges,yedges):
    return [np.digitize(x,xedges)-1,np.digitize(y,yedges)-1]

# function to get bin center
def get_bin_center(bin_idx,xedges,yedges):
    x_idx, y_idx = bin_idx
    center_x = (xedges[x_idx] + xedges[x_idx + 1]) / 2
    center_y = (yedges[y_idx] + yedges[y_idx + 1]) / 2
    return [center_x, center_y]

# function for counts given coordinate
def get_count(x, y, counts, xedges, yedges):
    x_bin_index = np.digitize(x, xedges) - 1
    y_bin_index = np.digitize(y, yedges) - 1

    x_bin_index = np.clip(x_bin_index, 0, counts.shape[0] - 1)
    y_bin_index = np.clip(y_bin_index, 0, counts.shape[1] - 1)

    return counts[x_bin_index, y_bin_index]

def parse_center_PCA(s):
    x_index = s.index('X')
    y_index = s.index('Y')
    locX = int(s[x_index + 1:y_index])
    locY = int(s[y_index + 1:])
    return [locX,locY]

df = pd.read_csv(snakemake.input[0])

df_pca = df[(df['used_in_pca']==1)]
# histogram
x = df_pca['PC1']
y = df_pca['PC2']
counts, xedges, yedges = np.histogram2d(x,y,bins=20)

# read in snakemake params
centers = snakemake.params.center_list#['center1','center2','center3']
w_list = snakemake.params.w_list#[25000,50000,100000]

# distances
center_bin_1 = parse_center_PCA(centers[0])
# center_bin_2 = parse_center_PCA(centers[1])
# center_bin_3 = parse_center_PCA(centers[2])
center_coord_1 = get_bin_center(center_bin_1,xedges,yedges)
# center_coord_2 = get_bin_center(center_bin_2,xedges,yedges)
# center_coord_3 = get_bin_center(center_bin_3,xedges,yedges)
df_pca[centers[0]]= np.sqrt(((df_pca.iloc[:,7:9] - center_coord_1) ** 2).sum(axis=1))
# df_pca[centers[1]]= np.sqrt(((df_pca.iloc[:,7:9] - center_coord_2) ** 2).sum(axis=1))
# df_pca[centers[2]]= np.sqrt(((df_pca.iloc[:,7:9] - center_coord_3) ** 2).sum(axis=1))

# binned counts
df_pca['binned_counts'] = df_pca.apply(lambda row: get_count(row['PC1'],row['PC2'],counts,xedges,yedges),axis=1)

# filter
thresh=15
df_pca_filt = df_pca[df_pca['binned_counts']>thresh]
N=len(df_pca_filt)

# weights
df_pca_filt['freq_binned']=df_pca_filt['binned_counts']/N

# uniform weights
df_pca_filt['uniformpca'] = 1/df_pca_filt['freq_binned']
df_pca_filt['uniformpca'] = df_pca_filt['uniformpca'].apply(lambda x: 0 if abs(x) < 1e-100 else x)


unif_weights = pd.DataFrame({
            'id1': df_pca_filt['id'].astype(int),
            'id2': df_pca_filt['id'].astype(int),
            'weights': df_pca_filt['uniformpca']
        })
unif_weights.to_csv('data/weights/uniformpca.weights', sep=" ", index=False, header=False)

# centers = ['center1','center2','center3']
# w_list = [25000,50000,100000]

for c in centers:
    for w in w_list:
        colname = f'{c}pca{w}'
        df_pca_filt[colname] = df_pca_filt.apply(lambda row: gaussian_2d(row[c],w),axis=1)/df_pca_filt['freq_binned']
        df_pca_filt[colname] = df_pca_filt[colname].apply(lambda x: 0 if abs(x) < 1e-100 else x)
        temp = pd.DataFrame({
            'id1': df_pca_filt['id'].astype(int),
            'id2': df_pca_filt['id'].astype(int),
            'weights': df_pca_filt[colname]
        })
        temp.to_csv(f'data/weights/{colname}.weights', sep=" ", index=False, header=False)


