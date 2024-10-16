import pandas as pd
import numpy as np

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

def get_count(x, y, counts, xedges, yedges):
    x_bin_index = np.digitize(x, xedges) - 1
    y_bin_index = np.digitize(y, yedges) - 1

    x_bin_index = np.clip(x_bin_index, 0, counts.shape[0] - 1)
    y_bin_index = np.clip(y_bin_index, 0, counts.shape[1] - 1)

    return counts[x_bin_index, y_bin_index]

df = pd.read_csv(snakemake.input[0])
seed = int(snakemake.params.seed)

# Sampling based on the label
label = snakemake.params.label

if 'geo' in label:
    ### geo sampling
    # subsample to individuals (1) with birthplace in UK data (2) meeting QC and (3) within 1eps PCA
    df_geo = df[(df['birth_east_coord'].notna()) & (df['birth_UKelsewhere']!='Elsewhere') & (df['used_in_pca']==1) & (df['within_1epsilon_pca']==True)]

    # histogram
    x = df_geo['birth_east_coord']
    y = df_geo['birth_north_coord']
    counts, xedges, yedges = np.histogram2d(x,y,bins=20)

    # calc distances
    center_bin_1 = [12,6]
    center_bin_2 = [9,9]
    center_bin_3 = [15,2]
    center_coord_1 = get_bin_center(center_bin_1,xedges,yedges)
    center_coord_2 = get_bin_center(center_bin_2,xedges,yedges)
    center_coord_3 = get_bin_center(center_bin_3,xedges,yedges)
    df_geo['dist_center_1']= np.sqrt(((df_geo.iloc[:,3:5] - center_coord_1) ** 2).sum(axis=1))
    df_geo['dist_center_2']= np.sqrt(((df_geo.iloc[:,3:5] - center_coord_2) ** 2).sum(axis=1))
    df_geo['dist_center_3']= np.sqrt(((df_geo.iloc[:,3:5] - center_coord_3) ** 2).sum(axis=1))

    # binned counts
    df_geo['binned_counts'] = df_geo.apply(lambda row: get_count(row['birth_east_coord'],row['birth_north_coord'],counts,xedges,yedges),axis=1)

    # filter
    thresh=15
    df_geo_filt = df_geo[df_geo['binned_counts']>thresh]
    N=len(df_geo_filt)

    # weights
    df_geo_filt['freq_binned']=df_geo_filt['binned_counts']/N
    if label=='SIRgeouniform':
        df_geo_filt['weights']=1/df_geo_filt['freq_binned']
    elif label=='SIRcenter1geo10000':
        df_geo_filt['weights']=df_geo_filt.apply(lambda row: gaussian_2d(row['dist_center_1'],10000),axis=1)/df_geo_filt['freq_binned']
    elif label=='SIRcenter1geo50000':
        df_geo_filt['weights']=df_geo_filt.apply(lambda row: gaussian_2d(row['dist_center_1'],50000),axis=1)/df_geo_filt['freq_binned']
    elif label=='SIRcenter1geo100000':
        df_geo_filt['weights']=df_geo_filt.apply(lambda row: gaussian_2d(row['dist_center_1'],100000),axis=1)/df_geo_filt['freq_binned']
    elif label=='SIRcenter2geo10000':
        df_geo_filt['weights']=df_geo_filt.apply(lambda row: gaussian_2d(row['dist_center_2'],10000),axis=1)/df_geo_filt['freq_binned']
    elif label=='SIRcenter2geo50000':
        df_geo_filt['weights']=df_geo_filt.apply(lambda row: gaussian_2d(row['dist_center_2'],50000),axis=1)/df_geo_filt['freq_binned']
    elif label=='SIRcenter2geo100000':
        df_geo_filt['weights']=df_geo_filt.apply(lambda row: gaussian_2d(row['dist_center_2'],100000),axis=1)/df_geo_filt['freq_binned']
    elif label=='SIRcenter3geo10000':
        df_geo_filt['weights']=df_geo_filt.apply(lambda row: gaussian_2d(row['dist_center_3'],10000),axis=1)/df_geo_filt['freq_binned']
    elif label=='SIRcenter3geo50000':
        df_geo_filt['weights']=df_geo_filt.apply(lambda row: gaussian_2d(row['dist_center_3'],50000),axis=1)/df_geo_filt['freq_binned']
    elif label=='SIRcenter3geo100000':
        df_geo_filt['weights']=df_geo_filt.apply(lambda row: gaussian_2d(row['dist_center_3'],100000),axis=1)/df_geo_filt['freq_binned']
    else:
        raise ValueError(f"Unexpected label: {label}")

    sampled_df = df_geo_filt.sample(n=snakemake.params.numsamp,replace=True,weights='weights',random_state=snakemake.params.seed)

elif 'pca' in label:
    # subsample to individuals meeting QC
    df_pca = df[(df['used_in_pca'] == 1)]
    # histogram
    x_pca = df_pca['PC1']
    y_pca = df_pca['PC2']
    counts_pca, xedges_pca, yedges_pca, = np.histogram2d(x_pca, y_pca, bins=20)

    # distances
    center_bin_1 = [19, 4]
    center_bin_2 = [0, 0]
    center_bin_3 = [12, 19]
    center_coord_1 = get_bin_center(center_bin_1, xedges_pca, yedges_pca)
    center_coord_2 = get_bin_center(center_bin_2, xedges_pca, yedges_pca)
    center_coord_3 = get_bin_center(center_bin_3, xedges_pca, yedges_pca)
    df_pca['dist_center_1'] = np.sqrt(((df_pca.iloc[:, 7:9] - center_coord_1) ** 2).sum(axis=1))
    df_pca['dist_center_2'] = np.sqrt(((df_pca.iloc[:, 7:9] - center_coord_2) ** 2).sum(axis=1))
    df_pca['dist_center_3'] = np.sqrt(((df_pca.iloc[:, 7:9] - center_coord_3) ** 2).sum(axis=1))

    # binned counts
    df_pca['binned_counts'] = df_pca.apply(
        lambda row: get_count(row['PC1'], row['PC2'], counts_pca, xedges_pca, yedges_pca), axis=1)
    # filter
    thresh = 15
    df_pca_filt = df_pca[df_pca['binned_counts'] > thresh]
    N = len(df_pca_filt)

    # weights
    df_pca_filt['freq_binned'] = df_pca_filt['binned_counts'] / N

    if label=='SIRpcauniform':
        df_pca_filt['weights']=1/df_pca_filt['freq_binned']
    elif label=='SIRcenter1pca0.0005':
        df_pca_filt['weights'] = df_pca_filt.apply(lambda row: gaussian_2d(row['dist_center_1'], 0.0005), axis=1) / \
                                 df_pca_filt['freq_binned']
    elif label=='SIRcenter1pca0.0025':
        df_pca_filt['weights'] = df_pca_filt.apply(lambda row: gaussian_2d(row['dist_center_1'], 0.0025), axis=1) / \
                                 df_pca_filt['freq_binned']
    elif label=='SIRcenter1pca0.005':
        df_pca_filt['weights'] = df_pca_filt.apply(lambda row: gaussian_2d(row['dist_center_1'], 0.005), axis=1) / \
                                 df_pca_filt['freq_binned']
    elif label=='SIRcenter2pca0.0005':
        df_pca_filt['weights'] = df_pca_filt.apply(lambda row: gaussian_2d(row['dist_center_2'], 0.0005), axis=1) / \
                                 df_pca_filt['freq_binned']
    elif label=='SIRcenter2pca0.0025':
        df_pca_filt['weights'] = df_pca_filt.apply(lambda row: gaussian_2d(row['dist_center_2'], 0.0025), axis=1) / \
                                 df_pca_filt['freq_binned']
    elif label=='SIRcenter2pca0.005':
        df_pca_filt['weights'] = df_pca_filt.apply(lambda row: gaussian_2d(row['dist_center_2'], 0.005), axis=1) / \
                                 df_pca_filt['freq_binned']
    elif label=='SIRcenter3pca0.0005':
        df_pca_filt['weights'] = df_pca_filt.apply(lambda row: gaussian_2d(row['dist_center_3'], 0.0005), axis=1) / \
                                 df_pca_filt['freq_binned']
    elif label=='SIRcenter3pca0.0025':
        df_pca_filt['weights'] = df_pca_filt.apply(lambda row: gaussian_2d(row['dist_center_3'], 0.0025), axis=1) / \
                                 df_pca_filt['freq_binned']
    elif label=='SIRcenter3pca0.005':
        df_pca_filt['weights'] = df_pca_filt.apply(lambda row: gaussian_2d(row['dist_center_3'], 0.005), axis=1) / \
                                 df_pca_filt['freq_binned']
    else:
        raise ValueError(f"Unexpected label: {label}")

    sampled_df = df_pca_filt.sample(n=snakemake.params.numsamp, replace=True, weights='weights',
                                    random_state=snakemake.params.seed)

else:
    raise ValueError(f"Unexpected label: {label}")

new_df = pd.DataFrame({"X1": sampled_df["id"], "X2": sampled_df["id"], "Clust": label})
new_df.to_csv(snakemake.output[0], sep="\t", index=False)

