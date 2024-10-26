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

def parse_center(s):
    e_index = s.index('E')
    n_index = s.index('N')
    locE = int(s[e_index + 1:n_index])
    locN = int(s[n_index + 1:])
    return [locE,locN]

df = pd.read_csv(snakemake.input[0])

df_geo = df[(df['birth_east_coord'].notna()) & (df['birth_UKelsewhere']!='Elsewhere') & (df['used_in_pca']==1) & (df['within_1epsilon_pca']==True)]
# histogram
x = df_geo['birth_east_coord']
y = df_geo['birth_north_coord']
counts, xedges, yedges = np.histogram2d(x,y,bins=20)

# read in snakemake params
centers = snakemake.params.center_list#['center1','center2','center3']
w_list = snakemake.params.w_list#[25000,50000,100000]

# distances
center_bin_1 = parse_center(centers[0])
center_bin_2 = parse_center(centers[1])
center_bin_3 = parse_center(centers[2])
center_coord_1 = get_bin_center(center_bin_1,xedges,yedges)
center_coord_2 = get_bin_center(center_bin_2,xedges,yedges)
center_coord_3 = get_bin_center(center_bin_3,xedges,yedges)
df_geo[centers[0]]= np.sqrt(((df_geo.iloc[:,3:5] - center_coord_1) ** 2).sum(axis=1))
df_geo[centers[1]]= np.sqrt(((df_geo.iloc[:,3:5] - center_coord_2) ** 2).sum(axis=1))
df_geo[centers[2]]= np.sqrt(((df_geo.iloc[:,3:5] - center_coord_3) ** 2).sum(axis=1))

# binned counts
df_geo['binned_counts'] = df_geo.apply(lambda row: get_count(row['birth_east_coord'],row['birth_north_coord'],counts,xedges,yedges),axis=1)

# filter
thresh=15
df_geo_filt = df_geo[df_geo['binned_counts']>thresh]
N=len(df_geo_filt)

# weights
df_geo_filt['freq_binned']=df_geo_filt['binned_counts']/N

# uniform weights
df_geo_filt['uniformgeo'] = 1/df_geo_filt['freq_binned']
# df_geo_filt['uniformgeo'] = df_geo_filt['uniformgeo']*scale_factor
total_weight = df_geo_filt['uniformgeo'].sum()
df_geo_filt['uniformgeo'] /= total_weight
df_geo_filt['uniformgeo'] = df_geo_filt['uniformgeo'].apply(lambda x: 0 if abs(x) < 1e-100 else x)
# re-normalize
total_weight = df_geo_filt['uniformgeo'].sum()
df_geo_filt['uniformgeo'] /= total_weight

unif_weights = pd.DataFrame({
            'id1': df_geo_filt['id'].astype(int),
            'id2': df_geo_filt['id'].astype(int),
            'weights': df_geo_filt['uniformgeo']
        })
#pd.DataFrame([df_geo_filt['id'].astype(int),df_geo_filt['id'].astype(int),df_geo_filt['uniformgeo']]).T
unif_weights.to_csv('data/weights/uniformgeo.weights', sep=" ", index=False, header=False)

# centers = ['center1','center2','center3']
# w_list = [25000,50000,100000]

for c in centers:
    for w in w_list:
        colname = f'{c}geo{w}'
        df_geo_filt[colname] = df_geo_filt.apply(lambda row: gaussian_2d(row[c],w),axis=1)/df_geo_filt['freq_binned']
        # df_geo_filt[colname] = df_geo_filt[colname] * scale_factor
        total_weight = df_geo_filt[colname].sum()
        df_geo_filt[colname] /= total_weight
        df_geo_filt[colname] = df_geo_filt[colname].apply(lambda x: 0 if abs(x) < 1e-100 else x)
        # re-normalize
        total_weight = df_geo_filt[colname].sum()
        df_geo_filt[colname] /= total_weight
        # total_weight = df_geo_filt[colname].sum()
        # df_geo_filt[colname] /= total_weight
        # print sum as check
        # tot = df_geo_filt[colname].sum()
        # print(f'abs(total-1) = {np.abs(tot-1.0)}')
        # generate data frame
        temp = pd.DataFrame({
            'id1': df_geo_filt['id'].astype(int),
            'id2': df_geo_filt['id'].astype(int),
            'weights': df_geo_filt[colname]
        })
        temp.to_csv(f'data/weights/{colname}.weights', sep=" ", index=False, header=False)


