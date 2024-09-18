import numpy as np
import pandas as pd

chroms = np.arange(1, 23)
scen = ['uniform', '5epsilonpca', '10epsilonpca', '50epsilonpca', '100epsilonpca',
        'closest10kpca', '5epsilongeo', '10epsilongeo', '50epsilongeo', 'closest10kgeo', 'uniformgeo']

seeds = np.arange(1, 11)

priv_pca = pd.read_csv('private_vars_pca_geq10_leq20.csv')
priv_geo = pd.read_csv('private_vars_geo_geq10_leq20.csv')

sc_list = []
se_list = []
priv_pca_list = []
priv_geo_list = []


overall_df = pd.DataFrame(columns=['scenario', 'seed', 'priv_pca_mean', 'priv_geo_mean'])

for sc in scen:
    for se in seeds:
        priv_pca_sum = 0
        priv_geo_sum = 0
        priv_pca_count = 0
        priv_geo_count = 0
        
        for c in chroms:
            filepath = f'data/frq/chr{c}_{sc}_seed{se}_n10000.frq.strat'
            print(f'Processing file {filepath}')
            df = pd.read_csv(filepath, delim_whitespace=True)
            df['MAF'] = pd.to_numeric(df['MAF'], errors='coerce')
            
            if 'geo' in sc:
                priv_geo_temp = pd.merge(df, priv_geo, on='SNP', how='inner')
                priv_geo_sum += priv_geo_temp['MAF'].sum()
                priv_geo_count += priv_geo_temp['MAF'].count()
            else:
                priv_pca_temp = pd.merge(df, priv_pca, on='SNP', how='inner')
                priv_pca_sum += priv_pca_temp['MAF'].sum()
                priv_pca_count += priv_pca_temp['MAF'].count()
        
        priv_pca_mean = priv_pca_sum / priv_pca_count if priv_pca_count > 0 else np.nan
        priv_geo_mean = priv_geo_sum / priv_geo_count if priv_geo_count > 0 else np.nan
        print(priv_pca_mean)
        print(priv_geo_mean)        
        overall_df = overall_df.append({
            'scenario': sc,
            'seed': se,
            'priv_pca_mean': priv_pca_mean,
            'priv_geo_mean': priv_geo_mean
        }, ignore_index=True)


overall_df.to_csv('private_var_mean_freq_MAC_10_20.csv', index=False)

