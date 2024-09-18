import numpy as np
import pandas as pd

chroms = np.arange(1,23)
scen = ['uniform','5epsilonpca','10epsilonpca','50epsilonpca','100epsilonpca',
        'closest10kpca','5epsilongeo','10epsilongeo','50epsilongeo','closest10kgeo','uniformgeo']
        
seeds = np.arange(1,11)

priv_pca = pd.read_csv('private_vars_pca_geq10_leq20.csv')
priv_geo = pd.read_csv('private_vars_geo_geq10_leq20.csv')

c_list = []
sc_list = []
se_list = []
priv_pca_list = []
priv_geo_list = []

for c in chroms:
    for sc in scen:
        for se in seeds:
            filepath = f'data/frq/chr{c}_{sc}_seed{se}_n10000.frq.strat'
            print(f'Processing file {filepath}')
            df = pd.read_csv(filepath,delim_whitespace=True)
            df['MAC'] = pd.to_numeric(df['MAC'],errors='coerce')
            
            c_list.append(c)
            sc_list.append(sc)
            se_list.append(se)
            
            if 'geo' in sc:
                priv_geo_temp = pd.merge(df,priv_geo,on='SNP',how='inner')
                priv_geo_list.append((priv_geo_temp['MAC']!=0).sum())
                priv_pca_list.append(np.nan)
                print((priv_geo_temp['MAC']!=0).sum())
            else:
                priv_pca_temp = pd.merge(df,priv_pca,on='SNP',how='inner')
                priv_pca_list.append((priv_pca_temp['MAC']!=0).sum())
                priv_geo_list.append(np.nan)
                print((priv_pca_temp['MAC']!=0).sum())
            
            
            
res = pd.DataFrame({
    'chrom':c_list,
    'scenario':sc_list,
    'seed':se_list,
    'priv_pca_count':priv_pca_list,
    'priv_geo_count':priv_geo_list
    })
      
res.to_csv('private_var_counts_MAC_10_20.csv',index=False)      
            
            
