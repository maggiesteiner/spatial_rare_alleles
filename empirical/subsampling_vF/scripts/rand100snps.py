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

snps_pca = priv_pca.sample(n=100,random_state=10)['SNP'].tolist()
snps_geo = priv_geo.sample(n=100,random_state=10)['SNP'].tolist()

overall_df = pd.DataFrame(columns=['snp','scenario', 'seed', 'freq'])

for sc in scen:
    for se in seeds:
        for snp in snps_pca:
            print(f'SNP: {snp}')
            chrom = snp.split(':')[0]
            filepath = f'data/frq/chr{chrom}_{sc}_seed{se}_n10000.frq.strat'
            print(f'Processing file {filepath}')
            
            df = pd.read_csv(filepath, delim_whitespace=True)
            df['MAF'] = pd.to_numeric(df['MAF'], errors='coerce')

            freq = df.loc[df['SNP'] == snp, 'MAF'].iloc[0]
            print(f'Frequency: {freq}')
            overall_df = overall_df.append({
                'snp': snp,
                'scenario': sc,
                'seed': se,
                'freq': freq,
            }, ignore_index=True)


overall_df.to_csv('random_100_snps_freq_MAC_10_20.csv', index=False)
