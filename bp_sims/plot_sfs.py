import numpy as np
from matplotlib import pyplot as plt
from scipy.stats import binom,nbinom
from scipy.special import exp1
from matplotlib import rcParams
import glob
import os
import json
from matplotlib.colors import LogNorm
from matplotlib.colors import Normalize

def get_lc_squared(sigma, s):
    return sigma ** 2 / s
# def sample_sfs(ps,zeros,n,max_allele_count=10000):
#     running_sfs = np.zeros(max_allele_count + 1)
#     j = np.arange(max_allele_count)
#     for p in ps:
#         # print(n)
#         # print(p)
#         p=p if p>1e-100 else 0.0
#         running_sfs[:-1] += binom.pmf(j, n, p)  # pmf, entries 0 through max_allele_count-1
#         running_sfs[-1] += binom.sf(max_allele_count - 1, n, p)  # 1 - cdf
#     running_sfs[0] += zeros
#     expected_sfs = running_sfs/np.sum(running_sfs)
#     return expected_sfs

def sample_sfs(ps, zeros, n, max_allele_count=10):
    running_sfs = np.zeros(max_allele_count + 1)
    ps = np.where(ps > 1e-100, ps, 0.0)
    j = np.arange(max_allele_count)
    pmf_matrix = binom.pmf(j[:, np.newaxis], n, ps)
    running_sfs[:-1] = np.sum(pmf_matrix, axis=1)
    running_sfs[-1] = np.sum(binom.sf(max_allele_count - 1, n, ps))
    running_sfs[0] += zeros
    expected_sfs = running_sfs / np.sum(running_sfs)
    return expected_sfs

def get_EP_theory(mu, s):
    return mu / s

def get_EPsquared_theory(mu, s, rho, sigma, w):
    lcs = get_lc_squared(sigma, s)
    return (mu/(s**2*rho*4*np.pi*lcs))*np.exp((w/np.sqrt(lcs))**2)*exp1((w/np.sqrt(lcs))**2)+mu**2/s**2

def get_sfs_theory(x,n,mu,s,rho,sigma,w):
    mean = get_EP_theory(mu,s)
    var = get_EPsquared_theory(mu,s,rho,sigma,w) - mean**2
    alpha = mean**2/var
    beta = mean/var
    return nbinom.pmf(x,alpha,beta/(beta+n))

def get_EPsquared_sim(ps, zeros):
    return (np.sum(ps ** 2) / (len(ps) + zeros))

def get_EP_sim(ps, zeros):
    return np.sum(ps) / (len(ps) + zeros)

def load_data(file):
    with open(file, "r") as json_file:
        data = json.load(json_file)
    return data
def concatenate_data(files):
    combined_sampled_p = []
    combined_zero_samples = 0
    for f in files:
        data = load_data(f)
        combined_sampled_p.extend(data['sampled_p_flattened'])
        combined_zero_samples += data['zero_samples']
    return np.array(combined_sampled_p), combined_zero_samples

def main():

    rcParams.update({'font.size': 14})

    s = 0.1
    mu = 1e-09
    dens = [2.5,20]
    L_list = [10000]
    sigma = 10
    n_list = [1000]
    max_x = 100
    w_list = np.logspace(1,4,10)

    cmap = plt.get_cmap('viridis')
    norm = LogNorm(vmin=min(w_list), vmax=max(w_list))

    for n in n_list:
        for rho in dens:
            for L in L_list:
                fig, ax = plt.subplots()
                for w in w_list[0:2]:
                    file_path = f"results/20240521/s{s}_mu{mu}_rho{rho}_L{L}_sigma{sigma}_time1000000.0_r0.1_burnin_wrapped_norm_gaussian_w{w}*"
                    files = glob.glob(os.path.join(file_path))
                    combined_sampled_p, combined_zero_samples = concatenate_data(files)
                    sfs = sample_sfs(combined_sampled_p, combined_zero_samples, n, max_x)


                    ax.loglog(np.arange(0, max_x), sfs[0:max_x], label=f'w={w} (simulation)', marker='^',linestyle='--',color=cmap(norm(w)))
                    ax.set_ylim(1e-16, 1e0)
                    ax.set_xlabel('allele count')
                    ax.set_ylabel('expected proportion of sites')

                    nb_dist = [get_sfs_theory(y,n,mu,s,rho,sigma,w) for y in np.arange(0, max_x)]
                    ax.loglog(np.arange(0, max_x), nb_dist, label=f'w={w} (theory)', marker='o',linestyle='-',color=cmap(norm(w)))
                # ax.legend()
                cbar = plt.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax)
                cbar.set_label(r'$w$')
                plt.title(f"n={n}, L={L}, density={rho}")
                plt.tight_layout()
                # plt.savefig(f"plots_20240522/sfs_rho{rho}_L{L}_sigma{sigma}_s{s}_mu{mu}.png")
                plt.show()



if __name__ == '__main__':
    main()
