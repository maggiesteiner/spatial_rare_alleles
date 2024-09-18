import numpy as np
from matplotlib import pyplot as plt
from scipy.stats import binom,nbinom
from scipy.special import exp1, factorial
from matplotlib import rcParams
import glob
import os
import json
from matplotlib.colors import LogNorm
from matplotlib.colors import Normalize

def get_lc_squared(sigma, s):
    return sigma ** 2 / s
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
    term = (w / np.sqrt(lcs)) ** 2
    if term <= 800:
        prod_term = np.exp(term) * exp1(term)
    else:
        prod_term = sum((factorial(k) / (-term)**k for k in range(7))) / term
    return (mu / (s ** 2 * rho * 4 * np.pi * lcs)) * prod_term + mu ** 2 / s ** 2

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
    rcParams = {'font.size': 14}
    plt.rcParams.update(rcParams)

    s_list = [0.1,0.01]
    mu = 1e-09
    dens = [2.5,20]
    L_list = [10000]#[1000]
    sigma = 10
    n_list = [10000]#[1000,10000,100000]
    max_x = 100
    w_list = np.logspace(1, 4, 10)[0:7]

    # cmap = plt.get_cmap('viridis')
    # norm = LogNorm(vmin=min(w_list), vmax=max(w_list))
    colors = ['steelblue','darkorange','mediumseagreen']

    for n in n_list:
        for rho in dens:
            for L in L_list:
                for s in s_list:
                    fig, ax = plt.subplots(figsize=(10, 8))
                    i=0
                    if L == 1000:
                        w_list_subset = [w_list[2], w_list[3], w_list[4]]
                    elif L == 10000:
                        w_list_subset = [w_list[3], w_list[-1]]
                    for w in w_list_subset:
                        file_path = f"results/20240521/s{s}_mu{mu}_rho{rho}_L{L}_sigma{sigma}_time1000000.0_r0.1_burnin_wrapped_norm_gaussian_w{w}*"
                        files = glob.glob(os.path.join(file_path))
                        if len(files) == 10:
                            combined_sampled_p, combined_zero_samples = concatenate_data(files)
                            sfs = sample_sfs(combined_sampled_p, combined_zero_samples, n, max_x)
                            ax.loglog(np.arange(0, max_x), sfs[0:max_x], label=f'w={round(w,1)} (simulation)', marker='X',
                                      linestyle='',
                                      markersize=8, color=colors[i])
                            nb_dist = [get_sfs_theory(y, n, mu, s, rho, sigma, w) for y in np.arange(0, max_x)]
                            ax.loglog(np.arange(0, max_x), nb_dist, label=f'w={round(w,1)} (theory)', marker=None, linestyle='-',
                                      linewidth=3, alpha=0.8, color=colors[i])
                            i+=1
                        else:
                            continue
                    ax.set_ylim(1e-16, 1e0)
                    ax.set_xlabel('allele count')
                    ax.set_ylabel('expected proportion of sites')
                    ax.set_title(f"n={n}, s={s}, L={L}, rho={rho}, l_c={round(np.sqrt(get_lc_squared(sigma,s)),2)}")
                    ax.legend()
                    # cbar = plt.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax)
                    # cbar.set_label(r'$w$')
                    plt.tight_layout()
                    plt.savefig(f"plots_20240604/n{n}_s{s}_sigma{sigma}_L{L}_rho{rho}_sfs_3vals.png")
                    plt.close()


if __name__ == '__main__':
    main()