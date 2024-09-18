import numpy as np
from matplotlib import pyplot as plt
from scipy.stats import binom,nbinom,gamma
from scipy.special import exp1
from matplotlib import rcParams
import glob
import os
import json
from matplotlib.colors import LogNorm
from matplotlib.colors import Normalize

def get_theta(ep,ep2):
    var = ep2-ep**2
    return ep**2/var

def get_gamma(ep,ep2):
    var = ep2 - ep ** 2
    return ep/var

def get_lc_squared(sigma, s):
    return sigma ** 2 / s

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

def print_moments(files):
    combined_sampled_p, combined_zero_samples = concatenate_data(files)
    f = files[0]
    with open(f,"r") as json_file:
        data = json.load(json_file)
    s = data["s"]
    mu = data["mu"]
    rho = data["rho"]
    sigma = data["sigma"]
    w = data["w"]

    # print(f"File path: {f}")
    # print(f"s: {s}")
    # print(f"mu: {mu}")
    # print(f"rho: {rho}")
    # print(f"sigma: {sigma}")
    # print(f"w: {w}\n")

    # first moment
    sim_EP = get_EP_sim(combined_sampled_p, combined_zero_samples)
    theory_EP = get_EP_theory(mu, s)

    # print(f"EP (sim): {sim_EP}")
    # print(f"EP (theory): {theory_EP}\n")

    # second moment
    sim_EP2 = get_EPsquared_sim(combined_sampled_p, combined_zero_samples)
    theory_EP2 = get_EPsquared_theory(mu=mu, s=s, rho=rho, sigma=sigma, w=w)

    return sim_EP, theory_EP, sim_EP2, theory_EP2

def main():
    rcParams = {'font.size': 14}
    plt.rcParams.update(rcParams)

    s_list = [0.1,0.01]
    mu = 1e-09
    dens = [2.5,20]
    L_list = [1000,10000]#[1000]
    sigma = 10
    n_list = [10000]#[1000,10000,100000]
    max_x = 100
    w_list = np.logspace(1, 4, 10)[0:7]

    colors = ['steelblue','darkorange','mediumseagreen']

    for n in n_list:
        for rho in dens:
            for L in L_list:
                for s in s_list:
                    fig, ax = plt.subplots(1,2,figsize=(16, 8))
                    i=0
                    if L == 1000:
                        w_list_subset = [w_list[2], w_list[3], w_list[4]]
                    elif L == 10000:
                        w_list_subset = [w_list[3], w_list[-2], w_list[-1]]
                    for w in w_list_subset:
                        file_path = f"results/20240521/s{s}_mu{mu}_rho{rho}_L{L}_sigma{sigma}_time1000000.0_r0.1_burnin_wrapped_norm_gaussian_w{w}*"
                        files = glob.glob(os.path.join(file_path))
                        if len(files) == 10:
                            combined_sampled_p, combined_zero_samples = concatenate_data(files)
                            combined_sampled_p = np.where(combined_sampled_p > 1e-100, combined_sampled_p, 0.0)
                            zeros_array = np.zeros(combined_zero_samples)
                            ps = np.concatenate([combined_sampled_p,zeros_array])
                            # print(ps)
                            ps_no_zero = ps[ps != 0]
                            # log_bins = np.logspace(np.log10(min(ps_no_zero)), np.log10(max(ps)), 50)
                            log_bins = np.logspace(-55,-3,1000)
                            # print(log_bins)

                            ax[0].hist(ps,bins=log_bins,color=colors[i],label=f'w={round(w,1)}',alpha=0.6)
                            ax[0].set_xscale('log')

                            ax[1].hist(ps, bins=log_bins, color=colors[i], label=f'w={round(w, 1)}', alpha=0.6)
                            ax[1].set_xscale('log')
                            ax[1].set_xlim(1e-7,1e-2)

                            # print(np.mean(ps))
                            # print(mu/s)

                            # sfs = sample_sfs(combined_sampled_p, combined_zero_samples, n, max_x)
                            # ax.loglog(np.arange(0, max_x), sfs[0:max_x], label=f'w={round(w,1)} (simulation)', marker='X',
                            #           linestyle='',
                            #           markersize=8, color=colors[i])
                            #
                            # fit_pmf = fitted_nb_pmf(sfs)
                            # ax.loglog(np.arange(0,max_x), fit_pmf, label=f'w={round(w,1)} (NB fit)', marker=None, linestyle='--',
                            #           linewidth=3, alpha=1, color=colors[i])
                            # nb_dist = [get_sfs_theory(y, n, mu, s, rho, sigma, w) for y in np.arange(0, max_x)]
                            # ax.loglog(np.arange(0, max_x), nb_dist, label=f'w={round(w,1)} (theory)', marker=None, linestyle='-',
                            #           linewidth=3, alpha=0.8, color=colors[i])
                            i+=1
                        else:
                            continue
                    # ax.set_ylim(1e-16, 1e0)
                    # ax.set_xlabel('allele count')
                    # ax.set_ylabel('expected proportion of sites')
                    plt.suptitle(f"s={s}, L={L}, rho={rho}, l_c={round(np.sqrt(get_lc_squared(sigma,s)),2)}")
                    ax[0].legend()
                    ax[1].legend()
                    ax[0].axvline(mu/s,linestyle='--',color='black',label=r'$mu/s$')
                    plt.tight_layout()
                    # plt.show()
                    plt.savefig(f"plots_20240604/s{s}_sigma{sigma}_L{L}_rho{rho}_dist_p.png")
                    plt.close()


if __name__ == '__main__':
    main()