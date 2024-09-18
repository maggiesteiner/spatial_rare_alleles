import numpy as np
from scipy.special import expi, exp1, factorial
import argparse
import os
import json
import re
import glob
from matplotlib import pyplot as plt
from matplotlib import rcParams

def get_lc_squared(sigma, s):
    return sigma ** 2 / s

def u2_exact(wt):
    u2 = (-1 * np.exp(wt ** 2 / (4 * np.pi ** 2)) / (8 * np.pi)) * expi(-1 * wt ** 2 / (4 * np.pi ** 2))
    return u2

def get_EPsquared_theory(mu, s, rho, sigma, w):
    lcs = get_lc_squared(sigma, s)
    term = (w / np.sqrt(lcs)) ** 2
    if term <= 800:
        prod_term = np.exp(term) * exp1(term)
    else:
        prod_term = sum((factorial(k) / (-term)**k for k in range(7))) / term
    return (mu / (s ** 2 * rho * 4 * np.pi * lcs)) * prod_term + mu ** 2 / s ** 2
def get_EPsquared_theory_swap(mu, s, rho, sigma, w):
    lcs = get_lc_squared(w, s)
    u2 = u2_exact(sigma / np.sqrt(lcs))
    return ((2 * mu / (s ** 2 * rho * lcs)) * u2 + mu ** 2 / (s ** 2))

def get_EP_theory(mu, s):
    return mu / s

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

    print(f"File path: {f}")
    print(f"s: {s}")
    print(f"mu: {mu}")
    print(f"rho: {rho}")
    print(f"sigma: {sigma}")
    print(f"w: {w}\n")

    # first moment
    sim_EP = get_EP_sim(combined_sampled_p, combined_zero_samples)
    theory_EP = get_EP_theory(mu, s)

    print(f"EP (sim): {sim_EP}")
    print(f"EP (theory): {theory_EP}\n")

    # second moment
    sim_EP2 = get_EPsquared_sim(combined_sampled_p, combined_zero_samples)
    theory_EP2 = get_EPsquared_theory(mu=mu, s=s, rho=rho, sigma=sigma, w=w)

    return sim_EP, theory_EP, sim_EP2, theory_EP2
    # print(f"EP2 (sim): {sim_EP2}")
    # print(f"EP2 (theory finite habitat): {theory_EP2}")
    # print("-----------------------\n")

def get_theta(ep,ep2):
    var = ep2-ep**2
    return ep**2/var

def get_gamma(ep,ep2):
    var = ep2 - ep ** 2
    return ep/var

def main():

    rcParams.update({'font.size': 14})

    s_list = [0.1,0.01]
    mu = 1e-09
    dens = [2.5,20]
    L_list = [1000,10000]
    sigma = 10
    w_list = np.logspace(1,4,10)

    for rho in dens:
        for L in L_list:
            for s in s_list:
                fig, axs = plt.subplots(1,2,figsize=(10,6))
                sim_theta_list = []
                sim_gamma_list = []
                theory_theta_list = []
                theory_gamma_list = []
                for w in w_list:
                    file_path = f"results/20240521/s{s}_mu{mu}_rho{rho}_L{L}_sigma{sigma}_time1000000.0_r0.1_burnin_wrapped_norm_gaussian_w{w}*"
                    files = glob.glob(os.path.join(file_path))
                    if len(files) == 10:
                        sim_EP, theory_EP, sim_EP2, theory_EP2 = print_moments(files)
                        sim_theta_list.append(get_theta(sim_EP,sim_EP2))
                        theory_theta_list.append(get_theta(theory_EP,theory_EP2))
                        sim_gamma_list.append(get_gamma(sim_EP,sim_EP2))
                        theory_gamma_list.append(get_gamma(theory_EP,theory_EP2))
                    else:
                        continue

                if len(sim_gamma_list) != len(w_list):
                    continue

                axs[0].loglog(w_list,theory_theta_list,label='Gaussian (theory)',linestyle='',marker='^',markersize=8,color='black')
                axs[0].loglog(w_list, sim_theta_list, label='Gaussian (simulation)', marker='x', linestyle='',color='red',markersize=8)
                axs[0].set_title(r"Effective mutation supply ($\theta_E$)")
                axs[0].set_xlabel(r"sampling width ($w$)")
                axs[0].set_ylabel(r"$\theta_E$")
                axs[0].legend(frameon=True,loc='best',title='Sampling kernel')
                axs[1].loglog(w_list, theory_gamma_list, label='Gaussian (theory)',linestyle='',marker='^',markersize=8,color='black')
                axs[1].loglog(w_list, sim_gamma_list, label='Gaussian (simulation)', marker='x', linestyle='',color='red',markersize=8)
                axs[1].set_title(r"Effective selection intensity ($\gamma_E$)")
                axs[1].set_xlabel(r"sampling width ($w$)")
                axs[1].set_ylabel(r"$\gamma_E$")

                axs[0].axvline(np.sqrt(get_lc_squared(sigma, s)), color='gray', linestyle='--')
                axs[1].axvline(np.sqrt(get_lc_squared(sigma, s)), color='gray', linestyle='--')

                # axs[1].legend()
                fig.suptitle(f"L={L}, density={rho}, s={s}, mu={mu}, sigma={sigma}")
                plt.tight_layout()
                plt.savefig(f"plots_20240604/effective_params_over_w_s{s}_rho{rho}_L{L}_sigma{sigma}.png")
            # plt.show()



if __name__ == '__main__':
    main()
