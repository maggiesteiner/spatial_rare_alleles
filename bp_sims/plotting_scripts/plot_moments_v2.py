import numpy as np
from scipy.special import expi, exp1, factorial
import argparse
import os
import json
import re
import glob
from matplotlib import pyplot as plt
from matplotlib import rcParams
import warnings
plt.rcParams['font.size'] = 18
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

def get_EPsquared_theory_approx(mu, s, rho, sigma, w):
    lcs = get_lc_squared(sigma, s)
    term = (w / np.sqrt(lcs)) ** 2
    prod_term = sum((factorial(k) / (-term)**k for k in range(1))) / term
    return (mu / (s ** 2 * rho * 4 * np.pi * lcs)) * prod_term + mu ** 2 / s ** 2

def get_EP_theory(mu, s):
    return mu / s

def get_EPsquared_sim(ps, zeros):
    return (np.sum(ps ** 2) / (len(ps) + zeros))

def get_EP_sim(ps, zeros):
    return np.sum(ps) / (len(ps) + zeros)

def load_data(file):
    try:
        with open(file, "r") as json_file:
            data = json.load(json_file)
        return data
    except FileNotFoundError:
        print(f"File not found: {file}")
        return None

def concatenate_data(files):
    combined_sampled_p = []
    combined_zero_samples = 0
    for f in files:
        data = load_data(f)
        if data is not None:
            combined_sampled_p.extend(data['sampled_p_flattened'])
            combined_zero_samples += data['zero_samples']
    return np.array(combined_sampled_p), combined_zero_samples

def print_moments(files):
    combined_sampled_p, combined_zero_samples = concatenate_data(files)
    if len(combined_sampled_p) == 0:
        print("No valid data found in files.")
        return None, None, None, None

    f = files[0]
    data = load_data(f)
    if data is None:
        return None, None, None, None

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
    theory_EP2_approx = get_EPsquared_theory_approx(mu=mu, s=s, rho=rho, sigma=sigma, w=w)


    return sim_EP, theory_EP, sim_EP2, theory_EP2, theory_EP2_approx

def main():

    rcParams.update({'font.size': 14})

    s_list = [0.1, 0.01, 0.001]#,0.01]
    mu = 1e-09
    dens = [20]#, 20]
    L_list = [1000]#, 10000]
    sigma = 10
    w_list = np.logspace(1, 4, 10)

    for rho in dens:
        for L in L_list:
            fig, axs = plt.subplots(3, 2, figsize=(12, 18))
            for i,s in enumerate(s_list):

                sim_ep_list = []
                sim_ep2_list = []
                theory_ep_list = []
                theory_ep2_list = []
                theory_ep2_list_approx = []
                for w in w_list:
                    file_path = f"results/20240521/s{s}_mu{mu}_rho{rho}_L{L}_sigma{sigma}_time1000000.0_r0.1_burnin_wrapped_norm_gaussian_w{w}*"
                    files = glob.glob(os.path.join(file_path))
                    if len(files)==10:
                        sim_EP, theory_EP, sim_EP2, theory_EP2, theory_EP2_approx = print_moments(files)
                        if sim_EP is not None:
                            sim_ep_list.append(sim_EP)
                            sim_ep2_list.append(sim_EP2)
                            theory_ep_list.append(theory_EP)
                            theory_ep2_list.append(theory_EP2)
                            theory_ep2_list_approx.append(theory_EP2_approx)
                    else:
                        continue

                if not sim_ep2_list:
                    continue

                axs[i,0].loglog(w_list, theory_ep_list, label='Wrapped Gaussian (theory)', linestyle='', marker='^', markersize=8, color='black')
                axs[i,0].loglog(w_list, sim_ep_list, label='Wrapped Gaussian (simulation)', marker='x', linestyle='', color='red', markersize=8)
                axs[i,0].loglog(w_list, np.repeat(mu/s, 10), label='Uniform limit (theory)', linestyle='--', color='black')
                axs[i,0].set_title(r"$\mathbb{E}[P]$, $s=$"+str(s))
                axs[i,0].set_xlabel(r"sampling width ($w$)")
                axs[i,0].set_ylabel(r"$\mathbb{E}[P]$")
                axs[i,0].set_ylim((mu/s)*1e-1, (mu/s)*1e1)
                axs[i,0].legend(frameon=False, loc='upper left')#, title='Sampling kernel')
                axs[i,1].loglog(w_list, theory_ep2_list, label='Wrapped Gaussian (theory)', linestyle='', marker='^', markersize=8, color='black')
                # axs[1].loglog(w_list, theory_ep2_list_approx, label='Gaussian (theory approx)', linestyle='', marker='+',
                #               markersize=8, color='blue')
                axs[i,1].loglog(w_list, sim_ep2_list, label='Wrapped Gaussian (simulation)', marker='x', linestyle='', color='red', markersize=8)
                axs[i,1].loglog(w_list, np.repeat(mu / (s**2*L**2*rho), 10), label='Uniform limit (theory)', linestyle='--', color='black')
                axs[i,1].set_title(r"$\mathbb{E}[P^2]$, $s=$"+str(s))
                axs[i,1].set_xlabel(r"sampling width ($w$)")
                axs[i,1].set_ylabel(r"$\mathbb{E}[P^2]$")
                # axs[1].axvline(np.sqrt(get_lc_squared(sigma,s)), color='gray', linestyle='--')
                # fig.suptitle(f"L={L}, density={rho}, s={s}, mu={mu}, sigma={sigma}")
            plt.tight_layout()
            plt.savefig(f"plots_20240703/moments_over_w_rho{rho}_L{L}_sigma{sigma}_mu{mu}.pdf")

    # for rho in dens:
    #     for w in w_list:
    #         fig, axs = plt.subplots(1, 1, figsize=(6, 6))
    #         sim_ep2_list = []
    #         theory_ep2_list = []
    #         for L in L_list:
    #             file_path = f"results/20240521/s{s}_mu{mu}_rho{rho}_L{L}_sigma{sigma}_time1000000.0_r0.1_burnin_wrapped_norm_gaussian_w{w}*"
    #             files = glob.glob(os.path.join(file_path))
    #             if files:
    #                 sim_EP, theory_EP, sim_EP2, theory_EP2 = print_moments(files)
    #                 if sim_EP is not None:
    #                     sim_ep2_list.append(sim_EP2)
    #                     theory_ep2_list.append(theory_EP2)
    #             else:
    #                 continue
    #
    #         if not sim_ep2_list:
    #             continue
    #
    #         axs.loglog(L_list, theory_ep2_list, label='Gaussian (theory)', linestyle='', marker='^', markersize=8, color='black')
    #         axs.loglog(L_list, sim_ep2_list, label='Gaussian (simulation)', marker='x', linestyle='', color='red', markersize=8)
    #         axs.set_title(r"$\mathbb{E}[P^2]$")
    #         axs.set_xlabel(r"habitat size ($L$)")
    #         axs.set_ylabel(r"$\mathbb{E}[P^2]$")
    #         axs.legend(loc='upper right', title='Sampling kernel')
    #         fig.suptitle(f"w={round(w)}, density={rho}, s={s}, mu={mu}, sigma={sigma}")
    #         plt.tight_layout()
    #         plt.savefig(f"plots_20240604/moments_over_L_s{s}_rho{rho}_w{round(w)}.png")

if __name__ == '__main__':
    main()
