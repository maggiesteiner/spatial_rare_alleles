import numpy as np
from scipy.special import expi
import argparse
import os
import json

def get_lc_squared(sigma, s):
    return sigma ** 2 / s

def u2_exact(wt):
    u2 = (-1 * np.exp(wt ** 2 / (4 * np.pi ** 2)) / (8 * np.pi)) * expi(-1 * wt ** 2 / (4 * np.pi ** 2))
    return u2

def get_EPsquared_theory(mu, s, rho, sigma, w):
    lcs = get_lc_squared(sigma, s)
    u2 = u2_exact(w / np.sqrt(lcs))
    return ((2 * mu / (s ** 2 * rho * lcs)) * u2 + mu ** 2 / (s ** 2))

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


def print_moments(f):
    with open(f,"r") as json_file:
        data = json.load(json_file)
    s = data["s"]
    mu = data["mu"]
    rho = data["rho"]
    sigma = data["sigma"]
    w = data["w"]
    sampled_p_flattened = np.array(data['sampled_p_flattened'])
    zero_samples = data['zero_samples']

    print(f"File path: {f}")
    print(f"s: {s}")
    print(f"mu: {mu}")
    print(f"rho: {rho}")
    print(f"sigma: {sigma}")
    print(f"w: {w}\n")

    # first moment
    sim_EP = get_EP_sim(sampled_p_flattened, zero_samples)
    theory_EP = get_EP_theory(mu, s)

    print(f"EP (sim): {sim_EP}")
    print(f"EP (theory): {theory_EP}\n")

    # second moment
    sim_EP2 = get_EPsquared_sim(sampled_p_flattened, zero_samples)
    theory_EP2 = get_EPsquared_theory(mu=mu, s=s, rho=rho, sigma=sigma, w=w)
    theory_EP2_2 = get_EPsquared_theory(mu=mu, s=s, rho=rho, sigma=w, w=sigma)
    theory_EP2_3 = get_EPsquared_theory_swap(mu=mu, s=s, rho=rho, sigma=sigma, w=w)

    print(f"EP2 (sim): {sim_EP2}")
    print(f"EP2 (theory): {theory_EP2}")
    print(f"EP2 (theory swapped): {theory_EP2_2}")
    print(f"EP2 (theory swapped v2): {theory_EP2_3}")
    print("-----------------------\n")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", type=str, help="file path (json)")
    args = parser.parse_args()
    print_moments(args.f)

if __name__ == '__main__':
    main()
