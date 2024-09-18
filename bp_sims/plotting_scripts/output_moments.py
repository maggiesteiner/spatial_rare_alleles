import numpy as np
from scipy.special import expi, exp1
import argparse
import os
import json
import re
import glob

def get_lc_squared(sigma, s):
    return sigma ** 2 / s

def u2_exact(wt):
    u2 = (-1 * np.exp(wt ** 2 / (4 * np.pi ** 2)) / (8 * np.pi)) * expi(-1 * wt ** 2 / (4 * np.pi ** 2))
    return u2

def get_EPsquared_theory(mu, s, rho, sigma, w):
    lcs = get_lc_squared(sigma, s)
    u2 = u2_exact(w / np.sqrt(lcs))
    return ((2 * mu / (s ** 2 * rho * lcs)) * u2 + mu ** 2 / (s ** 2))

def get_EPsquared_theory_finitehabitat(mu, s, rho, sigma, w):
    lcs = get_lc_squared(sigma, s)
    return (mu/(s**2*rho*4*np.pi*lcs))*np.exp((w/np.sqrt(lcs))**2)*exp1((w/np.sqrt(lcs))**2)+mu**2/s**2

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
    theory_EP2_swap = get_EPsquared_theory(mu=mu, s=s, rho=rho, sigma=w, w=sigma)
    theory_EP2_finitehabitat = get_EPsquared_theory_finitehabitat(mu=mu, s=s, rho=rho, sigma=sigma, w=w)
    theory_EP2_finitehabitat_swap = get_EPsquared_theory_finitehabitat(mu=mu, s=s, rho=rho, sigma=w, w=sigma)

    print(f"EP2 (sim): {sim_EP2}")
    # print(f"EP2 (theory infinite habitat): {theory_EP2}")
    # print(f"EP2 (theory infinite habitat) (swap): {theory_EP2_swap}")
    print(f"EP2 (theory finite habitat): {theory_EP2_finitehabitat}")
    # print(f"EP2 (theory finite habitat) (swap): {theory_EP2_finitehabitat_swap}")
    print("-----------------------\n")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", type=str, help="file path (json)")
    args = parser.parse_args()

    base_name = re.sub(r'_rep\d+\.json$', '', os.path.basename(args.f))
    directory = os.path.dirname(args.f)

    files = glob.glob(os.path.join(directory, base_name + "*"))
    print(files)
    print_moments(files)

if __name__ == '__main__':
    main()
