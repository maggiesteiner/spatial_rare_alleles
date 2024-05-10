import numpy as  np
from scipy.special import expi
import argparse
import os

def get_lc_squared(sigma, s):
    return sigma ** 2 / s

def u2_exact(wt):
    u2 = (-1 * np.exp(wt ** 2 / (4 * np.pi ** 2)) / (8 * np.pi)) * expi(-1 * wt ** 2 / (4 * np.pi ** 2))
    return u2

def get_EPsquared_theory(mu, s, rho, sigma, w):
    lcs = get_lc_squared(sigma, s)
    u2 = u2_exact(w / np.sqrt(lcs))
    wt = w / np.sqrt(lcs)
    return ((2 * mu / (s ** 2 * rho * lcs)) * u2 + mu ** 2 / (s ** 2))

def get_EP_theory(mu, s):
    return mu / s

def get_EPsquared_sim(ps, zeros):
    return (np.sum(ps ** 2) / (len(ps) + zeros))


def get_EP_sim(ps, zeros):
    return np.sum(ps) / (len(ps) + zeros)


def print_moments(f):
    # read params
    prefix = os.path.basename(f)
    parameters = prefix.split("_")
    s = float(parameters[0][1:])
    mu = float(parameters[1][2:])
    rho = int(parameters[2][3:])
    sigma = int(parameters[4][5:])
    w = int(parameters[11][1:])

    print(f"File prefix: {prefix}")
    print(f"s: {s}")
    print(f"mu: {mu}")
    print(f"rho: {rho}")
    print(f"sigma: {sigma}")
    print(f"w: {w}\n")

    # read files
    pfile = f + ".p"
    zfile = f + ".zero"
    ps = np.loadtxt(pfile)
    zeros = np.loadtxt(zfile)

    # first moment
    sim_EP = get_EP_sim(ps,zeros)
    theory_EP = get_EP_theory(mu,s)

    print(f"EP (sim): {sim_EP}")
    print(f"EP (theory): {theory_EP}\n")

    # second moment
    sim_EP2 = get_EPsquared_sim(ps,zeros)
    theory_EP2 = get_EPsquared_theory(mu, s, rho, sigma, w)

    print(f"EP2 (sim): {sim_EP2}")
    print(f"EP2 (theory): {theory_EP2}\n")
    print("-----------------------\n")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", type=str, help="file prefix")
    args = parser.parse_args()
    print_moments(args.f)

if __name__ == '__main__':
    main()
