import numpy as np
from scipy.special import expi, exp1
import argparse
import os
import json
import re
import glob
import matplotlib.pyplot as plt

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
    return (mu / (s ** 2 * rho * 4 * np.pi * lcs)) * np.exp((w / np.sqrt(lcs)) ** 2) * exp1(
        (w / np.sqrt(lcs)) ** 2) + mu ** 2 / s ** 2

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

def extract_w_from_filename(filename):
    match = re.search(r'_w(\d+)_', filename)
    if match:
        return int(match.group(1))
    else:
        return None

def calculate_moments(files):
    combined_sampled_p, combined_zero_samples = concatenate_data(files)
    sim_EP2 = get_EPsquared_sim(combined_sampled_p, combined_zero_samples)
    return sim_EP2

def plot_data(ws, sim_EP2s, uniform_EP2):
    plt.figure(figsize=(10, 6))
    plt.plot(ws, sim_EP2s, 'o-', label='Simulated EP^2')
    plt.axhline(y=uniform_EP2, color='r', linestyle='--', label='Uniform EP^2')
    plt.xlabel('w')
    plt.ylabel('EP^2 (simulated)')
    plt.title('Simulated EP^2 vs w')
    plt.legend()
    plt.grid(False)
    plt.show()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", type=str, help="directory path containing JSON files")
    args = parser.parse_args()

    files = glob.glob(os.path.join(args.d, "*.json"))

    ws = []
    sim_EP2s = []
    uniform_files = []

    w_files_map = {}

    for file in files:
        if "uniform" in file:
            uniform_files.append(file)
        else:
            w = extract_w_from_filename(file)
            if w is not None:
                if w not in w_files_map:
                    w_files_map[w] = []
                w_files_map[w].append(file)

    for w, files in w_files_map.items():
        sim_EP2 = calculate_moments(files)
        ws.append(w)
        sim_EP2s.append(sim_EP2)

    if uniform_files:
        uniform_EP2 = calculate_moments(uniform_files)
    else:
        uniform_EP2 = None

    ws, sim_EP2s = zip(*sorted(zip(ws, sim_EP2s)))
    plot_data(ws, sim_EP2s, uniform_EP2)

if __name__ == '__main__':
    main()
