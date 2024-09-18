import sys
# sys.path.insert(0, '/home/mianniravn/.local/lib/python3.7/site-packages')
#from matplotlib import pyplot as plt
import pyslim as ps
import numpy as np
from scipy.stats import multivariate_normal as mvn
import tskit
import pandas as pd
import argparse
import os

# Parse command-line arguments
parser = argparse.ArgumentParser(description='Calculate and save SFS data.')
parser.add_argument('input_file', help='File to be processed')
parser.add_argument('--num_samples', type=int, default=100, help='Number of individuals to sample')
parser.add_argument('--output', type=str, help='Output path')
parser.add_argument('--niter', type=int, default=1000, help='Number of iterations for sampling')
args = parser.parse_args()

input_file = args.input_file
num_to_sample = args.num_samples
output = args.output
niter = args.niter

# Create the output filename with the "sfs" directory included
num_samples_str = str(num_to_sample)

# Load tree sequence
slim_ts = tskit.load(input_file)
slim_ts = ps.update(slim_ts)
locs = slim_ts.individual_locations  # locations of  all individuals

samps_at_end = ps.individuals_alive_at(slim_ts, 0)  # individuals alive at end

# genotype matrix - actually a haplotype matrix
H = slim_ts.genotype_matrix(isolated_as_missing=False)

# Assuming G is a numpy array of shape (L, 2N)
L, N = H.shape[0], H.shape[1] // 2

# Reshape G to (L, N, 2) to separate the two genomes
G_reshape = np.reshape(H, (L, N, 2))

# Sum the two genomes in each individual to get a genotype matrix
G = np.sum(G_reshape, axis=2)
G[G < 1] = 0  # fixes an occasional mysterious -1 which appears


def get_weights(lcs, mean=25, var=5):
    # gets sampling weights of individuals
    wts = np.zeros(len(lcs[:, 0]))
    sampkernel = mvn(mean=[mean, mean], cov=[[var, 0], [0, var]])

    for i in range(len(lcs[:, 0])):
        wts[i] = sampkernel.pdf(lcs[i, 0:2])

    wts /= np.sum(wts)  # normalize weights

    return wts


# Gaussian samples with different sds
sd_range_disp = np.arange(5, 205, 20)
sd_range = sd_range_disp * np.sqrt(0.06)

variant_data={}
for sd_value in sd_range:
    print(f"Sampling sd {sd_value}.")
    variant_freqs = [0] * G.shape[0]
    for it in range(niter):
        var_value = sd_value ** 2 # get variance
        weights = get_weights(locs, 25, var_value) # get weights
        samples = np.random.choice(samps_at_end, size=num_to_sample, replace=True, p=weights) # sample individuals
        samples = list(samples.astype(int))
        # get a vector of variant frequencies
        # and average over sampling iterations
        variant_freqs += (np.sum(G[:,samples],axis=1)/np.shape(G)[1])/niter
    variant_data[sd_value] = variant_freqs.tolist()

print(f"St. dev of x coordinates of all individuals: {np.std(locs[:, 0])}")

# Create DataFrame from the SFS data
variant_df = pd.DataFrame(variant_data)

# Save SFS DataFrame to file
print("Saving file to")
print(output)
variant_df.to_csv(output, sep='\t', index=False, na_rep="NA")




