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

np.random.seed(1234)

# Parse command-line arguments
parser = argparse.ArgumentParser(description='Calculate and save SFS data.')
parser.add_argument('input_file', help='File to be processed')
parser.add_argument('--num_samples', type=int, default=100, help='Number of individuals to sample')
parser.add_argument('--output', type=str, help='Output path')
parser.add_argument('--niter', type=int, default=1000, help='Number of iterations for sampling')
parser.add_argument('--bins', type=int, default=50, help='Number of bins')
args = parser.parse_args()

input_file = args.input_file
num_to_sample = args.num_samples
output = args.output
niter = args.niter
n_bins = args.bins

# Create output directory for SFS files
output_dir = os.path.join(input_file, 'sfs', str(num_to_sample))


# Construct the file paths
# Get the directory and base filename from the input file path
file_dir = os.path.dirname(input_file)
base_filename = os.path.splitext(os.path.basename(input_file))[0]

# Create the output filename with the "sfs" directory included
num_samples_str = str(num_to_sample)
output_filename = os.path.join(file_dir, "sfs_log_niter", base_filename + '_' + num_samples_str + '.sfs.tsv')

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
sd_range_disp = [0,5,10,50,100,185]
sd_range_disp = np.array(sd_range_disp,dtype=float)
sd_range = sd_range_disp * np.sqrt(0.06)

sfs_data = {}

# Create an array 's' containing 25 equally spaced values in log space
s = np.linspace(np.log10(1), np.log10(num_to_sample*2), num=n_bins)

# Compute 'bin_edges' and add a zero at the beginning
bin_edges = np.concatenate(([0.1], np.power(10, s)))

# Calculate 'integer_bin_edges' and 'bin_widths'
integer_bin_edges = np.floor(bin_edges).astype(int)
bin_widths = np.where(np.diff(integer_bin_edges) == 0, 1, np.diff(integer_bin_edges))

# Add a zero at the beginning of array 's'
s = np.concatenate(([0], s))

# Calculate 'bin_centers'
bin_centers = np.power(10, (s[:-1] + s[1:]) / 2)

for sd_value in sd_range:
    var_value = sd_value ** 2
    if var_value == 0:
        weights = np.ones(len(locs[:, 0])) / len(locs[:, 0])  # Uniform weights for random sample
    else:
        weights = get_weights(locs, 25, var_value)

    for i in range(niter):
        samples = np.random.choice(samps_at_end, size=num_to_sample, replace=True, p=weights)
        samples = list(samples.astype(int))
        counts = np.sum(G[:, samples], 1)  # get counts by summing over individuals
        mut_afs, edg = np.histogram(counts, bins=bin_edges)
        mut_afs = mut_afs / bin_widths  # normalise by bin width
        mut_afs[np.isinf(mut_afs)] = np.nan  # changed
        if i==0:
            mut_arr = mut_afs
        else:
            mut_arr = np.vstack((mut_arr,mut_afs))
    sfs_data[sd_value] = np.mean(mut_arr,axis=0)
    print(f"SFS for width {sd_value} computed, peek {mut_afs[:10]}.")

# Calculate SFS for each sample

print(f"St. dev of x coordinates of all individuals: {np.std(locs[:, 0])}")

# Create DataFrame from the SFS data
sfs_data['allele counts'] = np.array(bin_centers)

sfs_df = pd.DataFrame(sfs_data)

# Save SFS DataFrame to file
print("Saving file to")
print(output)
sfs_df.to_csv(output, sep='\t', index=False, na_rep="NA")




