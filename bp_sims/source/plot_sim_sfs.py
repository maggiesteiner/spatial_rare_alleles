import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import argparse
from scipy.special import loggamma # type: ignore
import re

def finite_sfs_k_unif(k,n,s,mu,N):
    gammae = s*N
    thetae = mu*N
    logval = k*np.log(n)+thetae*np.log(gammae)-(k+thetae)*np.log(n+gammae)+loggamma(k+thetae)-loggamma(k+1)-loggamma(thetae)
    return(np.e**logval)

def plot_sim_sfs(input_file,output_file,plot_theory,n,s,mu,N,max_x=501):
    df = pd.read_csv(input_file)
    fig, axs = plt.subplots(1,1,figsize=(8,8))
    axs.loglog(np.arange(1,max_x),df[1:max_x],label='simulation',marker='o')
    axs.set_ylim(1e-8,1e-1)
    if plot_theory is True:
        sfs_theory = [finite_sfs_k_unif(k,n,s,mu,N) for k in range(1,max_x,1)]
        axs.loglog(np.arange(1,max_x),sfs_theory,label='theory')
    axs.set_xlabel('allele count')
    axs.set_ylabel('expected value')
    axs.set_title("n="+str(n)+", s="+str(s)+", mu="+str(mu)+", N="+str(N))
    plt.legend()
    plt.savefig(output_file)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--sfs_file', type=str, help='*.sfs file from simulation')
    parser.add_argument('--filename', type=str, help='output file for plot')
    parser.add_argument('--plot_theory',action='store_true',help='plot theory SFS on top of sim')
    args = parser.parse_args()

    pattern = r"s([\d.e-]+)_n(\d+)_mu([\de.-]+)_rho(\d+)_L(\d+)"
    match = re.search(pattern, args.sfs_file)
    s = float(match.group(1))
    n = int(match.group(2))
    mu = float(match.group(3))
    rho = int(match.group(4))
    L = int(match.group(5))

    plot_sim_sfs(args.sfs_file,args.filename,args.plot_theory,n,s,mu,rho*L**2)

if __name__ == '__main__':
    main()
