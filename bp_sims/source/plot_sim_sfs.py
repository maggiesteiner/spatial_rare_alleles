import argparse
import re

import numpy as np
import pandas as pd  # type: ignore
from matplotlib import pyplot as plt  # type: ignore
from scipy.special import loggamma  # type: ignore
from scipy.interpolate import interp1d

def finite_sfs_k_unif(k,n,s,mu,N):
    gammae = s*N
    thetae = mu*N
    logval = k*np.log(n)+thetae*np.log(gammae)-(k+thetae)*np.log(n+gammae)+loggamma(k+thetae)-loggamma(k+1)-loggamma(thetae)
    return(np.e**logval)

def plot_sim_sfs_unif(input_file,output_file,plot_theory,n,s,mu,N,max_x=501):
    df = pd.read_csv(input_file,header=None)
    #print(df)
    fig, axs = plt.subplots(1,1,figsize=(8,8))
    axs.loglog(np.arange(0,max_x),df[0:max_x],label='simulation',marker='o')
    axs.set_ylim(1e-8,1e-1)
    if plot_theory is True:
        sfs_theory = [finite_sfs_k_unif(k,n,s,mu,N) for k in range(1,max_x,1)]
        axs.loglog(np.arange(1,max_x),sfs_theory,label='theory')
    axs.set_xlabel('allele count')
    axs.set_ylabel('expected value')
    axs.set_title("n="+str(n)+", s="+str(s)+", mu="+str(mu)+", N="+str(N))
    plt.legend()
    plt.savefig(output_file)

### functions for theory ###
def poles(w,w_vals,pole_vals):
    f=interp1d(w_vals,pole_vals,fill_value="extrapolate")
    return(f(w))

def residues(w,w_vals, res_vals):
    res_vals=[-1*x for x in res_vals]
    f = interp1d(w_vals,res_vals,fill_value="extrapolate")
    return (f(w))

def get_gammae(w,s,w_vals,pole_vals,N=10000,sigma=1,d=1):
    l_c=np.sqrt(sigma/s)
    return(s*N*(l_c**d)*poles(w/l_c,w_vals,pole_vals))

def get_thetae(w,s,w_vals,res_vals,mu=1e-8,N=10000,sigma=1,d=1):
    l_c = np.sqrt(sigma / s)
    return(mu*N*(l_c**d)*residues(w/l_c,w_vals,res_vals))

def finite_sfs_k_gaussian(n,k,w,s,w_vals,pole_vals,res_vals,mu=1e-7,rho=100,sigma=1,d=2):
    gammae = get_gammae(w,s,w_vals,pole_vals,rho,sigma,d)
    thetae = get_thetae(w,s,w_vals,res_vals,mu,rho,sigma,d)
    logval = k*np.log(n)+thetae*np.log(gammae)-(k+thetae)*np.log(n+gammae)+loggamma(k+thetae)-loggamma(k+1)-loggamma(thetae)
    return(np.e**logval)

def plot_sim_sfs_gaussian(input_file,output_file,plot_theory,n,w,s,w_vals,pole_vals,res_vals,mu=1e-7,rho=100,sigma=1,d=2,max_x=501):
    df = pd.read_csv(input_file, header=None)
    # print(df)
    fig, axs = plt.subplots(1, 1, figsize=(8, 8))
    axs.loglog(np.arange(0, max_x), df[0:max_x], label='simulation', marker='o')
    axs.set_ylim(1e-8, 1e-1)
    if plot_theory is True:
        sfs_theory = [finite_sfs_k_gaussian(n,k,w,s,w_vals,pole_vals,res_vals,mu,rho,sigma,d=2) for k in range(1, max_x, 1)]
        axs.loglog(np.arange(1, max_x), sfs_theory, label='theory')
    axs.set_xlabel('allele count')
    axs.set_ylabel('expected value')
    axs.set_title("n=" + str(n) + ", s=" + str(s) + ", mu=" + str(mu) + ", w=" + str(w))
    plt.legend()
    plt.savefig(output_file)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--sfs_file', type=str, help='*.sfs file from simulation')
    parser.add_argument('--filename', type=str, help='output file for plot')
    parser.add_argument('--plot_theory',action='store_true',help='plot theory SFS on top of sim')
    parser.add_argument('--gaussian',action='store_true',help='gaussian sampling kernel')
    parser.add_argument('-w',type=float,help='sampling width',default=None)
    parser.add_argument('--sigma', type=float, help='diffusion param', default=None)
    args = parser.parse_args()

    pattern = r"s([\d.e-]+)_n(\d+)_mu([\de.-]+)_rho(\d+)_L(\d+)"
    match = re.search(pattern, args.sfs_file)
    s = float(match.group(1))
    n = int(match.group(2))
    mu = float(match.group(3))
    rho = int(match.group(4))
    L = int(match.group(5))

    if args.gaussian is not True:
        plot_sim_sfs_unif(args.sfs_file,args.filename,args.plot_theory,n,s,mu,rho*L**2)

    elif args.gaussian is True:

        data = pd.read_csv("../../theory/old_files/results/spatial_integrals_dim2.csv")
        data_pr = pd.read_csv("../../theory/old_files/results/cleaned_data_dim2_errorFalse.csv")
        data_pr = data_pr.loc[data_pr['poly_type'] == '1_1']
        w_vals = data['w'].tolist()
        res_vals = data_pr['residues']
        pole_vals = data_pr['poles']

        plot_sim_sfs_gaussian(args.sfs_file,args.filename,args.plot_theory,n,args.w,s,w_vals,pole_vals,res_vals,mu,rho,args.sigma)

if __name__ == '__main__':
    main()
