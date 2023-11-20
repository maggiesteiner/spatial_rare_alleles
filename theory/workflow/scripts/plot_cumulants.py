from matplotlib import pyplot as plt
import pandas as pd
import argparse
import numpy as np

def plot_cumulants(xvals,lists,labs,xlab="spatial dispersion of sample (sigma)",ylab="value",outname="plt.png",xtrans=None,ytrans=None):
    plt.figure()
    for i in range(len(lists)):
        plt.plot(xvals,lists[i],label=labs[i])
    plt.legend()
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    if xtrans is not None:
        plt.xscale(xtrans)
    if ytrans is not None:
        plt.yscale(ytrans)
    plt.title("Value of cumulants")
    plt.savefig(outname)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--filename", type=str, help="name of input file", default="spatial_integrals")
    parser.add_argument("--dim", type=int, help="dimension, default is 1", choices=[1, 2], default=1)
    parser.add_argument("--outname", type=str, help="name for output files (without extension)",
                        default="cumulants")
    args = parser.parse_args()

    data = pd.read_csv(args.filename+"_dim"+str(args.dim)+".csv")
    sigma_list = data["sigma"]

    if args.dim==1:
        plot_cumulants(sigma_list,
                       lists=[np.repeat(1,len(data["u2_GQ"])),data["u2_GQ"],data["u3_GQ"],data["u4_GQ"]],
                       labs=["u1","u2","u3","u4"],outname=args.outname+"_dim"+str(args.dim))

    elif args.dim==2:
        plot_cumulants(sigma_list,
                       lists=[np.repeat(1, len(data["u2_GQ"])), data["u2_GQ"], data["u3_GQ"]],
                       labs=["u1", "u2", "u3"], outname=args.outname + "_dim" + str(args.dim))

if __name__ == '__main__':
    main()