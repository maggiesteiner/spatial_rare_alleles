from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import argparse
from scipy.stats import gamma
from scipy.interpolate import interp1d

def poles(sigma,sigma_vals,pole_vals):
    f=interp1d(sigma_vals,pole_vals,fill_value="extrapolate")
    return(f(sigma))

def residues(sigma,sigma_vals, res_vals):
    res_vals=[-1*x for x in res_vals]
    f = interp1d(sigma_vals,res_vals,fill_value="extrapolate")
    return (f(sigma))

def rate_p(sigma,s,sigma_vals,pole_vals,N=10000,D=1,d=1):
    l_c=np.sqrt(D/s)
    return(s*N*(l_c**d)*poles(sigma/l_c,sigma_vals,pole_vals))

def shape_p(sigma,s,sigma_vals,res_vals,mu=1e-8,N=10000,D=1,d=1):
    l_c = np.sqrt(D / s)
    return(mu*N*(l_c**d)*residues(sigma/l_c,sigma_vals,res_vals))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--filename", type=str, help="name of input file", default="cleaned_data")
    parser.add_argument("--dim", type=int, help="dimension, default is 1", choices=[1, 2], default=1)
    parser.add_argument("--ymin",type=float,help="ymin for sfs plots",default=1e-5)
    parser.add_argument("--ymax", type=float, help="ymax for sfs plots", default=1e5)
    parser.add_argument("--pt",type=str,help="polynomial type to use",default="2_1")
    parser.add_argument("--s_list", type=float, nargs="+",help="values of selection coefficient, should have 6 values",default = [1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1])
    parser.add_argument("--sigma_list",type=float,nargs="+",help="values of sigma",default=[0.1,1,10,100])
    parser.add_argument("--N", type=float,help="population density",default=10000)
    parser.add_argument("--outname", type=str, help="name for output files (without extension)",
                        default="plots_sfs")
    parser.add_argument("--plot_both",action="store_true")
    args = parser.parse_args()
    data = pd.read_csv(args.filename + "_dim" + str(args.dim) + "_errorFalse.csv") #always use errorFalse version to get higher orders
    data = data.loc[data['poly_type'] == args.pt]

    # SFS with  fixed s, varying sigma
    sigma_vals_plot = args.sigma_list
    labs = [str(l) for l in sigma_vals_plot]
    
    

    s_vals = args.s_list
    fig, axs = plt.subplots(2, 3,figsize=(10,8))

    axs[0,0].set_title("s="+str(s_vals[0]))
    axs[0,0].set_xscale("log")
    axs[0, 0].set_yscale("log")
    axs[0, 0].set_ylim(args.ymin, args.ymax)
    colors=sns.color_palette("colorblind",len(sigma_vals_plot))
    for i in range(len(sigma_vals_plot)):
        axs[0,0].plot(x_range,gamma.pdf(x_range,
                                        a=shape_p(sigma=sigma_vals_plot[i],s=s_vals[0],sigma_vals=data['sigma'],res_vals=data['residues'],d=args.dim,N=args.N),
                                        scale=1 / rate_p(sigma=sigma_vals_plot[i], s=s_vals[0], sigma_vals=data['sigma'],
                                                      pole_vals=data['poles'],d=args.dim,N=args.N)),color=colors[i])
    axs[0, 0].legend(labels=labs, title="sigma")

    axs[0, 1].set_title("s=" + str(s_vals[1]))
    axs[0, 1].set_xscale("log")
    axs[0, 1].set_yscale("log")
    axs[0, 1].set_ylim(args.ymin, args.ymax)
    for i in range(len(sigma_vals_plot)):
        axs[0, 1].plot(x_range, gamma.pdf(x_range,
                                          a=shape_p(sigma=sigma_vals_plot[i], s=s_vals[1], sigma_vals=data['sigma'],
                                                   res_vals=data['residues'],d=args.dim),
                                          scale=1 / rate_p(sigma=sigma_vals_plot[i], s=s_vals[1],
                                                            sigma_vals=data['sigma'],
                                                            pole_vals=data['poles'],d=args.dim,N=args.N)), color=colors[i])
    axs[0, 1].legend(labels=labs, title="sigma")

    axs[0, 2].set_title("s=" + str(s_vals[2]))
    axs[0, 2].set_xscale("log")
    axs[0, 2].set_yscale("log")
    axs[0, 2].set_ylim(args.ymin, args.ymax)
    for i in range(len(sigma_vals_plot)):
        axs[0, 2].plot(x_range, gamma.pdf(x_range,
                                          a=shape_p(sigma=sigma_vals_plot[i], s=s_vals[2], sigma_vals=data['sigma'],
                                                   res_vals=data['residues'],d=args.dim,N=args.N),
                                          scale=1 / rate_p(sigma=sigma_vals_plot[i], s=s_vals[2],
                                                            sigma_vals=data['sigma'],
                                                            pole_vals=data['poles'],d=args.dim,N=args.N)), color=colors[i])
    axs[0, 2].legend(labels=labs, title="sigma")

    axs[1, 0].set_title("s=" + str(s_vals[3]))
    axs[1, 0].set_xscale("log")
    axs[1, 0].set_yscale("log")
    axs[1, 0].set_ylim(args.ymin, args.ymax)
    for i in range(len(sigma_vals_plot)):
        axs[1,0].plot(x_range, gamma.pdf(x_range,
                                          a=shape_p(sigma=sigma_vals_plot[i], s=s_vals[3], sigma_vals=data['sigma'],
                                                   res_vals=data['residues'],d=args.dim,N=args.N),
                                          scale=1 / rate_p(sigma=sigma_vals_plot[i], s=s_vals[3],
                                                            sigma_vals=data['sigma'],
                                                            pole_vals=data['poles'],d=args.dim,N=args.N)), color=colors[i])
    axs[1, 0].legend(labels=labs, title="sigma")

    axs[1, 1].set_title("s=" + str(s_vals[4]))
    axs[1, 1].set_xscale("log")
    axs[1, 1].set_yscale("log")
    axs[1, 1].set_ylim(args.ymin, args.ymax)
    for i in range(len(sigma_vals_plot)):
        axs[1, 1].plot(x_range, gamma.pdf(x_range,
                                          a=shape_p(sigma=sigma_vals_plot[i], s=s_vals[4], sigma_vals=data['sigma'],
                                                   res_vals=data['residues'],d=args.dim,N=args.N),
                                          scale=1 / rate_p(sigma=sigma_vals_plot[i], s=s_vals[4],
                                                            sigma_vals=data['sigma'],
                                                            pole_vals=data['poles'],d=args.dim,N=args.N)), color=colors[i])
    axs[1, 1].legend(labels=labs, title="sigma")

    axs[1, 2].set_title("s=" + str(s_vals[5]))
    axs[1, 2].set_xscale("log")
    axs[1, 2].set_yscale("log")
    axs[1, 2].set_ylim(args.ymin, args.ymax)
    for i in range(len(sigma_vals_plot)):
        axs[1, 2].plot(x_range, gamma.pdf(x_range,
                                          a=shape_p(sigma=sigma_vals_plot[i], s=s_vals[5], sigma_vals=data['sigma'],
                                                   res_vals=data['residues'],d=args.dim,N=args.N),
                                          scale=1 / rate_p(sigma=sigma_vals_plot[i], s=s_vals[5],
                                                            sigma_vals=data['sigma'],
                                                            pole_vals=data['poles'],d=args.dim,N=args.N)), color=colors[i])
    axs[1, 2].legend(labels=labs, title="sigma")

    plt.savefig("plots_sfs_dim" + str(args.dim) + "_polytype_" + args.pt + "_N"+str(args.N)+".png")

    # shape and rate parameters over s
    colors  = sns.color_palette("colorblind", 3)
    fig, axs = plt.subplots(1,2)
    s_range=np.linspace(10e-6,1)
    axs[0].plot(s_range,
             shape_p(s=s_range, sigma=1, sigma_vals=data['sigma'].tolist(), res_vals=data['residues'].tolist(),d=args.dim,N=args.N),color=colors[0])
    axs[0].plot(s_range,
             shape_p(s=s_range, sigma=10, sigma_vals=data['sigma'].tolist(), res_vals=data['residues'].tolist(),d=args.dim,N=args.N),color=colors[1])
    axs[0].plot(s_range,shape_p(s=s_range,sigma=100,sigma_vals=data['sigma'].tolist(),res_vals=data['residues'].tolist(),d=args.dim,N=args.N),color=colors[2])
    axs[0].legend(labels=['1','10','100'],title="sigma")
    axs[0].set_xscale("log")
    axs[0].set_yscale("log")
    axs[0].set_title("Shape parameter (mutation)")
    axs[0].set_xlabel("s")

    axs[1].plot(s_range,
                rate_p(s=s_range, sigma=1, sigma_vals=data['sigma'].tolist(), pole_vals=data['poles'].tolist(),d=args.dim,N=args.N),color=colors[0])
    axs[1].plot(s_range,
                rate_p(s=s_range, sigma=10, sigma_vals=data['sigma'].tolist(), pole_vals=data['poles'].tolist(),d=args.dim,N=args.N),color=colors[1])
    axs[1].plot(s_range,
                rate_p(s=s_range, sigma=100, sigma_vals=data['sigma'].tolist(), pole_vals=data['poles'].tolist(),d=args.dim,N=args.N),color=colors[2])
    axs[1].legend(labels=['1', '10', '100'], title="sigma")
    axs[1].set_xscale("log")
    axs[1].set_yscale("log")
    axs[1].set_title("Rate parameter (selection)")
    axs[1].set_xlabel("s")
    plt.savefig("plots_params_selection_dim"+str(args.dim)+"_polytype_"+args.pt+"_N"+str(args.N)+".png")

    if args.plot_both==True:
        data = pd.read_csv(args.filename + "_dim" + str(
            1) + "_errorFalse.csv")  # always use errorFalse version to get higher orders
        data = data.loc[data['poly_type'] == '2_1']

        data2 = pd.read_csv(args.filename + "_dim" + str(
            2) + "_errorFalse.csv")  # always use errorFalse version to get higher orders
        data2 = data2.loc[data2['poly_type'] == '1_1']

        # shape and rate parameters over s
        colors = sns.color_palette("colorblind", 3)
        fig, axs = plt.subplots(1, 2)
        s_range = np.linspace(10e-6, 1)
        axs[0].plot(s_range,
                    shape_p(s=s_range, sigma=1, sigma_vals=data['sigma'].tolist(), res_vals=data['residues'].tolist(),
                           d=1, N=10000), color=colors[0])
        axs[0].plot(s_range,
                    shape_p(s=s_range, sigma=10, sigma_vals=data['sigma'].tolist(), res_vals=data['residues'].tolist(),
                           d=1, N=10000), color=colors[1])
        axs[0].plot(s_range,
                    shape_p(s=s_range, sigma=100, sigma_vals=data['sigma'].tolist(), res_vals=data['residues'].tolist(),
                           d=1, N=10000), color=colors[2])

        axs[0].plot(s_range,
                    shape_p(s=s_range, sigma=1, sigma_vals=data2['sigma'].tolist(), res_vals=data2['residues'].tolist(),
                           d=2, N=100), color=colors[0],linestyle="--")
        axs[0].plot(s_range,
                    shape_p(s=s_range, sigma=10, sigma_vals=data2['sigma'].tolist(), res_vals=data2['residues'].tolist(),
                           d=2, N=100), color=colors[1],linestyle="--")
        axs[0].plot(s_range,
                    shape_p(s=s_range, sigma=100, sigma_vals=data2['sigma'].tolist(), res_vals=data2['residues'].tolist(),
                           d=2, N=100), color=colors[2],linestyle="--")

        axs[0].legend(labels=['1', '10', '100'], title="sigma")
        axs[0].set_xscale("log")
        axs[0].set_yscale("log")
        axs[0].set_title("Shape parameter (mutation)")
        axs[0].set_xlabel("s")

        axs[1].plot(s_range,
                    rate_p(s=s_range, sigma=1, sigma_vals=data['sigma'].tolist(), pole_vals=data['poles'].tolist(),
                            d=1, N=10000), color=colors[0])
        axs[1].plot(s_range,
                    rate_p(s=s_range, sigma=10, sigma_vals=data['sigma'].tolist(), pole_vals=data['poles'].tolist(),
                            d=1, N=10000), color=colors[1])
        axs[1].plot(s_range,
                    rate_p(s=s_range, sigma=100, sigma_vals=data['sigma'].tolist(), pole_vals=data['poles'].tolist(),
                            d=1, N=10000), color=colors[2])

        axs[1].plot(s_range,
                    rate_p(s=s_range, sigma=1, sigma_vals=data2['sigma'].tolist(), pole_vals=data2['poles'].tolist(),
                            d=2, N=100), color=colors[0],linestyle="--")
        axs[1].plot(s_range,
                    rate_p(s=s_range, sigma=10, sigma_vals=data2['sigma'].tolist(), pole_vals=data2['poles'].tolist(),
                            d=2, N=100), color=colors[1],linestyle="--")
        axs[1].plot(s_range,
                    rate_p(s=s_range, sigma=100, sigma_vals=data2['sigma'].tolist(), pole_vals=data2['poles'].tolist(),
                            d=2, N=100), color=colors[2],linestyle="--")

        axs[1].legend(labels=['1', '10', '100'], title="sigma")
        axs[1].set_xscale("log")
        axs[1].set_yscale("log")
        axs[1].set_title("Rate parameter (selection)")
        axs[1].set_xlabel("s")
        plt.savefig("plots_params_selection_both_N10000_N100.png")


if __name__ == '__main__':
    main()
