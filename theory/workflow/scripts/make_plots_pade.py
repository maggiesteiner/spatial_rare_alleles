from matplotlib import pyplot as plt
import pandas as pd
import seaborn as sns
import argparse

def plot_scatter(df,y,x="sigma",outname="fig.png",xtrans=None,ytrans=None,ylab=None,linthreshy=None):
    sns.lmplot(x, y, data=df, hue='poly_type', fit_reg=False,scatter_kws={"s": 2},legend=False,palette="colorblind")
    if xtrans is not None:
        plt.xscale(xtrans)
    if ytrans is not None and linthreshy is not None:
        plt.yscale(ytrans,linthreshy=linthreshy)
    elif ytrans is not None:
        plt.yscale(ytrans)
    if ylab is not None:
        plt.ylabel(ylab)
    plt.legend(markerscale=5)
    plt.savefig(outname,dpi=300)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--filename", type=str, help="name of input file", default="cleaned_data")
    parser.add_argument("--dim", type=int, help="dimension, default is 1", choices=[1, 2], default=1)
    parser.add_argument("--calc_error", action='store_true')
    parser.add_argument("--outname", type=str, help="name for output files (without extension)",
                        default="plots_pade")
    args = parser.parse_args()

    data = pd.read_csv(args.filename + "_dim" + str(args.dim) + "_error" + str(args.calc_error) + ".csv")

    keep_list=['1_0','1_1','1_2','2_0','2_1','3_0']
    data = data[data['poly_type'].isin(keep_list)]

    plot_scatter(data, y='poles',
                 outname=args.outname + "_poles_dim" + str(args.dim) + "_error" + str(args.calc_error) + ".png",
                 xtrans='log', ytrans='log')

    plot_scatter(data,y='residues',
                 outname=args.outname+"_residues_dim"+str(args.dim)+"_error"+str(args.calc_error)+".png",
                 xtrans='log',ytrans='symlog',linthreshy=1e-1)

    plot_scatter(data,y='remainder',
                 outname=args.outname+"_remainder_dim"+str(args.dim)+"_error"+str(args.calc_error)+".png",
                 xtrans='log')

    data['residues']=-1*data['residues']
    plot_scatter(data, x='poles', y='residues',
                 outname=args.outname + "_poles_vs_residues_dim" + str(args.dim) + "_error" + str(
                     args.calc_error) + ".png",
                 xtrans='log',ytrans='log',ylab='-1*residues')

    if args.calc_error==True:
        plot_scatter(data, y='error_next',
                     outname=args.outname + "_error_dim" + str(args.dim) + "_error" + str(
                         args.calc_error) + ".png")

        plot_scatter(data, y='rel_error',
                     outname=args.outname + "_rel_error_dim" + str(args.dim) + "_error" + str(
                         args.calc_error) + ".png")

        data['error_pole'] = abs(data['error_next']*data['poles']**(data['total_deg']+2))
        plot_scatter(data, y='error_pole',
                     outname=args.outname + "_pole_error_dim" + str(args.dim) + "_error" + str(
                         args.calc_error) + ".png")





if __name__ == '__main__':
    main()


