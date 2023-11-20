import numpy as np
import pandas as pd
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--filename", type=str, help="name of input file", default="pade_approx")
    parser.add_argument("--dim", type=int, help="dimension, default is 1", choices=[1, 2], default=1)
    parser.add_argument("--calc_error",action='store_true')
    parser.add_argument("--outname", type=str, help="name for output files (without extension)",
                        default="cleaned_data")
    args = parser.parse_args()

    data = pd.read_csv(args.filename+"_dim"+str(args.dim)+"_error"+str(args.calc_error)+".csv")

    sigma_list = [float(x) for x in data['sigma'].tolist()]
    residues = [np.real(complex(x)) for x in data['residue'].tolist()]
    poles = [np.real(complex(x)) for x in data['pole'].tolist()]
    rem=[np.real(complex(x)) for x in data['remainder'].tolist()]
    data["poly_type"] = data["m"].astype(str) + "_" + data["n"].astype(str)
    data["total_deg"] = data["m"].astype(int)+data["n"].astype(int)
    poly_type = data["poly_type"].tolist()
    pole_2 = data["pole_2"].tolist()
    res_2 = data["res_2"].tolist()
    total_deg = data["total_deg"].tolist()
    if args.calc_error==True:
        error_next = data["error_next"].tolist()
        rel_error = data["rel_err"].tolist()
        df = pd.DataFrame(
            dict(sigma=sigma_list, residues=residues, poles=poles, poly_type=poly_type, pole_2=pole_2, res_2=res_2,
                 remainder=rem,error_next=error_next,rel_error=rel_error,total_deg=total_deg))
    else:
        df = pd.DataFrame(dict(sigma=sigma_list, residues=residues, poles=poles,poly_type=poly_type,pole_2=pole_2,res_2=res_2,remainder=rem))#,error_next=error_next,rel_error=rel_error

    df.to_csv(args.outname+"_dim"+str(args.dim)+"_error"+str(args.calc_error)+".csv",index=False)

if __name__ == '__main__':
    main()


