import numpy as np
import numpy.linalg
from scipy.interpolate import pade
from scipy.signal import residue
import pandas as pd
import math
import argparse

def deriv_rational(p,q,n=1):
    """
    Generate nth derivative of rational p/q using quotient rule
    """
    p = np.poly1d(p)
    q = np.poly1d(q)
    for i in range(n):
        p_der = np.polyder(p)
        q_der = np.polyder(q)
        new_p = np.polysub(np.polymul(q,p_der),np.polymul(p,q_der))
        new_q = np.polymul(q,q)
        p = new_p
        q = new_q
    return p,q


def get_taylor_coefs(p,q,m,a=0):
    """
    Compute Taylor series coefficients for rational function p/q
    """
    coefs = []
    for i in range(m):
        p_i, q_i = deriv_rational(p=p,q=q,n=i)
        coefs.append((p_i(a)/q_i(a))/math.factorial(i))
    return coefs


def get_residues(p,q):
    """
    Given function p/q, returns pole closest to zero and its residue
    """
    rs, pl, k = residue(p,q)
    mult = []
    pl_all = []
    rs_all = []
    if len(pl)>0:
        pl_abs = [abs(x) for x in pl]
        pole = pl[np.argmin(pl_abs)]
        res = rs[np.argmin(pl_abs)]
        mult_temp = len([x for x in pl if x==pole])
        mult.append(mult_temp)
        pl_all.append(pl)
        rs_all.append(rs)
    else:
        pole = np.nan
        res = np.nan
    return pole, res, mult, pl_all, rs_all, k


def calc_pade_table(coefs,calc_error=False):
    """
    Creates table containing coefficients for Pade approximation
    given coefficients (coefs) over all possible pairs (m,n) and
    associated Taylor series coefficients
    """
    mvals = []
    nvals = []
    pvals = []
    qvals = []
    taylor_coefs = []
    if calc_error==True:
        for m in range(0,len(coefs)-1):#len(coefs)-1
            for n in range(0,len(coefs)-m-1):#len(coefs)-1-m)
                try:
                    p, q = pade(coefs,m,n)
                    pvals.append(p.coefficients)
                    qvals.append(q.coefficients)
                    mvals.append(m)
                    nvals.append(n)
                    taylor = get_taylor_coefs(p,q,len(coefs))
                    taylor_coefs.append(taylor)
                except numpy.linalg.LinAlgError as err:
                    print("Warning: error thrown for m = %s and n = %s.\n %s" % (m,n,err))
    else:
        for m in range(0,len(coefs)):#len(coefs)-1
            for n in range(0,len(coefs)-m):#len(coefs)-1-m)
                try:
                    p, q = pade(coefs,m,n)
                    pvals.append(p.coefficients)
                    qvals.append(q.coefficients)
                    mvals.append(m)
                    nvals.append(n)
                    taylor = get_taylor_coefs(p,q,len(coefs))
                    taylor_coefs.append(taylor)
                except numpy.linalg.LinAlgError as err:
                    print("Warning: error thrown for m = %s and n = %s.\n %s" % (m,n,err))

    pade_tab = pd.DataFrame(list(zip(mvals,nvals,pvals,qvals,taylor_coefs)),columns=['m','n','p','q','tc'])
    return pade_tab

def calc_pole_res(tab):
    """
    Calculates pole nearest to zero and its residue for each Pade
    approximant given in tab and appends to data frame
    """
    polevals = []
    resvals = []
    multvals = []
    pl2vals = []
    res2vals = []
    remvals=[]
    kvals=[]
    for row in tab.itertuples():
        if row[tab.columns.get_loc('m')+1]>0:
            p, q = get_pade_poly(tab,row[tab.columns.get_loc('m')+1],row[tab.columns.get_loc('n')+1])
            pole, res, mult, pl_all, rs_all, k = get_residues(p, q)
            polevals.append(pole)
            resvals.append(res)
            multvals.append(mult)
            remainder=0
            if len(pl_all[0])>1: # if more than 1 pole, add remainder terms (other than k)
                remainder=0
                for i in range(len(pl_all[0])): # for  each  pole
                    p_temp=pl_all[0][i]
                    if p_temp != pole: # if current pole isn't pole of interest
                        remainder+=np.real(rs_all[0][i])/(pole-p_temp) # add to remainder
            temp = np.poly1d(k)
            remainder += temp(pole)  # add  k  term
            remvals.append(remainder)
            kvals.append(k)
            if len(pl_all[0])>1:
                pl2 = [x for x in pl_all[0] if x != pole][0]
                res2 = [x for x in rs_all[0] if x != res][0]
                pl2vals.append(pl2)
                res2vals.append(res2)
            else:
                pl2vals.append(np.nan)
                res2vals.append(np.nan)
        else:
            polevals.append(np.nan)
            resvals.append(np.nan)
            multvals.append(np.nan)
            pl2vals.append(np.nan)
            res2vals.append(np.nan)
            remvals.append(np.nan)
            kvals.append(np.nan)
    tab_new = tab
    tab_new['pole'] = polevals
    tab_new['residue'] = resvals
    tab_new['multiplicity_pole'] = multvals
    tab_new['pole_2'] = pl2vals
    tab_new['res_2'] = res2vals
    tab_new['remainder']=remvals
    tab_new['k']=kvals
    return tab_new

def get_pade_poly(tab,m,n):
    """
    Returns numpy polynomial object for a given entry in tab
    """
    tab_sub = tab[(tab['m'] == m) & (tab['n'] ==n )]
    if type(tab_sub.iloc[0]['p'])==str:
        p = np.poly1d([float(y) for y in tab_sub.iloc[0]['p'].strip('][').split()])
        q = np.poly1d([float(y) for y in tab_sub.iloc[0]['q'].strip('][').split()])

    else:
        p = np.poly1d(tab_sub.iloc[0]['p'])
        q = np.poly1d(tab_sub.iloc[0]['q'])
    return p,q

def calcError(tab,coefs):
    """
    Calculates error between unused coefficients (cumulants) and Talor coefficients
    of Pade approx
    """
    err_vals = []
    coefs_vals = []
    err_next_vals = []
    rel_err_vals = []
    for row in tab.itertuples():
        tc = row[5]
        m = row[1]
        n = row[2]
        err = [a-b for a,b in zip(coefs, tc)]
        err_vals.append(err)
        coefs_vals.append(coefs)
        err_next_vals.append(err[m+n+1])
        rel_err = abs(err[m+n+1])/coefs[m+n+1]
        rel_err_vals.append(rel_err)


    tab_new = tab
    tab_new['coef_input_vals'] = coefs_vals
    tab_new['error_vals_all'] = err_vals
    tab_new['error_next'] = err_next_vals
    tab_new['rel_err'] = rel_err_vals
    return tab_new

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--filename", type=str, help="name of input file", default="spatial_integrals")
    parser.add_argument("--dim", type=int, help="dimension, default is 1", choices=[1, 2], default=1)
    parser.add_argument("--calc_error",action='store_true')
    parser.add_argument("--outname", type=str, help="name for output files (without extension)",
                        default="pade_approx")
    args = parser.parse_args()

    data = pd.read_csv(args.filename+"_dim"+str(args.dim)+".csv")
    sigma_list = data['sigma'].tolist()
    res = pd.DataFrame()
    cols = ["u2_GQ", "u3_GQ"]
    if args.dim==1:
        cols.append("u4_GQ")

    for j in range(len(sigma_list)):
        coefs_sigma = [1]
        sigma = sigma_list[j]
        for i in range(len(cols)): # 3 with u4
            coefs_sigma.append(data.loc[data["sigma"]==sigma,cols].values.tolist()[0][i]*(i+1))
        temp = calc_pade_table(coefs_sigma,args.calc_error)
        if args.calc_error==True:
            temp = calcError(temp,coefs_sigma)
        temp = calc_pole_res(temp)
        temp.insert(loc=0, column='sigma', value=np.repeat(sigma,temp.shape[0]))
        res = pd.concat([res, temp], ignore_index=True, sort=False)
    res.to_csv(args.outname+'_dim'+str(args.dim)+'_error'+str(args.calc_error)+'.csv', index=False)


if __name__ == '__main__':
    main()