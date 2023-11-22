#!python3

"""Subsampling allele frequencies across populations from large catalogs."""

import numpy as np 
import pandas as pd
import click 


def resample_alleles(ac, an, min_n = 100, n=5000, seed=42):
    """Resample alleles for a given sample-size.
    
    We exclude variants where there are fewer than min_n called genotypes...
    """
    assert ac.size == an.size
    assert min_n > 1
    assert seed > 0   
    np.random.seed(seed)
    af = np.nan_to_num(ac/an)
    ac_resamp = np.random.binomial(n=n, p=af)
    af_resamp = ac_resamp / n
    ac_resamp[an < min_n] = 0
    af_resamp[an < min_n] = 0
    return ac_resamp, af_resamp


def resamp_alleles_multipop(acs, ans, props=np.array([1.0, 0.0, 0.0, 0.0, 0.0]), n=5000, seed=42):
    """Resample multipopulation alleles ..."""
    assert acs.shape[0] == ans.shape[0]
    assert acs.shape[1] == ans.shape[1]
    assert props.size == acs.shape[1]
    assert np.all(props >= 0.0)
    if ~np.isclose(np.sum(props), 1.0):
        # Rescale if necessary here ... 
        props = props / np.sum(props)

    joint_acs = np.zeros(shape=acs.shape, dtype=np.uint16)
    joint_afs = np.zeros(shape=acs.shape, dtype=np.float16)
    new_ns = np.zeros(props.size, dtype=np.uint32)
    for i in range(props.size):
        # Iterate through the populations now 
        cur_n = int(n*props[i])
        new_ns[i] = cur_n
        if cur_n < 1:
            pass
        else:
            assert cur_n > 0
            ac_pop, af_pop = resample_alleles(acs[:,i], ans[:,i], n=cur_n, seed=seed)
            joint_acs[:,i] = ac_pop
            joint_afs[:,i] = af_pop
    meta_acs = np.sum(joint_acs, axis=1)
    meta_afs = meta_acs / n
    return meta_acs, meta_afs, joint_acs, joint_afs, new_ns


# @click.command()
# @click.option(
#     "--sfstsv", "-i", required=True, type=str, help="TSV File detailing SFS."
# )
# @click.option(
#     "--poplist", "-s", required=True, type=str, help="Comma separated list of population definitions."
# )
# @click.option(
#     "--proportions", "-p", required=True, type=str, help="Sampling proportions for each population for allele freqs."
# )
# @click.option(
#     "--seed", required=False, type=int, default=42, help="Random seed for subsampling alleles."
# )
# @click.option(
#     "--n", "-n", required=False, type=int, default=5000, help="Number of samples to subsample."
# )
# @click.option(
#     "--out",
#     "-o",
#     required=True,
#     type=str,
#     default="karyohmm",
#     help="Output file prefix.",
# )
# def main(sfs_tsv, poplist, proportions, n=5000, seed=42, out="test.subsamp.sfs.tsv.gz"):
#     raise NotImplementedError("This function is not implemented yet!")
#     pops = poplist.split(',')
#     props = [float(x) for x in proportions.split(',')]
#     sfs_df = pd.read_csv(sfs_tsv, sep="\t")
#     # Any variants not called in a single population will be dropped ...
#     sfs_df[sfs_df == '.'] = np.nan
#     an_pops = [f'AN_{p}' for p in pops]
#     ac_pops = [f'AC_{p}' for p in pops]
#     sfs_df.dropna(subset=an_pops, inplace=True)
#     joint_acs = sfs_df[ac_pops].astype(int).values
#     joint_ans = sfs_df[an_pops].astype(int).values
#     max_pop_n = joint_ans.max(axis=0)
#     print(max_pop_n)
#     subsamp_acs, subsamp_afs, _, _, ns1 = resamp_alleles_multipop(acs=joint_acs, ans=joint_ans, props=props, n=n)
#
#     subsamp_sfs_df = sfs_df[['Annot', 'Effect']].iloc[subsamp_acs > 0,:]
#     subsamp_sfs_df['AC'] = subsamp_acs[subsamp_acs > 0]
#     subsamp_sfs_df['AF'] = subsamp_afs[subsamp_acs > 0]
#     subsamp_sfs_df['N'] = n
#     subsamp_sfs_df.to_csv(out, sep="\t", index=None)


if __name__ == "__main__":
    try:
        sfs_df = pd.read_csv(snakemake.input['input_jsfs'], sep="\t")
        print("Finished Reading input SFS!")
        # Any variants not called in a single population will be dropped ... 
        sfs_df[sfs_df == '.'] = np.nan
        poplist = snakemake.params['poplist']
        an_pops = [f'AN_{p}' for p in poplist]
        ac_pops = [f'AC_{p}' for p in poplist]
        pop_props = snakemake.params["props"]
        sfs_df.dropna(subset=an_pops, inplace=True)
        # Obtain the joint allele counts for subsampling ... 
        joint_acs = sfs_df[ac_pops].astype(int).values
        joint_ans = sfs_df[an_pops].astype(int).values
        max_pop_n = joint_ans.max(axis=0)
        print("Finished subsetting SFS!")
        print(max_pop_n)
        subsamp_acs, subsamp_afs, _, _, ns1 = resamp_alleles_multipop(acs=joint_acs, ans=joint_ans, props=pop_props, n=int(snakemake.params['n']), seed=int(snakemake.wildcards['seed']))
        
        #if 'barts' in poplist:
        #    subsamp_sfs_df = sfs_df['Annot']#.iloc[subsamp_acs > 0,:]
        #else:
        subsamp_sfs_df = sfs_df[['Annot', 'Effect']]#.iloc[subsamp_acs > 0,:]
        subsamp_sfs_df['AC'] = subsamp_acs#[subsamp_acs > 0]
        subsamp_sfs_df['AF'] = subsamp_afs#[subsamp_acs > 0]
        subsamp_sfs_df['CHROM'] = sfs_df['CHROM']
        subsamp_sfs_df['POS'] = sfs_df['POS']
        subsamp_sfs_df['N'] = int(snakemake.params['n']) 
        subsamp_sfs_df.to_csv(snakemake.output['subsamp_sfs_tsv'], sep="\t", index=None)
    except:
        main()
