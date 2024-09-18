##### Script to cut and recapitate trees, with option to change past population time at time x

# import modules
import pyslim as ps
import numpy as np
import tskit
import argparse
import os
import msprime

np.random.seed(1234)
gen_time = 5.6 # from simulations

# Parse command-line arguments
parser = argparse.ArgumentParser(description='Cut, recapitate and mutate trees.')
parser.add_argument('input_file', help='File to be processed')
parser.add_argument('--past_Ne', type=int, default=100, help='Effective population size in past (in coalescent time units)')
parser.add_argument('--output', type=str, help='Output')
parser.add_argument('--time_bottle', type=int, default=1000, help='Number of generations in past at which population size reduced')
args = parser.parse_args()

input_file = args.input_file
past_Ne = args.past_Ne
output = args.output
time_bottle = args.time_bottle

# Get the directory and base filename from the input file path
file_dir = os.path.dirname(input_file)
base_filename = os.path.splitext(os.path.basename(input_file))[0]

# Load tree sequence
slim_ts = tskit.load(input_file)
slim_ts = ps.update(slim_ts)

# compute mutation times
tables = slim_ts.dump_tables()
tables.compute_mutation_times()
tables.sort()
slim_ts = tables.tree_sequence()

# we decapitate
new_ts = slim_ts.decapitate(2000 * gen_time)

# specify ancestry:
N_past = past_Ne

# construct demography
demography = msprime.Demography.from_tree_sequence(slim_ts,
                                                   initial_size = slim_ts.num_individuals)

# adjust population metadata for demographic model
i=0
for pop in demography.populations:
    # must set their effective population sizes
    pop.initial_size = slim_ts.num_individuals * i * gen_time
    pop.initially_active = i
    i+=1

# Population size change from 2000 generations ago to the beginning of the simulation
demography.add_population_parameters_change(time = time_bottle * gen_time, # in past. Scaling needed im order to match SLiM time
                                            initial_size= N_past * gen_time, #
                                            population="pop_1")

# recapitate tree sequence
new_ts = msprime.sim_ancestry(initial_state = new_ts,  # currently this works only with slim_ts but it should be dec_ts
                            recombination_rate = 1e-8 / gen_time,
                              demography=demography) # populations halves in size

# mutate
new_ts = msprime.sim_mutations(new_ts,
                       rate = 1e-10 / gen_time, # fewer mutations per tick since 1 tick ~= 1/5 generations
                       model=msprime.SLiMMutationModel(type=0, slim_generation = 10000), # time ensures all the mutation times match up
                       keep=True,
                       random_seed=7)


# Write out the trees to the specified file
new_ts.dump(output)

