import numpy as np
from numpy.random import exponential
import argparse

def time_to_next(k,s,mu,rho,r,L):
    return exponential(k*(1-s)+k+mu*rho*L*L+r)

def choose_event(k,s,mu,rho,r,L):
    tot = k*(1-s)+k+mu*rho*L*L+r
    event = np.random.choice(['b','d','m','s'],p=[(k*(1-s)/tot),(k/tot),(mu*rho*L*L/tot),(r/tot)])
    return event

def get_alive(locations):
    return np.where(~np.isnan(locations[:, 0]))[0]

def get_nan(locations):
    return np.where(np.isnan(locations[:, 0]))[0]

def wrap_locations(locations,L):
    wrapped_locations = locations.copy()
    wrapped_locations[:, 0] = (wrapped_locations[:, 0]) % L
    wrapped_locations[:, 1] = (wrapped_locations[:, 1]) % L
    return wrapped_locations

def run_sim_spatial(s, mu, rho, r, sigma, num_iter, max_ind, L=50, sfs_len=100):
    """
    * Carriers appear de novo with rate `mu`*`rho`
    * Carriers give birth (split) with rate `1-s`
    * Carriers die at rate `1`
    * Take samples at rate `r`
    * Habitat is a square (0,L)x(0,L) with periodic boundary conditions

    Implement via Gillespie algorithm. At each step:
    1. Draw waiting time to next event
    2. Update locations of living carriers
    3. Select event type and apply it

    Data structures:
    * During simulation, store individual level-data (continuously updated)
    * At each sample point update the SFS distribution (array w/ pre-specified length)
    * Output SFS distribution
    """

    # initialize array for SFS distribution
    sfs_dist = np.zeros(sfs_len)
    # keep track of time steps

    # keep track of individual level data
    # [x coord, y coord]
    locations = np.full((max_ind, 2),np.nan)
    # initialize current time at 0
    for _ in range(num_iter):
        alive_rows = get_alive(locations)
        k = len(alive_rows)  # number of alive particles
        # draw time to next event
        t_next = time_to_next(k, s, mu, rho, r, L)
        # draw event type
        e_type = choose_event(k, s, mu, rho, r, L)

        ### update spatial coordinates
        if len(alive_rows) > 0:
            locations[alive_rows] += np.random.normal(loc=0, scale=sigma * np.sqrt(t_next), size=(len(alive_rows),2))

        ### mutation
        if e_type == 'm':
            # find next empty row
            empty_row_indices = get_nan(locations)
            if len(empty_row_indices) > 0:  # check that there is a row available
                next_row = empty_row_indices[0]  # choose the first available row
            else:
                locations_new = np.full((locations.shape[0]*2,2),np.nan)
                locations_new[:locations.shape[0],:] = locations
                next_row = locations.shape[0]+1
                locations = locations_new.copy()
            # add row for new lineage at random location
            locations[next_row, :] = [np.random.uniform(0, L), np.random.uniform(0, L)]

        ### death
        elif e_type == 'd':
            ## choose individual who dies
            random_index = np.random.choice(alive_rows)  # choose one at random
            locations[random_index, 0] = np.nan  # mark coord to nan
            locations[random_index, 1] = np.nan  # mark coord to nan

        ### birth
        elif e_type == 'b':
            # find next empty row
            empty_row_indices = get_nan(locations)
            if len(empty_row_indices) > 0:  # check that there is a row available
                next_row = empty_row_indices[0]  # choose the first available row
            else:
                print("Doubling array size")
                locations_new = np.full((locations.shape[0] * 2, 2), np.nan)
                locations_new[:locations.shape[0], :] = locations
                next_row = locations.shape[0] + 1
                locations = locations_new
            # choose parent at random from alive individuals
            random_index = np.random.choice(alive_rows)
            # add row for new (split) lineage with location of parent
            locations[next_row, :] = [locations[random_index, 0], locations[random_index, 1]]

        ### sample NOTE: WILL WANT TO UPDATE TO SPATIAL SAMPLING
        ### currently counts number of extant lineages & updates SFS
        elif e_type == 's':
            if int(k) < sfs_len:
                sfs_dist[int(k)] += 1
                # print for debugging
                if k > 0:
                    print(k)
            else:
                print("Error: SFS entry out of bounds ("+str(k)+")")

        # check if locations out of bounds
        if (np.any(locations < 0)) or (np.any(locations > L)):
           locations = wrap_locations(locations, L)

    return sfs_dist, locations

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s',type=float,help='selection coefficient',default=1e-2)
    parser.add_argument('--mu',type=float,help='mutation rate',default=1e-4)
    parser.add_argument('--dens',type=float,help='population density',default=2)
    parser.add_argument('-r',type=float,help='sampling rate',default=0.1)
    parser.add_argument('--sigma',type=float,help='diffusion coefficient',default=0.2)
    parser.add_argument('--num_iter',type=int,help='number of iterations',default=1000)
    parser.add_argument('--max_ind',type=int,help='max number of individuals',default=1000)
    parser.add_argument('-L',type=float,help='habitat width',default=50)
    parser.add_argument('--sfs_length',type=int,help='max length of sfs to output',default=1000)
    parser.add_argument('--sfs_out',type=str,help='output file name for sfs',default='sfs.csv')
    parser.add_argument('--loc_out',type=str,help='output file name for locations',default='loc.csv')
    parser.add_argument('--seed',type=int,help='random string',default=2024)
    args = parser.parse_args()

    # set seed
    np.random.seed(args.seed)

    # run simulation
    counts, df = run_sim_spatial(s=args.s, mu=args.mu, rho=args.dens, r=args.r, sigma=args.sigma, num_iter=args.num_iter,
                                 max_ind=args.max_ind, L=args.L, sfs_len=args.sfs_length)

    # save output as CSV
    np.savetxt(args.sfs_out,counts,delimiter=',')
    np.savetxt(args.loc_out, df, delimiter=',')

if __name__ == '__main__':
    main()


