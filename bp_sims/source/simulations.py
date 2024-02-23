import numpy as np
from numpy.random import exponential, poisson
import argparse

def time_to_next(k: int, s: float, theta: float, r: float) -> float:
    """Generate the waiting time to the next event."""
    if k == 0:
        # If nothing is alive, jump to the next mutation; don't sample
        total_rate = theta
    else:
        total_rate = k * (1 - s) + k + theta + r
    return exponential(scale=1 / total_rate)

def choose_event(k: int, s: float, theta: float, r: float) -> str:
    """Choose the type of the next event."""
    # If nothing is alive, jump to the next mutation; don't sample
    if k == 0:
        return "m"
    tot = k*(1-s)+k+theta+r
    event = np.random.choice(['b','d','m','s'],p=[(k*(1-s)/tot),(k/tot),(theta/tot),(r/tot)])
    return event

def generate_zeros(t_zero: float, r: float) -> int:
    """Simulate the number of samples taken during times when zero carriers are alive"""
    return poisson(t_zero * r)

def get_alive(locations):
    return np.where(~np.isnan(locations[:, 0]))[0]

def get_nan(locations):
    return np.where(np.isnan(locations[:, 0]))[0]

def wrap_locations(locations,L):
    return locations % L

def extend_locations(locations):
    return np.pad(locations, pad_width=((0, locations.shape[0]), (0, 0)), constant_values=np.nan)

def update_locations(locations,sigma,t_next,L):
    alive_rows = get_alive(locations)
    if len(alive_rows) > 0:
        locations[alive_rows] += np.random.normal(loc=0, scale=sigma * np.sqrt(t_next), size=(len(alive_rows), 2))
    locations = wrap_locations(locations, L)
    return locations

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
    * During simulation, store locations of individuals at the current time
    * At each sample point update the SFS distribution (array w/ pre-specified length)
    * Output SFS distribution
    """

    # parameter for total mutation rate
    theta = mu*rho*(L**2)

    # initialize array for SFS distribution
    sfs_dist = np.zeros(sfs_len)
    # keep track of time steps

    # keep track of individual level data
    # [x coord, y coord]
    locations = np.full((max_ind, 2),np.nan)

    # keep a running total of the time with zero carriers alive
    t_zero = 0

    # initialize current time at 0
    for _ in range(num_iter):
        alive_rows = get_alive(locations)
        k = len(alive_rows)  # number of alive particles
        # draw time to next event
        t_next = time_to_next(k,s,theta,r)
        if k == 0:
            t_zero += t_next
        # draw event type
        e_type = choose_event(k,s,theta,r)

        ### update spatial coordinates
        locations = update_locations(locations,sigma,t_next,L)

        ### mutation
        if e_type == 'm':
            # find next empty row
            empty_row_indices = get_nan(locations)
            if len(empty_row_indices) > 0:  # check that there is a row available
                next_row = empty_row_indices[0]  # choose the first available row
            else:
                print("Doubling array size")
                next_row = locations.shape[0] + 1
                locations = extend_locations(locations)

            # add row for new lineage at random location
            locations[next_row] = np.random.uniform(0, L, size=2)

        ### death
        elif e_type == 'd':
            ## choose individual who dies
            random_index = np.random.choice(alive_rows)  # choose one at random
            locations[random_index] = np.nan  # mark coords to nan

        ### birth
        elif e_type == 'b':
            # find next empty row
            empty_row_indices = get_nan(locations)
            if len(empty_row_indices) > 0:  # check that there is a row available
                next_row = empty_row_indices[0]  # choose the first available row
            else:
                print("Doubling array size")
                next_row = locations.shape[0] + 1
                locations = extend_locations(locations)

            # choose parent at random from alive individuals
            parent_index = np.random.choice(alive_rows)
            # add row for new (split) lineage with location of parent
            locations[next_row] = locations[parent_index]

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

    # Simulate the zero count SFS bin
    sfs_dist[0] += generate_zeros(t_zero, r)

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


