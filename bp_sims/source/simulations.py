import numpy as np
from numpy.random import exponential
import pandas as pd
import argparse

def time_to_next(k,s,mu,rho,r):
    return exponential(k*(1-s)+k+mu*rho+r)

def choose_event(k,s,mu,rho,r):
    tot = k*(1-s)+k+mu*rho+r
    event = np.random.choice(['b','d','m','s'],p=[(k*(1-s)/tot),(k/tot),(mu*rho/tot),(r/tot)])
    return event


def run_sim_spatial(s, mu, rho, r, sigma, num_iter, max_ind, L=50, sfs_len=100):
    # initialize array for SFS distribution
    sfs_dist = np.zeros(sfs_len)
    # keep track of time steps
    counter = 1  # initialize at 1
    # keep track of individual level data
    # [alive/dead, x coord, y coord, time @ birth, time @ death]
    ind_data = np.zeros((max_ind, 5))
    # initialize current time at 0
    curr_time = 0
    while counter <= num_iter:
        k = np.sum(ind_data[:, 0])  # number of alive particles
        # draw time to next event
        t_next = time_to_next(k, s, mu, rho, r)
        # update time
        curr_time += t_next
        # draw event type
        e_type = choose_event(k, s, mu, rho, r)

        ### update spatial coordinates
        indices = np.where(ind_data[:, 0] == 1)[0]
        if len(indices) > 0:
            ind_data[indices, 1] += np.random.normal(loc=0, scale=0.2 * t_next, size=len(indices))  # update x
            ind_data[indices, 2] += np.random.normal(loc=0, scale=0.2 * t_next, size=len(indices))  # update y

        ### mutation
        if e_type == 'm':
            # find next empty row
            empty_row_indices = np.where(np.all(ind_data == 0, axis=1))[0]
            if len(empty_row_indices) > 0:  # check that there is a row available
                next_row = empty_row_indices[0]  # choose the first available row
            else:
                print("ERROR: ran out of room in array!")  # print error message and stop
                break
            # add row for new lineage at random location
            ind_data[next_row, :] = [1, np.random.uniform(0, L), np.random.uniform(0, L), curr_time, 0]

        ### death
        elif e_type == 'd':
            ## choose individual who dies
            indices = np.where(ind_data[:, 0] == 1)  # all alive individuals
            random_index = np.random.choice(indices[0])  # choose one at random
            ind_data[random_index, 0] = 0  # mark first column to dead
            ind_data[random_index, 4] = curr_time  # mark time of death

        ### birth
        elif e_type == 'b':
            # find next empty row
            empty_row_indices = np.where(np.all(ind_data == 0, axis=1))[0]
            if len(empty_row_indices) > 0:  # check that there is a row available
                next_row = empty_row_indices[0]  # choose the first available row
            else:
                print("ERROR: ran out of room in array!")  # print error message and stop
                break
            # choose parent at random from alive individuals
            indices = np.where(ind_data[:, 0] == 1)
            random_index = np.random.choice(indices[0])
            # add row for new (split) lineage with location of parent
            ind_data[next_row, :] = [1, ind_data[random_index, 1], ind_data[random_index, 2], curr_time, 0]

        ### sample NOTE: WILL WANT TO UPDATE TO SPATIAL SAMPLING
        ### currently counts number of extant lineages & updates SFS
        elif e_type == 's':
            if int(k) < sfs_len:
                sfs_dist[int(k)] += 1
                # print for debugging
                if k > 0:
                    print(k)
            else:
                print("Error: SFS entry out of bounds")

        # update counter
        counter += 1

    return sfs_dist, ind_data

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
    args = parser.parse_args()

    counts, df = run_sim_spatial(s=args.s, mu=args.mu, rho=args.dens, r=args.r, sigma=args.sigma, num_iter=args.num_iter,
                                 max_ind=args.max_ind, L=args.L, sfs_len=args.sfs_length)
    print(counts[0:10])
if __name__ == '__main__':
    main()


