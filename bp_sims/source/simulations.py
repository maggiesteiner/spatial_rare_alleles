import argparse
import sys
from enum import Enum
from typing import TypeAlias

import numpy as np
from numpy.random import exponential, poisson
from numpy.typing import NDArray
from scipy.stats import binom, truncnorm  # type: ignore

Locations: TypeAlias = NDArray[np.float64]
SFS: TypeAlias = NDArray[np.float64]


class Event(Enum):
    BIRTH = 0
    DEATH = 1
    MUTATION = 2
    SAMPLE = 3


def event_rates(k: int, s: float, theta: float, r: float) -> dict[Event, float]:
    """Calculate the rate of each event type."""
    return {
        Event.BIRTH: k * (1 - s),
        Event.DEATH: k,
        Event.MUTATION: theta,
        Event.SAMPLE: r,
    }


def time_to_next(k: int, s: float, theta: float, r: float) -> float:
    """Generate the waiting time to the next event."""
    if k == 0:
        # If nothing is alive, jump to the next mutation; don't sample
        total_rate = theta
    else:
        total_rate = sum(event_rates(k, s, theta, r).values())
    return exponential(scale=1 / total_rate)


def choose_event(k: int, s: float, theta: float, r: float) -> Event:
    """Choose the type of the next event."""
    # If nothing is alive, jump to the next mutation; don't sample
    if k == 0:
        return Event.MUTATION
    rates = event_rates(k, s, theta, r)
    tot = sum(rates.values())
    choice = np.random.choice(len(Event), p=[rates[e] / tot for e in Event])
    return Event(choice)


def generate_zeros(t_zero: float, r: float) -> int:
    """Simulate the number of samples taken during times when zero carriers are alive"""
    return poisson(t_zero * r)


def get_alive(locations: Locations) -> NDArray[np.intp]:
    return np.where(~np.isnan(locations[:, 0]))[0]


def get_nan(locations: Locations) -> NDArray[np.intp]:
    return np.where(np.isnan(locations[:, 0]))[0]


def get_free_row(locations: Locations) -> int:
    """
    Get the index of the next free row in locations.

    Side effect: If there are no free rows, double the length of locations.
    """
    empty_row_indices = get_nan(locations)
    if len(empty_row_indices) > 0:
        return empty_row_indices[0]
    else:
        sys.stderr.write("No free rows. Doubling array size.\n")
        next_row = locations.shape[0]
        extend_locations(locations)
        return next_row


def extend_locations(locations: Locations) -> None:
    """Double the length of locations in-place."""
    old_length = locations.shape[0]
    locations.resize((2 * old_length, locations.shape[1]), refcheck=False)
    locations[old_length:] = np.nan


def wrap_locations(locations: Locations, L: float) -> Locations:
    return locations % L


def update_locations(
    locations: Locations, sigma: float, t_next: float, L: float
) -> Locations:
    alive_rows = get_alive(locations)
    if len(alive_rows) > 0:
        locations[alive_rows] += np.random.normal(
            loc=0, scale=sigma * np.sqrt(t_next), size=(len(alive_rows), 2)
        )
    locations = wrap_locations(locations, L)
    return locations


def sampling_probability_gaussian(
    locations: Locations, w: float, L: float, rho: float
) -> float:
    locations = locations[~np.isnan(locations).any(axis=1)]
    x1_dens = truncnorm.pdf(
        locations[:, 0], loc=L / 2, scale=w, a=(-L / 2) / w, b=(L / 2) / w
    )
    x2_dens = truncnorm.pdf(
        locations[:, 1], loc=L / 2, scale=w, a=(-L / 2) / w, b=(L / 2) / w
    )
    prod_dens = x1_dens * x2_dens
    return np.sum(prod_dens) / rho


def sample_sfs(
    k: int,
    N: float,
    n: int,
    max_allele_count: int,
    gaussian=None,
    w=None,
    locations=None,
    L=None,
    rho=None,
) -> NDArray[np.float64]:
    sfs_temp = np.zeros(max_allele_count + 1)
    j = np.arange(max_allele_count)
    if gaussian is True:
        sampling_prob = sampling_probability_gaussian(locations, w, L, rho)
        if sampling_prob > 1:
            sys.stderr.write("Warning: p>1, parameters violate model assumptions")
        p = min(sampling_prob, 1)
    else:
        p = k / N
    sfs_temp[:-1] += binom.pmf(j, n, p)  # pmf, entries 0 through max_allele_count-1
    sfs_temp[-1] += binom.sf(max_allele_count - 1, n, p)  # 1 - cdf
    return sfs_temp


def run_sim_spatial(
    n,
    s: float,
    mu: float,
    rho: float,
    r: float,
    sigma: float,
    num_iter: int,
    max_ind: int,
    L: float = 50,
    max_allele_count: int = 100,
    gaussian: bool = False,
    w: float = 1.0,
) -> tuple[Locations, SFS]:
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
    N = rho * (L**2)
    theta = mu * N

    # initialize array for SFS distribution
    running_sfs = np.zeros(max_allele_count + 1)

    # keep track of individual level data
    # [x coord, y coord]
    locations = np.full((max_ind, 2), np.nan)

    # keep a running total of the time with zero carriers alive
    t_zero = 0.0

    # initialize current time at 0
    for _ in range(num_iter):
        alive_rows = get_alive(locations)
        k = len(alive_rows)  # number of alive particles
        # draw time to next event
        t_next = time_to_next(k, s, theta, r)
        if k == 0:
            t_zero += t_next
        # draw event type
        event = choose_event(k, s, theta, r)

        ### update spatial coordinates
        locations = update_locations(locations, sigma, t_next, L)

        if event is Event.MUTATION:
            # add a new lineage at a random location
            next_row = get_free_row(locations)
            locations[next_row] = np.random.uniform(0, L, size=2)

        elif event is Event.DEATH:
            ## choose a random individual to die
            random_index = np.random.choice(alive_rows)
            locations[random_index] = np.nan

        elif event is Event.BIRTH:
            # choose parent at random from alive individuals
            parent_index = np.random.choice(alive_rows)
            # add new lineage with location of parent
            next_row = get_free_row(locations)
            locations[next_row] = locations[parent_index]

        elif event is Event.SAMPLE:
            running_sfs += sample_sfs(
                k, N, n, max_allele_count, gaussian, w, locations, L, rho
            )

    # Simulate the zero count SFS bin
    running_sfs[0] += generate_zeros(t_zero, r)

    # Normalize expected SFS to one
    expected_sfs = running_sfs / np.sum(running_sfs)

    return expected_sfs, locations


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-n", type=int, help="sample size", default=1e3)
    parser.add_argument("-s", type=float, help="selection coefficient", default=1e-2)
    parser.add_argument("--mu", type=float, help="mutation rate", default=1e-4)
    parser.add_argument("--dens", type=float, help="population density", default=2)
    parser.add_argument("-r", type=float, help="sampling rate", default=0.1)
    parser.add_argument(
        "--sigma", type=float, help="diffusion coefficient", default=0.2
    )
    parser.add_argument(
        "--num_iter", type=int, help="number of iterations", default=1000
    )
    parser.add_argument(
        "--max_ind", type=int, help="max number of individuals", default=1000
    )
    parser.add_argument("-L", type=float, help="habitat width", default=50)
    parser.add_argument(
        "--max_allele_count", type=int, help="max allele count to track", default=1000
    )
    parser.add_argument(
        "--sfs_out", type=str, help="output file name for sfs", default="sfs.csv"
    )
    parser.add_argument(
        "--loc_out", type=str, help="output file name for locations", default="loc.csv"
    )
    parser.add_argument("--seed", type=int, help="random string", default=2024)
    parser.add_argument(
        "--gaussian",
        action="store_true",
        help="implement Gaussian sampling kernel",
        default=False,
    )
    parser.add_argument("-w", type=float, help="width for sampling kernel", default=1)
    args = parser.parse_args()

    # set seed
    np.random.seed(args.seed)

    # run simulation
    counts, df = run_sim_spatial(
        n=args.n,
        s=args.s,
        mu=args.mu,
        rho=args.dens,
        r=args.r,
        sigma=args.sigma,
        num_iter=args.num_iter,
        max_ind=args.max_ind,
        L=args.L,
        max_allele_count=args.max_allele_count,
        gaussian=args.gaussian,
        w=args.w,
    )

    # save output as CSV
    np.savetxt(args.sfs_out, counts, delimiter=",")
    np.savetxt(args.loc_out, df, delimiter=",")


if __name__ == "__main__":
    main()
