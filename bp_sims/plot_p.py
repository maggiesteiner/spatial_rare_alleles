import numpy as np
import matplotlib.pyplot as plt

N = 1e5
mu = 1e-8
prefix = "results/uniform/s{s}_mu{mu}/s{s}_mu{mu}_rho100000.0_L1.0_sigma10_iter1000000_r0.1_uniform_rep{rep}"

for s in [0.1, 0.01, 0.001]:
    for rep in range(10):
        pfile = prefix.format(mu=mu, s=s, rep=rep) + ".p"
        zfile = prefix.format(mu=mu, s=s, rep=rep) + ".zero"

        ps = np.loadtxt(pfile)
        zeros = np.loadtxt(zfile)

        plt.figure()
        plt.subplot(311)
        plt.plot(ps)
        plt.ylabel("$P$")
        plt.title(f"s={s}")
        plt.subplot(312)
        n_samples = np.arange(1, len(ps) + 1)
        n_zeros = zeros * n_samples / n_samples[-1]
        plt.plot(np.cumsum(ps) / (n_samples + n_zeros))
        plt.hlines(mu / s, 0, n_samples[-1], color="0.5")
        plt.ylabel("running mean $P$")
        plt.subplot(313)
        plt.plot(np.cumsum(ps**2) / (n_samples + n_zeros))
        plt.hlines(mu / s**2 * (mu + 1 / N), 0, n_samples[-1], color="0.5")
        plt.ylabel("running mean $P^2$")
        plt.xlabel("sample number")
        plt.savefig(f"example_s{s}_rep{rep}.png")
        plt.close()
