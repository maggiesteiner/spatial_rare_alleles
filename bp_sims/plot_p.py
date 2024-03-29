import numpy as np
import matplotlib.pyplot as plt

L = 10000
rho = 1
N = rho * L**2
mu = 1e-9
print(f"N mu = {N * mu}")
prefix = "results/uniform/s{s}_mu{mu}/s{s}_mu{mu}_rho{rho}_L{L}_sigma10_time50000.0_r0.1_uniform_{rep}"

for s in [0.1, 0.01, 0.001]:
    for rep in [f"rep{i}" for i in range(10)] + ["all"]:
        pfile = prefix.format(mu=mu, s=s, rho=rho, L=L, rep=rep) + ".p"
        zfile = prefix.format(mu=mu, s=s, rho=rho, L=L, rep=rep) + ".zero"

        ps = np.loadtxt(pfile)
        zeros = np.loadtxt(zfile)
        n_samples = np.arange(1, len(ps) + 1)
        n_zeros = zeros * n_samples / n_samples[-1]

        if rep == "all":
            print(f"s = {s}")
            print(f"mu / s = {mu / s}")
            print(f"E[P] = {np.sum(ps) / (len(ps) + n_zeros)}")

        plt.figure()
        plt.subplot(311)
        plt.plot(N*ps)
        plt.hlines([1 / s, 2 / s], 0, len(ps), color="0.5")
        plt.ylabel("$NP$")
        plt.ylim([0, 3/s])
        plt.title(f"s={s}")
        plt.subplot(312)
        plt.plot(np.cumsum(ps) / (n_samples + n_zeros))
        plt.hlines(mu / s, 0, n_samples[-1], color="0.5")
        plt.ylabel("running mean $P$")
        plt.subplot(313)
        plt.plot(np.cumsum(ps**2) / (n_samples + n_zeros))
        plt.hlines(mu / s**2 * (mu + 1 / N), 0, n_samples[-1], color="0.5")
        plt.ylabel("running mean $P^2$")
        plt.xlabel("sample number")
        plt.savefig(f"example_plots/example_s{s}_{rep}.png")
        plt.close()
