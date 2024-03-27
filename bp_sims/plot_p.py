import numpy as np
import matplotlib.pyplot as plt

mu = 1e-8
prefix = "results/uniform/s{s}_mu{mu}/s{s}_mu{mu}_rho100000.0_L1.0_sigma10_iter100000_r0.1_uniform_all"

for s in [0.1, 0.01, 0.001]:
    pfile = prefix.format(mu=mu, s=s) + ".p"
    zfile = prefix.format(mu=mu, s=s) + ".zero"

    ps = np.loadtxt(pfile)
    zeros = np.loadtxt(zfile)
    ep1 = np.sum(ps) / (len(ps) + zeros)
    print(ep1)
    print(mu / s)

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
    plt.ylabel("running mean $P^2$")
    plt.xlabel("sample number")
    plt.savefig(f"example_s={s}.png")
