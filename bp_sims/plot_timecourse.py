import numpy as np
from matplotlib import pyplot as plt
# ../results/uniform/s0.1_mu1e-09/s0.1_mu1e-09_rho1_L10000_sigma10_time100000.0_r0.1_uniform_rep0.p
L = 10000
rho = 1
N = rho * L**2
mu = 1e-9
print(f"N mu = {N * mu}")
prefix = "results/uniform/s{s}_mu{mu}/s{s}_mu{mu}_rho1_L10000_sigma10_time1000000.0_r0.1_burnin_uniform_{rep}"
for s in [0.1,0.01,0.001]:
    for rep in [f"rep{i}" for i in range(10)] + ["all"]:

        pfile = prefix.format(mu=mu, s=s, rho=rho, L=L, rep=rep) + ".p"
        zfile = prefix.format(mu=mu, s=s, rho=rho, L=L, rep=rep) + ".zero"


        ps = np.loadtxt(pfile)
        zeros = np.loadtxt(zfile)
        n_samples = np.arange(1, len(ps) + 1)
        n_zeros = zeros * n_samples / n_samples[-1]

        if rep == "all":
            print((np.cumsum(ps)[-1] / (n_samples[-1] + n_zeros))[-1])
            # print(f"s = {s}")
            # print(f"mu / s = {mu / s}")
            # print(f"E[P] = {np.sum(ps) / (len(ps) + n_zeros)}")

        plt.figure(figsize=(20,24))
        plt.subplot(311)
        plt.plot(N*ps)
        plt.hlines([1 / s, 2 / s], 0, len(ps), color="black")
        plt.ylabel("$NP$")
        plt.ylim([0, 10/s])
        plt.title(f"s={s}")
        plt.subplot(312)
        plt.plot(np.cumsum(ps)/ (n_samples + n_zeros))
        plt.hlines(mu / s, 0, n_samples[-1], color="black")
        plt.ylabel("running mean $P$")
        plt.subplot(313)
        plt.plot(np.cumsum(ps**2) / (n_samples + n_zeros))
        plt.hlines(mu / s**2 * (mu + 1 / N), 0, n_samples[-1], color="black")
        plt.ylabel("running mean $P^2$")
        plt.xlabel("sample number")
        plt.savefig(f"test_plots/test_s{s}_{rep}.png")

        plt.close()