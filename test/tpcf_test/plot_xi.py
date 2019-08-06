import matplotlib
matplotlib.use('TkAgg')
import numpy as np
import matplotlib.pylab as plt
import iiiraven2 as ir2
ir2.set_plot_default()

data = np.loadtxt('xi.txt', unpack=True)
r = data[0]
xi = data[4]

data = np.loadtxt('xi_weighted.txt', unpack=True)
rw = data[0]
xiw = data[4]

plt.figure(figsize=(6, 6))
plt.plot(r, xi)
plt.plot(rw, xiw)
plt.xlabel(r'$r$', fontsize=16)
plt.ylabel(r'$\xi(r)$', fontsize=16)
plt.xscale('log')
plt.yscale('symlog')
plt.tight_layout()
plt.savefig('plot_r_vs_xi.pdf')
plt.close()
