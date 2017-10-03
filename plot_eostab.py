import numpy as np
import matplotlib.pyplot as plt

with open('eostab.dat', 'r', encoding='utf8') as f:

    n_eos_bins_ee, n_eos_bins_r = tuple(map(int, f.readline().split()))
    ln_eem_cgs = np.array(list(map(float, f.readline().split())))
    ln_rm_cgs = np.array(list(map(float, f.readline().split())))
    ln_Pgm_cgs = np.array(list(map(float, f.readline().split()))).reshape((n_eos_bins_ee, n_eos_bins_r), order='F')

fig = plt.figure()
ax = fig.add_subplot(111)

im = ax.imshow(np.log10(np.exp(ln_Pgm_cgs)),
               extent=[np.log10(np.exp(ln_eem_cgs[0])),
                       np.log10(np.exp(ln_eem_cgs[-1])),
                       np.log10(np.exp(ln_rm_cgs[0])),
                       np.log10(np.exp(ln_rm_cgs[-1]))],
               aspect='auto',
               origin='lower')

plt.colorbar(im, ax=ax, label=r'$\log_{10}{\left(P_\mathrm{g}\;[\mathrm{dyn}]\right)}$')

ax.set_xlabel(r'$\log_{10}{\left(e\;[\mathrm{erg}/\mathrm{g}]\right)}$')
ax.set_ylabel(r'$\log_{10}{\left(\rho\;[\mathrm{g}/\mathrm{cm}^3]\right)}$')

ax.set_title('EOS table')

plt.show()
