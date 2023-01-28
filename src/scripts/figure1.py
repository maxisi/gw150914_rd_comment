#! /usr/bin/env python

from matplotlib import pyplot as plt

import os
import paths

import arviz as az
import seaborn as sns

from utils.constants import *
import utils.plots as pu
import utils.io as io
import utils.stats as us

sns.set(context='notebook', palette='colorblind', style='ticks')

keys = [
    "pyring_t0-good_timestamps-good_seglen-good",
    "ringdown_t0-good_timestamps-good_seglen-good",
    "pyring_t0-bad_timestamps-bad_seglen-bad",
    "ringdown_t0-bad_timestamps-good_seglen-bad",
    "ringdown_t0-bad_timestamps-good_seglen-good",
]

path_template = os.path.join(paths.data, "figure1/{}/result.nc")
results = {k: az.from_netcdf(path_template.format(k)) for k in keys}

bad_chains = {
    "ringdown_t0-bad_timestamps-good_seglen-good": [2],
}

plot_scale_exp = -20

fig, (ax, axins) = plt.subplots(2, figsize=pu.figsize_column)
plot_kws = {
    "pyring_t0-bad_timestamps-bad_seglen-bad": dict(alpha=0.8, color='k', ls=':', label=r'Cotesta+2022'+' "$-0.72M$"\n'+r'($t_0 = t_{\rm ref} - 0.5M,T=0.1\,{\rm s}$)'),
    "ringdown_t0-good_timestamps-good_seglen-good": dict(lw=2, color=sns.color_palette()[0], label=r'ringdown'+'\n'+r'($t_0 = t_{\rm ref},T=0.2\,{\rm s}$)'),
    "pyring_t0-good_timestamps-good_seglen-good": dict(lw=2, color=sns.color_palette()[0], ls='--', label=r'PyRing'+'\n'+r'($t_0 = t_{\rm ref},T=0.2\,{\rm s}$)'),
}

for k, kws in plot_kws.items():
    r = results[k]
    c = [i for i in range(r.posterior.dims['chain']) if i not in bad_chains.get(k, [])]
    if 'A1' in r.posterior:
        x = A_scale_pr*r.posterior.A1[c,:].values.flatten()
    else:
        x = A_scale*r.posterior.A[c,:,1].values.flatten()
    sns.histplot(x/10**plot_scale_exp, element='step', fill=False,
                 stat='density', ax=ax, **kws);

inset_kws = {
    "pyring_t0-bad_timestamps-bad_seglen-bad": dict(alpha=0.8, ls=':', color='k'),
    "ringdown_t0-bad_timestamps-good_seglen-bad": dict(alpha=0.5, lw=1.5, color=sns.color_palette(desat=0.25)[2], label='fixed time-stamps\n'+r'($t_0 = t_{\rm ref} - 0.5M,T=0.1\,{\rm s}$)'),
    "ringdown_t0-bad_timestamps-good_seglen-good": dict(lw=2, ls='-', color=sns.color_palette(desat=1)[2], label='fixed time-stamps & duration\n'+r'($t_0 = t_{\rm ref} - 0.5M,T=0.2\,{\rm s}$)'),
}

for k, kws in inset_kws.items():
    r = results[k]
    c = [i for i in range(r.posterior.dims['chain']) if i not in bad_chains.get(k, [])]
    if 'A1' in r.posterior:
        x = A_scale_pr*r.posterior.A1[c,:].values.flatten()
    else:
        x = A_scale*r.posterior.A[c,:,1].values.flatten()
    sns.histplot(x/10**plot_scale_exp, element='step', fill=False,
                 stat='density', ax=axins, **kws);

plt.subplots_adjust(hspace=0)
ax.tick_params(direction='in', labelbottom=False)
axins.set_yticks(axins.get_yticks()[:-2])
axins.set_xlabel(r'$A_1\, /\, 10^{{{}}}$'.format(plot_scale_exp))

for a in [ax, axins]:
    a.legend()
    a.set_xlim(-0.02, 1.75)
    a.set_ylim(0, 2)
figpath = os.path.join(paths.figures, "a1_posterior_at_reference_time.pdf")
plt.savefig(figpath, bbox_inches='tight')
