#! /usr/bin/env python

# -*- coding: utf-8 -*-
#
#       Copyright 2023
#       Maximiliano Isi <max.isi@ligo.org>
#       Will M. Farr <will.farr@ligo.org>
#
#       This program is free software; you can redistribute it and/or modify
#       it under the terms of the GNU General Public License as published by
#       the Free Software Foundation; either version 2 of the License, or
#       (at your option) any later version.
#
#       This program is distributed in the hope that it will be useful,
#       but WITHOUT ANY WARRANTY; without even the implied warranty of
#       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#       GNU General Public License for more details.
#
#       You should have received a copy of the GNU General Public License
#       along with this program; if not, write to the Free Software
#       Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#       MA 02110-1301, USA.

from matplotlib import pyplot as plt
import os
import arviz as az
import seaborn as sns
from utils.constants import *
import utils.plots as pu
import paths

sns.set(context='notebook', palette='colorblind', style='ticks')

# ----------------------------------------------------------------------------
# LOAD DATA

keys = [
    "pyring_t0-good_timestamps-good_seglen-good",
    "ringdown_t0-good_timestamps-good_seglen-good",
    "pyring_t0-bad_timestamps-bad_seglen-bad",
    "ringdown_t0-bad_timestamps-good_seglen-bad",
    "ringdown_t0-bad_timestamps-good_seglen-good",
]

path_template = os.path.join(paths.data, "figure1/figure1/{}/result.nc")
results = {k: az.from_netcdf(path_template.format(k)) for k in keys}

# sometime HMC chains get stuck and must be removed from the result;
# indicate such chains here, if any
bad_chains = {
    "ringdown_t0-bad_timestamps-good_seglen-good": [2],
}

# ----------------------------------------------------------------------------
# PLOT

plot_scale_exp = -20

def plot_curves(plot_kws, ax):
    for k, kws in plot_kws.items():
        r = results[k]
        c = [i for i in range(r.posterior.dims['chain'])
             if i not in bad_chains.get(k, [])]
        if 'A1' in r.posterior:
            x = A_scale_pr*r.posterior.A1[c,:].values.flatten()
        else:
            x = A_scale*r.posterior.A[c,:,1].values.flatten()
        sns.histplot(x/10**plot_scale_exp, element='step', fill=False,
                     stat='density', ax=ax, **kws)

fig, (ax, axins) = plt.subplots(2, figsize=pu.figsize_column)
plot_kws = {
    "pyring_t0-bad_timestamps-bad_seglen-bad":
        dict(alpha=0.8, color='k', ls=':',
             label='Cotesta+2022 "$-0.72M$"\n'+r'($t_0 = t_{\rm ref} - 0.5M,T=0.1\,{\rm s}$)'),
    "ringdown_t0-good_timestamps-good_seglen-good":
        dict(lw=2, color=sns.color_palette()[0],
             label='ringdown\n'+r'($t_0 = t_{\rm ref},T=0.2\,{\rm s}$)'),
    "pyring_t0-good_timestamps-good_seglen-good":
        dict(lw=2, color=sns.color_palette()[0], ls='--',
             label='PyRing\n'+r'($t_0 = t_{\rm ref},T=0.2\,{\rm s}$)'),
}
plot_curves(plot_kws, ax)

inset_kws = {
    "pyring_t0-bad_timestamps-bad_seglen-bad":
        dict(alpha=0.8, ls=':', color='k'),
    "ringdown_t0-bad_timestamps-good_seglen-bad":
        dict(alpha=0.5, lw=1.5, color=sns.color_palette(desat=0.25)[2],
             label='fixed time-stamps\n'+r'($t_0 = t_{\rm ref} - 0.5M,T=0.1\,{\rm s}$)'),
    "ringdown_t0-bad_timestamps-good_seglen-good":
        dict(lw=2, ls='-', color=sns.color_palette(desat=1)[2],
             label='fixed time-stamps & duration\n'+r'($t_0 = t_{\rm ref} - 0.5M,T=0.2\,{\rm s}$)'),
}
plot_curves(inset_kws, axins)

plt.subplots_adjust(hspace=0)
fig.canvas.draw()
yl = [tick.get_text() for tick in axins.get_yticklabels()]
axins.set_yticklabels(yl[:-2])
axins.tick_params(direction='in')
ax.tick_params(direction='in', labelbottom=False)
axins.set_xlabel(r'$A_1\, /\, 10^{{{}}}$'.format(plot_scale_exp))

for a in [ax, axins]:
    a.legend()
    a.set_xlim(-0.02, 1.75)
    a.set_ylim(0, 2)
figpath = os.path.join(paths.figures, "a1_posterior_at_reference_time.pdf")
plt.savefig(figpath, bbox_inches='tight')
