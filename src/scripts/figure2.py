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

import numpy as np
from matplotlib import pyplot as plt
import os
import pandas as pd
import seaborn as sns
from utils.constants import *
import utils.plots as pu
import utils.stats as us
import paths

sns.set(context='notebook', palette='colorblind', style='ticks')

# ----------------------------------------------------------------------------
# LOAD DATA

def load_summary_samples(result_paths, all_results=None, force=False,
                         data_key="samples"):
    all_results = all_results or {}
    for run, rundir in result_paths.items():
        if run not in all_results or force:
            path = os.path.join(rundir, "summary_samples.hdf5")
            if os.path.exists(path):
                all_results[run] = pd.read_hdf(path, data_key)
                all_results[run]["run"] = run
            else:
                print(f"WARNING: did not find {path}")
    return all_results

keys = ["gw150914", "ds-o1-fake", "ds-o1-none"]

path_template = os.path.join(paths.data, "figure2/figure2/{}/")
result_paths = {k: path_template.format(k) for k in keys}

all_results = load_summary_samples(result_paths)
all_results_df = pd.concat(all_results.values(), ignore_index=True)

all_results_df[r"$A_0$"] *= A_scale
all_results_df[r"$A_1$"] *= A_scale


kde_grids_path = os.path.join(paths.data, f"fm_t0_kde_grid.h5")
kde_grids_df = pd.read_hdf(kde_grids_path, "data")

injkws = {
    'model': 'mchi_aligned',
    't0': 1126259463.4083147,
    'm': 70.60913274,
    'chi': 0.71379398,
    'cosi': -0.99470868,
    'a': [2.20945342e-21, 3.45682856e-21],
    'phi': [1.20019219, -1.52542676],
    'f': [246.66715929, 241.66408143],
    'tau': [0.00433489, 0.0014342],
}


# ----------------------------------------------------------------------------
# PLOT

ref_key = 'gw150914'
tlim = -5, 7
qs = [0.68]

plot_scale_exp = -20
plot_scale = 10**plot_scale_exp

sr_kws = {
    2048: {'c': sns.color_palette("colorblind")[0]},
    16384: {'c': "gray"},
}

fig, axs = plt.subplots(2, 1, figsize=(pu.fig_width, 4))
for sr, lkws in sr_kws.items():
    df = all_results_df[(all_results_df['run'] == ref_key) &
                        (all_results_df['sr'] == sr)][['t0M', r'$A_1$']].copy()
    df[r'$A_1$'] /= plot_scale
    df = df[['t0M', r'$A_1$']][(df['t0M'] >= tlim[0]) & (df['t0M'] <= tlim[1])]
    l, = axs[0].plot(df.groupby('t0M').median(), label=f'{sr} Hz', **lkws)
    for q in qs:
        hpds = df.groupby('t0M').apply(lambda x: pd.Series(us.compute_hpd(x[r'$A_1$'], q)))
        axs[0].fill_between(hpds.index, hpds[0], hpds[1], color=l.get_color(),
                        alpha=0.1, lw=0)

# add F&M
t0_choices = [-2, 0, 2, 4, 6]
df = kde_grids_df[t0_choices]
m = df.apply(us.get_quantile_from_grid)*A_scale_fm / plot_scale
for q in qs:
    hpd_df = df.apply(us.get_hpd_from_grid, q=q)*A_scale_fm / plot_scale
    axs[0].errorbar(m.index, m, np.abs(hpd_df-m), fmt='.', color='k', alpha=0.5,
                label='Finch & Moore')

# add simulated result
sr = 16384

run_kws = {
    "ds-o1-none": dict(label="no noise", c=sns.color_palette('Set2', desat=0.5)[0]),
    "ds-o1-fake": dict(label="fake noise", c='gray'),
}

for run, lkws in run_kws.items():
    df = all_results_df[(all_results_df['run'] == run) &
                        (all_results_df['sr'] == sr)][['t0M', r'$A_1$']]
    df = df[['t0M', r'$A_1$']][(df['t0M'] >= tlim[0]) & (df['t0M'] <= tlim[1])].copy()
    df[r'$A_1$'] /= plot_scale
    m = df.groupby('t0M').median()
    l, = axs[1].plot(m[m.index >= 0], **{k:v for k,v in lkws.items()})
    axs[1].plot(m[m.index <= 0], alpha=0.5, **{k:v for k,v in lkws.items() if k not in ['alpha', 'label']})
    for q in qs:
        hpds = df.groupby('t0M').apply(lambda x: pd.Series(us.compute_hpd(x[r'$A_1$'], q)))
        t = hpds.index
        axs[1].fill_between(t[t>=0], hpds[t>=0][0], hpds[t>=0][1], color=l.get_color(),
                            alpha=lkws.get('alpha', 0.1), lw=0)
        axs[1].fill_between(t[t<=0], hpds[t<=0][0], hpds[t<=0][1], color=l.get_color(),
                            alpha=lkws.get('alpha', 0.05), lw=0)

t = np.linspace(0, tlim[1])
axs[1].plot(t, A_scale*injkws['a'][1]*np.exp(-t*tM/injkws['tau'][1])/plot_scale, ls='--',
             lw=2, label='truth', zorder=-50, c='k')

# plot incorrect time used in Cotesta et al
axs[0].axvline((1126259462.422828 - 1126259462.423)/tM, ls='--', c=pu.line_color,
                alpha=0.5, zorder=-100)

# plot peak time according to Cotesta et al
axs[0].axvline(0.68, ls='-.', c=pu.line_color, alpha=0.5, zorder=-100)

bbox = dict(boxstyle="round,pad=0.3", fc="white", alpha=0.7, lw=0)
kws = dict(fontsize=16, color='0.2', bbox=bbox)
axs[0].annotate("GW150914", (-5, 8.5E-21/plot_scale), **kws)
axs[1].annotate("simulation", (-5, 8.5E-21/plot_scale), **kws)

axs[0].set_xticklabels([])
plt.subplots_adjust(hspace=0.15)

ax_top = axs[0].twiny()
ax_top.set_xticks(axs[0].get_xticks());
ax_top.set_xlim(axs[0].get_xlim());
ax_top.set_xticklabels(['{:.1f}'.format(t) for t in axs[0].get_xticks()*tM/1E-4]);

ax_in = axs[1].inset_axes([0.05, 0.1, 0.2, 0.3])
t = np.linspace(-80, 80, 300)*tM
h = np.sum([injkws['a'][n]*np.exp(-np.abs(t)/injkws['tau'][n])*np.cos(2*np.pi*injkws['f'][n]*t + injkws['phi'][n])/plot_scale
            for n in range(2)], axis=0)
l, = ax_in.plot(t, h, lw=0.5, ls='-', c='k')
ax_in.plot(t[t>0], h[t>0], lw=1, c=l.get_color())
ax_in.set_axis_off()

axs[1].set_xlabel(r'$(t_0 - t_{\rm ref})/M$');
ax_top.set_xlabel(r'$(t_0 - t_{\rm ref})/10^{-4}\,\mathrm{s}$');
for ax in axs:
    ax.axhline(0, c='gray', ls=':', alpha=0.5, lw=1, zorder=-100)
    ax.axvline(0, c=pu.line_color, alpha=0.5, zorder=-100)
    ax.set_ylabel(r"$A_1\, /\, 10^{{{}}}$".format(plot_scale_exp))
    ax.legend();

plt.subplots_adjust(hspace=0)
for ax in list(axs) + [ax_top]:
    ax.tick_params(direction='in')

figpath = os.path.join(paths.figures, "a1_time_variation.pdf")
plt.savefig(figpath, bbox_inches='tight')
