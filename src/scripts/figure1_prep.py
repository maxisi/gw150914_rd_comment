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

import pickle as pkl
import os
import arviz as az
import seaborn as sns
from utils.constants import *
import paths
from tqdm import tqdm

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
# indicate such chains in a list per key here, if any
bad_chains = {}

# ----------------------------------------------------------------------------
# EXPORT DATA

samples = {}
for k, r in tqdm(results.items()):
    c = [i for i in range(r.posterior.dims['chain'])
         if i not in bad_chains.get(k, [])]
    if 'A1' in r.posterior:
        x = A_scale_pr*r.posterior.A1[c,:].values.flatten()
    else:
        x = A_scale*r.posterior.A[c,:,1].values.flatten()
    samples[k] = x

fname = paths.data / 'fig1_samples_dict.pkl'
with open(fname, 'wb') as f:
    pkl.dump(samples, f, protocol=-1)
print(f"Saved: {fname}")
