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
from scipy.stats import gaussian_kde

def compute_hpd(samples, p=0.68, out="both"):
    # compute the minimum-width p%-credible interval 
    sorted_samples = np.sort(samples)
    n = len(samples)
    # number of samples that will be enclosed in p%-credible interval
    n_in = int(np.floor(p*n))
    # number of samples that will be outside p%-credible interval
    # this is also how many locations there are for the first sample in the 
    # interval, since we can slide the `n_in` long window over `n_out` slots
    n_out = n - n_in
    # compute p%-credible interval widths for different starting locations
    # the first term `sorted_samples[n_in]` is the lowest-possible location
    # of the high-end of the CI; the second term `sorted_samples[n_out]` is
    # the highest-possible location of the low_end of the CI; others are in
    # between with steps of 1
    widths = sorted_samples[n_in-1:-1] - sorted_samples[0:n_out]
    # find location of first sample in tightest CI
    i = np.argmin(widths)
    if out.lower() == "both":
        return sorted_samples[i], sorted_samples[i+n_in]
    elif out.lower() == "high":
        return sorted_samples[i+n_in]
    elif out.lower() == "low":
        return sorted_samples[i]
    
def get_hpd_from_grid(x, q=0.68, A_min=0):
    xs = x.sort_values(ascending=False)
    d = x.index[1] - x.index[0]
    lebesgue_integral = np.cumsum(xs)*d
    # normalize it
    lebesgue_integral /= max(lebesgue_integral)
    l, h = np.sort(np.abs(lebesgue_integral-q).sort_values().index.values[:2])
    if np.abs(l-h) < 2*d:
        # the interval has effectively zero width
        # assume that's because we've hit an edge on the LHS and set
        # that lower value to the minimum
        l = A_min
    #print(q, lebesgue_integral[l], lebesgue_integral[h])
    return l, h

def get_quantile_from_grid(x, q=0.5, return_q=False):
    xs = x.sort_index(ascending=True)
    d = x.index[1] - x.index[0]
    cdf = np.cumsum(xs)*d
    cdf /= max(cdf)
    i = (np.abs(cdf - q)).idxmin()
    if np.abs(cdf[i] - q) > 0.01:
        print(f"WARNING: quantile error greater than 1% ({cdf[i]})")
    if return_q:
        return i, cdf[i]
    else:
        return i
    
