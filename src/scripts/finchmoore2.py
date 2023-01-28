#! /usr/bin/env python
#
# The following code is adapted directly from Finch & Moore (2022)
#
# overtone_amplitude_at_tref.ipynb
# https://zenodo.org/record/6949492#.Y9WaTC-B2Zw
#
# which was released under the following public license:
# Creative Commons Attribution 4.0 International
#
# Adapted by: Max Isi (max.isi@ligo.org) [Jan, 2023]

import paths
import os
import pandas as pd
from tqdm import tqdm
import multiprocessing as mp

import numpy as np
from scipy.stats import gaussian_kde

import qnm
import scipy.constants as consts
G, c = consts.G, consts.c
Msun = 1.9884e30

# GW150914
Mf_det = 68.779
conversion = Mf_det*Msun*G/c**3

t_ref = 1126259462.423/conversion

# We convert all times to be in units of Mf

# Overtone analysis
# -----------------

usecols = ['t_0 - t_ref [s] [Hanford]', 'A_rd_1', 'M_f [M_sun]', 'chi_f']

N1_dir = os.path.join(paths.data, 'finch_moore/posterior_samples/GW150914/3W220221')
posterior = pd.read_csv(f'{N1_dir}/posterior_samples.dat', usecols=usecols)
t0_Hanford_posterior = posterior['t_0 - t_ref [s] [Hanford]'].to_numpy()/conversion
A1_posterior = posterior['A_rd_1'].to_numpy()*1e21
Mf_posterior = posterior['M_f [M_sun]'].to_numpy()
chif_posterior = posterior['chi_f'].to_numpy()

qnm_func = qnm.modes_cache(-2, 2, 2, 1)

A1_posterior_tref = []

A1_ref_cache = os.path.join(paths.data, 'A1_posterior_tref.dat')
A1_posterior_tref = np.loadtxt(A1_ref_cache)

def gaussian(x, mu=0, sigma=1):
    return (1/(sigma*np.sqrt(2*np.pi)))*np.exp(-0.5*((x-mu)/sigma)**2)

# Overtone amplitude KDEs
# -----------------------

# Choices for the Gaussian prior means
t0_choices = [-2, 0, 2, 4, 6]

sigma = 1

A_kde_list = []
for i, t0_choice in enumerate(t0_choices):
    # Calculate weights for each sample by evaluating the t0 prior at each t0
    # posterior sample
    weights = gaussian(t0_Hanford_posterior, mu=t0_choice, sigma=sigma)
    A_kde = gaussian_kde(A1_posterior, bw_method=0.4, weights=weights)
    A_kde_list.append(A_kde)

kde_grids_path = os.path.join(paths.data, "fm_t0_kde_grid.h5")

A_min, A_max = 0, 100
# evaluate KDEs over a grid and cache
x = np.linspace(A_min, A_max, 5000)
def get_kde_grid(i):
    kde = A_kde_list[i]
    return pd.Series(data=kde(x), index=x)
with mp.Pool(4) as pool:
    kde_grids = pool.map(get_kde_grid, range(len(A_kde_list)))
kde_grids_df = pd.DataFrame(dict(zip(t0_choices, kde_grids)))
kde_grids_df.to_hdf(kde_grids_path, "data")
