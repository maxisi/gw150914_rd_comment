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

for i in tqdm(range(len(A1_posterior))):

    # Get corresponding params
    A1 = A1_posterior[i]
    t0 = t0_Hanford_posterior[i]
    Mf = Mf_posterior[i]
    chif = chif_posterior[i]

    # Calculate damping time
    omega, A, C = qnm_func(chif)
    tau = -1/np.imag(omega)

    # Calculate A1 at t_ref
    A1_posterior_tref.append(A1*np.exp(t0/tau))

A1_posterior_tref = np.array(A1_posterior_tref)
np.savetxt(A1_ref_cache, A1_posterior_tref)
