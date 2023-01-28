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


# ----------------------------------------------------------------------------
# FUNDAMENTAL CONSTANTS

GMSUN_SI = 1.32712442099000e+20
C_SI = 2.99792458000000e+08
T_MSUN = GMSUN_SI / C_SI**3

gwq150914_remnant_mass = 69
tM = gwq150914_remnant_mass*T_MSUN

# ----------------------------------------------------------------------------
# CODE CONVERSIONS
#
# Different codes define the ringdown tone amplitudes differently. We will
# rescale them all to the used in the original PRL by Isi et al. That is as
# follows
# 
# $$
# h_+ = \sum A_n \frac{1}{2} \left(1 + \cos^2 \iota\right)
# \cos (\omega_n t + \phi_n)\,  \exp(-t/\tau_n),
# $$
# 
# $$
# h_\times = \sum A_n \cos\iota\,\sin (\omega_n t + \phi_n)\,  \exp(-t/\tau_n) ,
# $$
# for each $\ell = |m| = 2$ tone, labeled by the overtone number $n$ such that
# $\tau_n < \tau_{n+1}$.  The outputs of the `ringdown`, `PyRing` and Finch &
# Moore codes differ from the above definition by the following factors:

# ringdown scaling
A_scale = 2

# pyring scaling
A_scale_pr = 0.63*1E-21

# Finch & Moore scaling
A_scale_fm = 1E-21

