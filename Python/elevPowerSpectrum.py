#!/usr/bin/env python
# encoding: utf-8
"""
elevPowerSpectrum.py

Created by Geoff Cureton on 2011-03-06.
Copyright (c) 2011-2013 Geoff Cureton. All rights reserved.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
                                                          
"""

file_Date = '$Date$'
file_Revision = '$Revision$'
file_Author = '$Author$'
file_HeadURL = '$HeadURL$'
file_Id = '$Id$'

__author__ = 'G.P. Cureton <geoff.cureton@physics.org>'
__version__ = '$Id$'
__docformat__ = 'Epytext'

import os, sys, logging, traceback
from os import path,uname,environ

import string, copy
import numpy as np

# every module should have a LOG object
LOG = logging.getLogger(__file__)

######################################################################
###   Function to calculate a Phillips elevation power spectrum,   ###
###   and the desired phase relationships                          ###
######################################################################

def phillips_elev_spectrum(Scale,Power,NLCoupling,powExponent):
	
    LOG.info("Calculating the Phillips power spectrum...")

    N = Scale.N
    dk = Scale.delta_k
    k_max = Scale.k_max
    k = Scale.k

    LOG.info("N     = {:12d} ".format(N))
    LOG.info("k_max = {:12.6f} meters^{{-1}}".format(k_max))
    LOG.info("dk    = {:12.6f} meters^{{-1}}".format(dk))
    LOG.info("k_N   = {:12.6f} meters^{{-1}}".format(np.double(N/2)*dk))

    B      =  0.005		# Dimensionless constant
    g      =  9.81      # Acceleration due to gravity in m * s^{-2}
    v      =  5.        # Wind velocity in m * s^{-1}
    gammaT =  0.073     # Water surface tension in N * m^{-1}
    rho    =  1000.     # Water density in kg * m^{3}

    kc0 = g/(v*v)										# Lower wavenumber cutoff for wind waves
    k_gamma = np.minimum(np.sqrt(rho*g/gammaT),k_max)	# Upper wavenumber cutoff for 
                                                        #the capillary wave region

    nlCutoff   = 1.5	### Upper limit of source waves is nlCutoff*kc0
    elevStDev  = np.sqrt(0.25*B*(kc0**(-2.)-k_gamma**(-2.)))
    elevVar    = 0.25*B*(kc0**(-2.)-k_gamma**(-2.))
    slopeStDev = np.sqrt(0.5*B*np.log(k_gamma/kc0))
    slopeVar   = 0.5*B*np.log(k_gamma/kc0)

    LOG.info("Theoretical elevation stdev:    {} meters".format(elevStDev))
    LOG.info("Theoretical elevation variance: {} meters^{{2}}".format(elevVar))
    LOG.info("Theoretical slope stdev:        {}".format(slopeStDev))
    LOG.info("Theoretical slope variance:     {}".format(slopeVar))

    LOG.info("kc0 cutoff is:     {:12.6f} meters^{{-1}}".format(kc0))
    LOG.info("nlCutoff is:       {:12.6f}".format(nlCutoff))
    LOG.info("(nlCutoff*kc0) is: {:12.6f} meters^{{-1}}".format(nlCutoff*kc0))
    LOG.info("k_gamma cutoff is: {:12.6f} meters^{{-1}}".format(k_gamma))

    ### Determine the indices of components in the wavenumber range for the free waves
    sourceIndex = np.squeeze(np.where((k > kc0)*(k < nlCutoff*kc0)))

    LOG.info("Indices between kc0 and nlCutoff*kc0 (free1 and free2, the primary power indicies)...")
    LOG.info("sourceindex {} at wavenumber {}".format(sourceIndex,k[sourceIndex]))

    ### Define the structure containing the phase relationships between the free
    ### and bound waves
    NLCoupling.Nbound = np.shape(sourceIndex)[0]    # Number of coupling relationships
    NLCoupling.bound = NLCoupling.free1 = NLCoupling.free2 \
        = np.zeros(NLCoupling.Nbound)

    LOG.info("There are {:3d} bound-free pairs".format((NLCoupling.Nbound)))

    ### Determine the indices of the bound waves
    coupleIndex = 2*sourceIndex

    LOG.info("Indices between 2*kc0 and 2*nlCutoff*kc0 (bound, the coupled power indicies)...")
    LOG.info("coupleIndex {} at wavenumber {}".format(coupleIndex,k[coupleIndex]))

    ### Assign the free and bound wave indicies to the NLCoupling structure
    NLCoupling.free1 = sourceIndex
    NLCoupling.free2 = sourceIndex
    NLCoupling.bound = coupleIndex

    powExponent = np.double(powExponent)   # power law spectra are proportional 
                                           # to k^(-powExponent)

    nlProp = 0.65	# proportion of the total wave power due to bound waves

    ### Compute the total power S(k)
    totalPower = np.zeros(N)
    totalPower[0] = 0.
    totalPower[1:] = (k[1:] > kc0)*(k[1:] < k_gamma)*B/(2.*(k[1:]**powExponent))
    Power.totalPower = totalPower

    ### Set the bound power at the bound wavenumbers to 35% of the total power
    ### at those wavenumbers
    nlPower = np.zeros(N)
    nlPower[0] = 0.
    nlPower[coupleIndex] = nlProp*B/(2.*(k[coupleIndex]**powExponent))
    Power.nlPower = nlPower

    ### Define the primary power spectrum
    primaryPower = np.zeros(N)
    primaryPower = totalPower - nlPower
    Power.primaryPower = primaryPower

def phillips_elev_slope_variances_theory(k):
	LOG.info("Calculating the Elevation and Slope variances...")

	B      =  0.005		# Dimensionless constant
	g      =  9.81      # Acceleration due to gravity in m * s^{-2}
	v      =  5.        # Wind velocity in m * s^{-1}
	gammaT =  0.073     # Water surface tension in N * m^{-1}
	rho    =  1000.     # Water density in kg * m^{3}

	k_max = k[-1]
	kc0 = g/(v*v)									# Lower wavenumber cutoff for wind waves
	k_gamma = np.min(np.sqrt(rho*g/gammaT),k_max)	# Upper wavenumber cutoff for the 
                                                    # capillary wave region

	nlCutoff   = 1.5	### Upper limit of source waves is nlCutoff*kc0
	elevStDev  = np.sqrt(0.25*B*(kc0**(-2.)-k_gamma**(-2.)))
	elevVar    = 0.25*B*(kc0**(-2.)-k_gamma**(-2.))
	slopeStDev = np.sqrt(0.5*B*np.log(k_gamma/kc0))
	slopeVar   = 0.5*B*np.log(k_gamma/kc0)

	return elevVar,slopeVar

