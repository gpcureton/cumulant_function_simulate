#!/usr/bin/env python
# encoding: utf-8
"""
cumulantFunctionSimulate.py


Purpose: This program simulates multiple elevation and slope realisations with
    given spectra, introducing quadratic phase coupling at selected wavenumbers.
    Glint realisations are computed from the slope realisations. The average 
    elevation, slope and glint moments and cumulants, and moment and cumulant 
    functions are computed. The results are written to a HDF5 file.
    
Input:
    N            : Data length                                
    NN           : Bispectrum data length                     
    delta_x      : Spatial increment in meters                
    N_r          : Number of realisations                     
    spectrumType : Form of the elevation power spectrum       
    specExp      : Elevation power spectrum is proportional to
                   k^{-specExp}                               
    nlSwitch     : Elevation phase coupling on/off            

Output:
    Output of results are written to HDF5 file.               

Details:
    * 

Preconditions:
    * 

Optional:
    * 

Minimum commandline:

    python progname.py  --input_files=INPUTFILES

where...

    INPUTFILES: The fully qualified path to the input files. May be a directory or a file glob.


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

import sys, time, string, getopt, copy
import numpy as np
from numpy import pi,sin,cos,tan,sqrt,abs,exp
from numpy.fft import fft,ifft
from numpy import float64 as double

from scipy import stats as stats
import time

#import matplotlib
#import matplotlib.cm as cm
#from matplotlib.colors import ListedColormap
#from matplotlib.figure import Figure

#matplotlib.use('Agg')
#from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas

#matplotlib.use('WXAgg')
#from matplotlib.backends.backend_wxagg import FigureCanvasAgg as FigureCanvas

# This must come *after* the backend is specified.
#import matplotlib.pyplot as ppl

import optparse as optparse

import tables as pytables
from tables import exceptions as pyEx

from elevPowerSpectrum import phillips_elev_spectrum

# Initialise the scale structure
class ScaleStruct :
    def __init__(self, N, delta_x) :
        self.N = N
        self.delta_x = delta_x
        self.x_max = double(N-1)*delta_x
        self.x = np.arange(N)*delta_x

        self.k_max = 2.*pi/delta_x
        self.delta_k = self.k_max/double(N-1)
        self.k_N = double(N/2)*self.delta_k # Nyquist "wavenumber"
        self.k = np.arange(N)*self.delta_k

class GeomStruct :
    def __init__(self, N_angles, angleLo, angleHi) :
        self.N_angles = N_angles
        self.source_angle = np.arange(N_angles)
        self.detector_angle = np.arange(N_angles)
        self.xi_min = np.arange(N_angles)
        self.xi_0 = np.arange(N_angles)
        self.xi_max = np.arange(N_angles)

        self.angleRange=(angleHi-angleLo)

        if ((N_angles-1) == 0) :
            self.d_angle = 0
        else :
            self.d_angle = self.angleRange/(N_angles-1)

        self.start_angle = angleLo
        d2r = pi/180.
        r2d = 180./pi
        beta = 0.68*d2r

        self.source_angle = ((self.start_angle + np.arange(N_angles,dtype=double)*self.d_angle))*d2r
        self.detector_angle = 0.0*d2r
        gamma = (self.source_angle - self.detector_angle)/2.
        self.xi_0 = tan(gamma)
        self.xi_min = self.xi_0 - (1.0 + self.xi_0**2.)*(beta/4.)
        self.xi_max = self.xi_0 + (1.0 + self.xi_0**2.)*(beta/4.)
        self.dxi = self.xi_max - self.xi_min

class PowerStruct :
    def __init__(self,N,spectrumType) :
        self.spectrumType = spectrumType
        self.primaryPower = np.zeros(N,double)
        self.nlPower = np.zeros(N,double)
        self.totalPower = np.zeros(N,double)

class NLCouplingStruct :
    def __init__(self,N) :
        self.Nbound = 0
        self.bound = np.zeros(N, np.long)
        self.free1 = np.zeros(N, np.long)
        self.free2 = np.zeros(N, np.long)

class DataStatsStruct :
    def __init__(self,numMoments) :
        self.mean = 0.
        self.variance = 0.
        self.skewness = 0.
        self.moments = np.zeros(numMoments, double)
        self.cumulants = np.zeros(numMoments, double)

    def cumulantsFromMoments(self):
        self.cumulants[0] = self.moments[0]
        self.cumulants[1] = self.moments[1] - self.moments[0]**2.
        self.cumulants[2] = self.moments[2] - 3.*self.moments[0]*self.moments[1] + 2.*self.moments[0]**3.
        self.mean     = self.cumulants[0]
        self.variance = self.cumulants[1]
        self.skewness = self.cumulants[2]/((self.cumulants[1])**1.5)


def cumulantFunctionSimulate(N,NN,delta_x,N_r,spectrumType,specExp,nlSwitch):

    """
    Determine the various scale parameters and populate the 
    Scale class 
    """

    Scale = ScaleStruct(N, delta_x)

    print 'Scale.N       = %15d' % (Scale.N)
    print 'Scale.delta_x = %15.6f meters' % (Scale.delta_x)
    print 'Scale.x_max   = %15.6f meters' % (Scale.x_max)
    print 'Scale.k_max   = %15.6f meters^{-1}' % (Scale.k_max)
    print 'Scale.delta_k = %15.6f meters^{-1}' % (Scale.delta_k)
    print 'Scale.k_N     = %15.6f meters^{-1}' % (Scale.k_N)

    # Make local copies of Scale attributes
    x_max   = Scale.x_max
    x       = Scale.x
    k_max   = Scale.k_max
    delta_k = Scale.delta_k
    k_N     = Scale.k_N
    k       = Scale.k

    """
       Populate the GEOM structure with angular quantities
    """

    N_angles = 5L
    angleLo=10.
    angleHi=30.
    Geom = GeomStruct(N_angles,angleLo,angleHi)

    """
        Populate the elevation power spectrum structures ElevPower, 
        and NLCoupling                                          
    """

    ElevPower = PowerStruct(N,spectrumType)

    NLCoupling = NLCouplingStruct(N)

    phillips_elev_spectrum(Scale,ElevPower,NLCoupling,specExp)

    print '\nFirst component indicies for free waves ...',NLCoupling.free1,Scale.k[NLCoupling.free1]/Scale.k_N

    print 'Second component indicies for free waves ...',NLCoupling.free2,Scale.k[NLCoupling.free2]/Scale.k_N

    print 'Indicies for bound waves...',NLCoupling.bound,Scale.k[NLCoupling.bound]/Scale.k_N

    totalElevPower = ElevPower.totalPower
    primaryElevPower = ElevPower.primaryPower
    nlElevPower = ElevPower.nlPower

    print "\nElevation stdev from power vector: %10.6f meters " % \
        (sqrt(np.sum(totalElevPower)*delta_k))
    print "Elevation variance from power vector: %10.6f meters^{2} " % \
        (np.sum(totalElevPower)*delta_k)

    print "\nTotal elevation power at the bound wavenumbers...",totalElevPower[NLCoupling.bound]
    print "Free elevation power at the bound wavenumbers...",ElevPower.primaryPower[NLCoupling.bound]
    print "Bound elevation power at the bound wavenumbers...",ElevPower.nlPower[NLCoupling.bound]
    print "Ratio of bound to free elevation power at the bound wavenumbers...",\
            ElevPower.nlPower[NLCoupling.bound]/totalElevPower[NLCoupling.bound]

    """
        Initialise the slope power spectrum structure
    """
    SlopePower = copy.deepcopy(ElevPower)
    SlopePower.primaryPower = k*k*ElevPower.primaryPower
    SlopePower.nlPower = k*k*ElevPower.nlPower
    SlopePower.totalPower = k*k*ElevPower.totalPower

    totalSlopePower = SlopePower.totalPower
    primarySlopePower = SlopePower.primaryPower
    nlSlopePower = SlopePower.nlPower

    print "\nSlope stdev from power vector: %10.6f meters " % \
        (sqrt(np.sum(totalSlopePower)*delta_k))
    print "Slope variance from power vector: %10.6f meters^{2} " % \
        (np.sum(totalSlopePower)*delta_k)

    print "\nTotal slope power at the bound wavenumbers...",totalSlopePower[NLCoupling.bound]
    print "Free slope power at the bound wavenumbers...",SlopePower.primaryPower[NLCoupling.bound]
    print "Bound slope power at the bound wavenumbers...",SlopePower.nlPower[NLCoupling.bound]
    print "Ratio of bound to free slope power at the bound wavenumbers...",\
            SlopePower.nlPower[NLCoupling.bound]/totalSlopePower[NLCoupling.bound]

    """
        Initialise the curvature power spectrum structure
    """
    CurvaturePower = copy.deepcopy(ElevPower)
    CurvaturePower.primaryPower = k*k*k*k*ElevPower.primaryPower
    CurvaturePower.nlPower = k*k*k*k*ElevPower.nlPower
    CurvaturePower.totalPower = k*k*k*k*ElevPower.totalPower

    totalCurvaturePower = CurvaturePower.totalPower
    primaryCurvaturePower = CurvaturePower.primaryPower
    nlCurvaturePower = CurvaturePower.nlPower

    print "\nCurvature stdev from power vector: %10.6f meters^{-1}" % \
        (sqrt(np.sum(totalCurvaturePower)*delta_k))
    print "Curvature variance from power vector: %10.6f meters^{-2}" % \
        (np.sum(totalCurvaturePower)*delta_k)

    print "\nTotal curvature power at the bound wavenumbers...",totalCurvaturePower[NLCoupling.bound]
    print "Free curvature power at the bound wavenumbers...",CurvaturePower.primaryPower[NLCoupling.bound]
    print "Bound curvature power at the bound wavenumbers...",CurvaturePower.nlPower[NLCoupling.bound]
    print "Ratio of bound to free curvature power at the bound wavenumbers...",\
            CurvaturePower.nlPower[NLCoupling.bound]/totalCurvaturePower[NLCoupling.bound]

    """
        Compute the total elevation amplitude, phase and spectrum,
        and the second moment function 
    """

    totalElevAmplitude = np.zeros(N,dtype=double)
    totalElevAmplitude = sqrt(0.5*totalElevPower*delta_k)
    totalElevAmplitude[N/2+1 :] = totalElevAmplitude[1L : N/2][::-1]
    totalElevSpectrum = np.zeros(N,dtype=np.complex64)

    primaryElevAmplitude = np.zeros(N,dtype=double)
    primaryElevAmplitude = sqrt(0.5*primaryElevPower*delta_k)
    primaryElevAmplitude[N/2+1 :] = primaryElevAmplitude[1 : N/2][::-1]

    nlElevAmplitude = np.zeros(N,dtype=double)
    nlElevAmplitude = sqrt(0.5*nlElevPower*delta_k)
    nlElevAmplitude[N/2+1 :] = nlElevAmplitude[1 : N/2][::-1]

    print "\nElevation stdev from amplitude vector: %10.6f meters " % \
            (sqrt(np.sum(totalElevAmplitude**2.)))
    print "Elevation variance from amplitude vector: %10.6f meters^{2}" % \
            (np.sum(totalElevAmplitude**2.))

    testElevPhase = np.random.rand(N)*2.*pi - pi
    totalElevSpectrum = totalElevAmplitude*(cos(testElevPhase) + 1j*sin(testElevPhase))
    totalElevSpectrum[N/2+1 :] = np.conjugate(totalElevSpectrum[1 : N/2][::-1])
    totalElevSurface = fft(totalElevSpectrum)

    print "\nElevation mean from surface:    %10.6f meters " % \
            np.mean(totalElevSurface.real)
    print "Elevation stdev from surface:    %10.6f meters " % \
            np.std(totalElevSurface.real)
    print "Elevation variance from surface: %10.6f meters^{2} " % \
            np.var(totalElevSurface.real)

    totalElevAvgPower = np.zeros(N,dtype=double)
    primaryElevAvgPower = np.zeros(N,dtype=double)
    nlElevAvgPower = np.zeros(N,dtype=double)

    """
        Compute the total slope amplitude, phase and spectrum
    """

    totalSlopeAmplitude = np.zeros(N,dtype=double)
    totalSlopeAmplitude = sqrt(0.5*totalSlopePower*delta_k)
    totalSlopeAmplitude[N/2+1 :] = totalSlopeAmplitude[1L : N/2][::-1]
    totalSlopeSpectrum = np.zeros(N,dtype=np.complex64)
    totalSlopeSurface = np.zeros(N,dtype=np.complex64)

    primarySlopeAmplitude = np.zeros(N,dtype=double)
    primarySlopeAmplitude = sqrt(0.5*primarySlopePower*delta_k)
    primarySlopeAmplitude[N/2+1 :] = primarySlopeAmplitude[1L : N/2][::-1]
    primarySlopeSpectrum = np.zeros(N,dtype=np.complex64)
    primarySlopeSurface = np.zeros(N,dtype=np.complex64)

    nlSlopeAmplitude = np.zeros(N,dtype=double)
    nlSlopeAmplitude = sqrt(0.5*nlSlopePower*delta_k)
    nlSlopeAmplitude[N/2+1 :] = nlSlopeAmplitude[1L : N/2][::-1]
    nlSlopeSpectrum = np.zeros(N,dtype=np.complex64)
    nlSlopeSurface = np.zeros(N,dtype=np.complex64)

    print "\nSlope stdev from amplitude vector: %10.6f" % \
            (sqrt(np.sum(totalSlopeAmplitude**2.)))
    print "Slope variance from amplitude vector: %10.6f" % \
            (np.sum(totalSlopeAmplitude**2.))

    totalSlopeSpectrum = totalSlopeAmplitude*(+sin(testElevPhase) - 1j*cos(testElevPhase))
    totalSlopeSpectrum[N/2+1 :] = np.conjugate(totalSlopeSpectrum[1L : N/2][::-1])
    totalSlopeSurface = fft(totalSlopeSpectrum)

    print "\nSlope mean from surface:    %10.6f" % \
            np.mean(totalSlopeSurface.real)
    print "Slope stdev from surface:    %10.6f" % \
            np.std(totalSlopeSurface.real)
    print "Slope variance from surface: %10.6f" % \
            np.var(totalSlopeSurface.real)

    totalSlopeAvgPower = np.zeros(N,dtype=double)
    primarySlopeAvgPower = np.zeros(N,dtype=double)
    nlSlopeAvgPower = np.zeros(N,dtype=double)

    """
        Compute the total curvature amplitude, phase and spectrum
    """

    totalCurvatureAmplitude = np.zeros(N,dtype=double)
    totalCurvatureAmplitude = sqrt(0.5*totalCurvaturePower*delta_k)
    totalCurvatureAmplitude[N/2+1 :] = totalCurvatureAmplitude[1L : N/2][::-1]
    totalCurvatureSpectrum = np.zeros(N,dtype=np.complex64)
    totalCurvatureSurface = np.zeros(N,dtype=np.complex64)

    primaryCurvatureAmplitude = np.zeros(N,dtype=double)
    primaryCurvatureAmplitude = sqrt(0.5*primaryCurvaturePower*delta_k)
    primaryCurvatureAmplitude[N/2+1 :] = primaryCurvatureAmplitude[1L : N/2][::-1]
    primaryCurvatureSpectrum = np.zeros(N,dtype=np.complex64)
    primaryCurvatureSurface = np.zeros(N,dtype=np.complex64)

    nlCurvatureAmplitude = np.zeros(N,dtype=double)
    nlCurvatureAmplitude = sqrt(0.5*nlCurvaturePower*delta_k)
    nlCurvatureAmplitude[N/2+1 :] = nlCurvatureAmplitude[1L : N/2][::-1]
    nlCurvatureSpectrum = np.zeros(N,dtype=np.complex64)
    nlCurvatureSurface = np.zeros(N,dtype=np.complex64)

    print "\nCurvature stdev from amplitude vector: %10.6f meters^{-1} " % \
            (sqrt(np.sum(totalCurvatureAmplitude**2.)))
    print "Curvature variance from amplitude vector: %10.6f meters^{-2}" % \
            (np.sum(totalCurvatureAmplitude**2.))

    totalCurvatureSpectrum = totalCurvatureAmplitude*(-cos(testElevPhase) - 1j*sin(testElevPhase))
    totalCurvatureSpectrum[N/2+1 :] = np.conjugate(totalCurvatureSpectrum[1L : N/2][::-1])
    totalCurvatureSurface = fft(totalCurvatureSpectrum)

    print "\nCurvature mean from surface:    %10.6f meters^{-1} " % \
            np.mean(totalCurvatureSurface.real)
    print "Curvature stdev from surface:    %10.6f meters^{-1}" % \
            np.std(totalCurvatureSurface.real)
    print "Curvature variance from surface: %10.6f meters^{-2} " % \
            np.var(totalCurvatureSurface.real)

    totalCurvatureAvgPower = np.zeros(N,dtype=double)
    primaryCurvatureAvgPower = np.zeros(N,dtype=double)
    nlCurvatureAvgPower = np.zeros(N,dtype=double)

    """
	    Define the glint, glint spectrum and glint power
    """

    glint = np.zeros(N,dtype=double)
    glintSpectrum = np.zeros(N,dtype=np.complex64)
    totalGlintAvgPower = np.zeros((Geom.N_angles,N),dtype=double) # DBLARR(N,GEOM.N_angles)

    """
	    Define the various point estimators for the elevation,
	    slope and glint
    """

    numMoments = 3

    elevStats = DataStatsStruct(numMoments)
    slopeStats = DataStatsStruct(numMoments)
    curvatureStats = DataStatsStruct(numMoments)
    glintStats = [DataStatsStruct(numMoments) for geoms in np.arange(Geom.N_angles) ]

    """
        Loop through the surface realisations for the quadratically
        coupled oscillations
    """

    seed = 30
    N_r_cum = 0
    angleRuns = np.zeros(Geom.N_angles,np.long)
    angleRunsCum = np.zeros(Geom.N_angles,np.long)

    #time.sleep(3.)
    t1 = time.time()

    while (angleRuns.sum() < N_r*Geom.N_angles) :

        N_r_cum += 1

        t2 = time.time()
        #print "Elapsed time = ",t2-t1
        if ((t2-t1) > 0.5):
            print "\n>>>>>>>>>>>>>>>>>>>>>\n"
            print "Computing realisation: %d at time %f" % (N_r_cum,(t2-t1))
            t1 = time.time()


		### Compute the independent phases for this realisation
        primaryElevPhase = np.random.rand(N)*2.*pi - pi # RANDOMU(seed,N)*2.D*!DPI - !DPI
        nlElevPhase      = np.random.rand(N)*2.*pi - pi # RANDOMU(seed,N)*2.D*!DPI - !DPI

        ### Apply the phase correlations between the free and bound wavenumbers for the nonlinear
        ### component, if (nlSwitch==1)
        if (nlSwitch==1) :
            nlElevPhase[NLCoupling.bound] = primaryElevPhase[NLCoupling.free1] + \
                    primaryElevPhase[NLCoupling.free2]

        """
		    Compute the elevation realisations from the elevation spectra
		    and the synthesised phases
        """
		
		### Calculate the elevation spectrum for the free waves
        primaryElevSpectrum = primaryElevAmplitude*(cos(primaryElevPhase) + 1j*sin(primaryElevPhase))
        primaryElevSpectrum[N/2+1 :] = np.conjugate(primaryElevSpectrum[1 : N/2][::-1])

        ### Calculate the elevation spectrum for the bound waves
        nlElevSpectrum = nlElevAmplitude*(cos(nlElevPhase) + 1j*sin(nlElevPhase))
        nlElevSpectrum[N/2+1 :] = np.conjugate(nlElevSpectrum[1 : N/2][::-1])

		### Compute specific realisation of the free and bound waves. Nonlinear elevation
		### (totalElevSurface) is sum of free and bound waves.
        primaryElevSurface = fft(primaryElevSpectrum)                      ### Free waves
        nlElevSurface = fft(nlElevSpectrum)                                ### Bound waves
        totalElevSurface = primaryElevSurface + nlElevSurface              ### Total surface
 
		### Compute the average power spectrum for free, bound and total elevation waves
        primaryElevAvgPower += abs(ifft(primaryElevSurface))**2.
        nlElevAvgPower      += abs(ifft(nlElevSurface))**2.
        totalElevAvgPower   += abs(ifft(totalElevSurface))**2.

        #print "\tElevation stdev from power vector:    %10.6e meters" % \
            #(sqrt(np.sum(abs(ifft(totalElevSurface.real))**2.)))
        #print "\tElevation variance from power vector: %10.6e meters^{2}\n" % \
            #(np.sum(abs(ifft(totalElevSurface.real))**2.))

        ### Compute the elevation moments

        elevStats.mean     += np.mean(totalElevSurface.real)
        elevStats.variance += np.var(totalElevSurface.real)
        elevStats.skewness += stats.skew(totalElevSurface.real)

        elevStats.moments += [ np.sum(totalElevSurface.real    )/double(N), \
                               np.sum(totalElevSurface.real**2.)/double(N), \
                               np.sum(totalElevSurface.real**3.)/double(N)]

		### Compute the Fourier spectrum of the total surface
        totalElevSpectrum = ifft(totalElevSurface)

        """
		    Compute the slope realisations from the slope spectra
		    and the synthesised phases
        """
		
        ### Calculate the slope spectrum for the free waves
        primarySlopeSpectrum = primarySlopeAmplitude*(sin(primaryElevPhase) - 1j*cos(primaryElevPhase))
        #primarySlopeSpectrum += 0.00001*MAX(totalSlopeAmplitude)*RANDOMN(seed,N)
        primarySlopeSpectrum[N/2+1 :] = np.conjugate(primarySlopeSpectrum[1 : N/2][::-1])

        ### Calculate the slope spectrum for the bound waves
        nlSlopeSpectrum = nlSlopeAmplitude*(sin(nlElevPhase) - 1j*cos(nlElevPhase))
        nlSlopeSpectrum[N/2+1 :] = np.conjugate(nlSlopeSpectrum[1 : N/2][::-1])

        ### Compute specific realisation of the free and bound waves. Nonlinear slope
        ### (totalSlopeSurface) is sum of free and bound waves.
        primarySlopeSurface = fft(primarySlopeSpectrum)                    ### Free waves
        nlSlopeSurface = fft(nlSlopeSpectrum)                              ### Bound waves
        totalSlopeSurface = primarySlopeSurface + nlSlopeSurface           ### Total surface

        ### Compute the average power spectrum for free, bound and total elevation waves
        primarySlopeAvgPower += abs(ifft(primarySlopeSurface))**2.
        nlSlopeAvgPower += abs(ifft(nlSlopeSurface))**2.
        totalSlopeAvgPower += abs(ifft(totalSlopeSurface))**2.

        #print "\n\tSlope stdev from power vector:    %10.6e meters" % \
            #(sqrt(np.sum(abs(ifft(totalSlopeSurface.real))**2.)))
        #print "\tSlope variance from power vector: %10.6e meters^{2}\n" % \
            #(np.sum(abs(ifft(totalSlopeSurface.real))**2.))

        ### Compute the slope moments

        slopeStats.mean     += np.mean(totalSlopeSurface.real)
        slopeStats.variance += np.var(totalSlopeSurface.real)
        slopeStats.skewness += stats.skew(totalSlopeSurface.real)

        slopeStats.moments += [ np.sum(totalSlopeSurface    )/double(N), \
                                np.sum(totalSlopeSurface**2.)/double(N), \
                                np.sum(totalSlopeSurface**3.)/double(N) ]

		### Compute the Fourier spectrum of the total surface
        totalSlopeSpectrum = ifft(totalSlopeSurface)

        """
		    Compute the curvature realisations from the curvature spectra
		    and the synthesised phases
        """
		
        ### Calculate the curvature spectrum for the free waves
        primaryCurvatureSpectrum = primaryCurvatureAmplitude*(-cos(primaryElevPhase) - 1j*sin(primaryElevPhase))
        #primaryCurvatureSpectrum += 0.00001*MAX(totalCurvatureAmplitude)*RANDOMN(seed,N)
        primaryCurvatureSpectrum[N/2+1 :] = np.conjugate(primaryCurvatureSpectrum[1 : N/2][::-1])

        ### Calculate the curvature spectrum for the bound waves
        nlCurvatureSpectrum = nlCurvatureAmplitude*(-cos(nlElevPhase) - 1j*sin(nlElevPhase))
        nlCurvatureSpectrum[N/2+1 :] = np.conjugate(nlCurvatureSpectrum[1 :N/2][::-1])

        ### Compute specific realisation of the free and bound waves. Nonlinear curvature
        ### (totalCurvatureSurface) is sum of free and bound waves.
        primaryCurvatureSurface = fft(primaryCurvatureSpectrum).real                    ### Free waves
        #primaryCurvatureSurface *= 1./((1. + primarySlopeSurface**2.)**1.5)

        nlCurvatureSurface = fft(nlCurvatureSpectrum).real                              ### Bound waves
        #nlCurvatureSurface *= 1./((1. + nlSlopeSurface**2.)**1.5)

        totalCurvatureSurface = primaryCurvatureSurface + nlCurvatureSurface       ### Total surface

        #totalCurvatureSurface *= 1./(    (1. + totalSlopeSurface**2.)**1.5)
        #totalCurvatureSurface *= 1./(sqrt(1. + totalSlopeSurface**2.)**3.)       ### Total surface

        ### Compute the average power spectrum for free, bound and total elevation waves
        primaryCurvatureAvgPower += abs(ifft(primaryCurvatureSurface))**2.
        nlCurvatureAvgPower += abs(ifft(nlCurvatureSurface))**2.
        totalCurvatureAvgPower += abs(ifft(totalCurvatureSurface))**2.

        #print "\n\tCurvature stdev from power vector:    %10.6e meters" % \
            #(sqrt(np.sum(abs(ifft(totalCurvatureSurface.real))**2.)))
        #print "\tCurvature variance from power vector: %10.6e meters^{2}\n" % \
            #(np.sum(abs(ifft(totalCurvatureSurface.real))**2.))

        ### Compute the curvature moments

        curvatureStats.mean     += np.mean(totalCurvatureSurface.real)
        curvatureStats.variance += np.var(totalCurvatureSurface.real)
        curvatureStats.skewness += stats.skew(totalCurvatureSurface.real)

        curvatureStats.moments += [ np.sum(totalCurvatureSurface    )/double(N), \
                                    np.sum(totalCurvatureSurface**2.)/double(N), \
                                    np.sum(totalCurvatureSurface**3.)/double(N) ]

		### Compute the Fourier spectrum of the total surface
        totalCurvatureSpectrum = ifft(totalCurvatureSurface)

        """
		    Loop through the geometries in the GEOM structure
        """

        for angle in np.arange(Geom.N_angles) :

            """
            Check if we have finished processing for this
            angle.
            """

            if (angleRuns[angle] < N_r) :

                #print "\n\tProcessing angle ",angle," for run ",angleRuns[angle]+1, \
                #" --> attempt ",angleRunsCum[angle]+1

                """
                    Compute the glint realisation from the slope
                    realisations
                """

                slopeMin = Geom.xi_min[angle]
                slopeMax = Geom.xi_max[angle]

                glint = np.double(totalSlopeSurface.real > slopeMin) * np.double(totalSlopeSurface.real < slopeMax)

                ### Check if all glint elements vanish
                result = np.where(glint)

                #if (glint.sum() == 0.) :
                if (np.shape(np.squeeze(np.where(glint))) == (0,)) :

                    #print "\tZero-glint realisation angle ",angle, \
                    #" for run ",angleRuns[angle]+1, \
                    #" --> attempt ",angleRunsCum[angle]+1

                    ### There are no glints, add to attempts count
                    ### If this angle fails and there are no glints, then steeper 
                    ### angles will fail also, so break out of the angle loop and 
                    ### proceed to the next realisation...

                    angleRunsCum[angle:Geom.N_angles] += 1
                    break

                else :

                    #print "\tSuccessful realisation angle ",angle, \
                    #" for run ",angleRuns[angle] + 1, \
                    #" --> attempt ",angleRunsCum[angle]+1

                    angleRuns[angle] += 1
                    angleRunsCum[angle] += 1

                    ### Compute the glint moments

                    glintStats[angle].mean     += np.mean(glint.real)
                    glintStats[angle].variance += np.var( glint.real)
                    glintStats[angle].skewness += stats.skew(glint.real)

                    glintStats[angle].moments += [ np.sum(glint.real    )/double(N), \
                                                   np.sum(glint.real**2.)/double(N), \
                                                   np.sum(glint.real**3.)/double(N) ]

                    ### Compute the Fourier spectrum of this glint realisation (using ifft 
                    ### which has same normalisation as IDL FFT routine).
                    glintSpectrum = ifft(glint)

                    """
                    Compute the average glint power spectrum
                    """

                    totalGlintAvgPower[angle] += abs(glintSpectrum)**2.

                ### End checking for zero-glint of this angle

            ### End checking for completion of this angle

        ### End angle loop

    ### End realisation loop
    
    print ""
    print "AngleRuns:    ",angleRuns," ... for total of ", \
        int(np.sum(angleRuns))," / ",N_r*Geom.N_angles
    print "AngleRunsCum: ",angleRunsCum
	
    N_runs = N_r
    N_r = N_r_cum
    print  "Final N_r_cum = ",N_r_cum

    """
        Have a look at the scipy calculated stats
    """
    elevStats.mean     /= double(N_runs)
    elevStats.variance /= double(N_runs)
    elevStats.skewness /= double(N_runs)

    slopeStats.mean     /= double(N_runs)
    slopeStats.variance /= double(N_runs)
    slopeStats.skewness /= double(N_runs)

    curvatureStats.mean     /= double(N_runs)
    curvatureStats.variance /= double(N_runs)
    curvatureStats.skewness /= double(N_runs)

    for geoms in np.arange(Geom.N_angles) :
        glintStats[geoms].mean     /= double(angleRunsCum[geoms])
        glintStats[geoms].variance /= double(angleRunsCum[geoms])
        glintStats[geoms].skewness /= double(angleRunsCum[geoms])

    print "\nPython Elevation mean %10.6e" % elevStats.mean
    print "Python Elevation stdev %10.6e" % sqrt(elevStats.variance)
    print "Python Elevation variance %10.6e" % elevStats.variance
    print "Python Elevation skewness %10.6e" % elevStats.skewness

    print "\nPython Slope mean %10.6e" % slopeStats.mean
    print "Python Slope stdev %10.6e" % sqrt(slopeStats.variance)
    print "Python Slope variance %10.6e" % slopeStats.variance
    print "Python Slope skewness %10.6e" % slopeStats.skewness

    print "\nPython Curvature mean %10.6e" % curvatureStats.mean
    print "Python Curvature stdev %10.6e" % sqrt(curvatureStats.variance)
    print "Python Curvature variance %10.6e" % curvatureStats.variance
    print "Python Curvature skewness %10.6e" % curvatureStats.skewness

    print "\nPython Glint mean, variance and skewness ...\n"
    for geoms in np.arange(Geom.N_angles) :
        print "\tAngle %1d:\t\t%10.6e\t%10.6e\t%10.6e" % (geoms,\
            glintStats[geoms].mean,\
            glintStats[geoms].variance,\
            glintStats[geoms].skewness)

    """
        Compute the moments and cumulants our way
    """
    elevStats.moments /= double(N_runs)
    elevStats.cumulantsFromMoments()
    slopeStats.moments /= double(N_runs)
    slopeStats.cumulantsFromMoments()
    curvatureStats.moments /= double(N_runs)
    curvatureStats.cumulantsFromMoments()
    for geoms in np.arange(Geom.N_angles) :
        glintStats[geoms].moments /= double(angleRunsCum[geoms])
        glintStats[geoms].cumulantsFromMoments()

    print "\nElevation first moment %10.6e" % elevStats.moments[0]
    print "Elevation second moment %10.6e" % elevStats.moments[1]
    print "Elevation third moment %10.6e" % elevStats.moments[2]

    print "\nElevation first cumulant %10.6e" % elevStats.cumulants[0]
    print "Elevation second cumulant %10.6e" % elevStats.cumulants[1]
    print "Elevation third cumulant %10.6e" % elevStats.cumulants[2]

    print "\nElevation mean %10.6e" % elevStats.mean
    print "Elevation stdev %10.6e" % sqrt(elevStats.variance)
    print "Elevation variance %10.6e" % elevStats.variance
    print "Elevation skewness %10.6e" % elevStats.skewness

    print "\nSlope first moment %10.6e" % slopeStats.moments[0]
    print "Slope second moment %10.6e" % slopeStats.moments[1]
    print "Slope third moment %10.6e" % slopeStats.moments[2]

    print "\nSlope first cumulant %10.6e" % slopeStats.cumulants[0]
    print "Slope second cumulant %10.6e" % slopeStats.cumulants[1]
    print "Slope third cumulant %10.6e" % slopeStats.cumulants[2]

    print "\nSlope mean %10.6e" % slopeStats.mean
    print "Slope stdev %10.6e" % sqrt(slopeStats.variance)
    print "Slope variance %10.6e" % slopeStats.variance
    print "Slope skewness %10.6e\n" % slopeStats.skewness

    print "\nCurvature first moment %10.6e" % curvatureStats.moments[0]
    print "Curvature second moment %10.6e" % curvatureStats.moments[1]
    print "Curvature third moment %10.6e" % curvatureStats.moments[2]

    print "\nCurvature first cumulant %10.6e" % curvatureStats.cumulants[0]
    print "Curvature second cumulant %10.6e" % curvatureStats.cumulants[1]
    print "Curvature third cumulant %10.6e" % curvatureStats.cumulants[2]

    print "\nCurvature mean %10.6e" % curvatureStats.mean
    print "Curvature stdev %10.6e" % sqrt(curvatureStats.variance)
    print "Curvature variance %10.6e" % curvatureStats.variance
    print "Curvature skewness %10.6e\n" % curvatureStats.skewness

    print "Glint moments ...\n"
    for geoms in np.arange(Geom.N_angles) :
        print "\tAngle %1d:\t\t%10.6e\t%10.6e\t%10.6e" % (geoms,\
            glintStats[geoms].moments[0],\
            glintStats[geoms].moments[1],\
            glintStats[geoms].moments[2])

    print "\nGlint cumulants ...\n"
    for geoms in np.arange(Geom.N_angles) :
        print "\tAngle %1d:\t\t%10.6e\t%10.6e\t%10.6e" % (geoms,\
            glintStats[geoms].cumulants[0],\
            glintStats[geoms].cumulants[1],\
            glintStats[geoms].cumulants[2])

    """
        compute the average elevation moment and cumulant functions.
    """
    elevSecondMomentFunction =  fft(totalElevAvgPower)
    elevSecondMomentFunction /= elevSecondMomentFunction.real[0]
	#elevSecondCumulantFunction = (elevMoments[1]*elevSecondMomentFunction - elevMoments[0]^2.D)/elevCumulants[1]

    """
        compute the average slope moment and cumulant functions.
    """
    slopeSecondMomentFunction =  fft(totalSlopeAvgPower)
    slopeSecondMomentFunction /= slopeSecondMomentFunction.real[0]

    """
        compute the average glint moment and cumulant functions.
    """
    glintSecondMomentFunction = np.zeros((Geom.N_angles,N),dtype=double)
    for geoms in np.arange(Geom.N_angles) :
        glintSecondMomentFunction[geoms] = fft(totalGlintAvgPower[geoms]).real
        glintSecondMomentFunction[geoms] /= glintSecondMomentFunction[geoms][0]

"""
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;;   Compute the elevation estimators   ;;;
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	;;; Compute the estimators of the IDL moments
	elevMean /= DOUBLE(N_r)
	elevVariance /= DOUBLE(N_r)
	elevSkewness /= DOUBLE(N_r)

	;;; Compute the average elevation moments
	elevMoments /= DOUBLE(N_r)

	;;; Compute the average elevation Cumulants
	elevCumulants  = DBLARR(numMoments)
	cumulantsFromMoments,elevMoments,elevCumulants
	
	;;; Compute the second order elev moment function and spectra
	totalElevAvgPower /= DOUBLE(N_r)
	elevSecondMomentFunction  = DBLARR(N)
	elevSecondMomentFunction =  FFT(totalElevAvgPower,/DOUBLE,/INVERSE)
	elevSecondMomentFunction /= elevSecondMomentFunction[0]
	
	;;; Compute the second order elev cumulant function
	elevSecondCumulantFunction  = DBLARR(N)
	elevSecondCumulantFunction = $
		(elevMoments[1]*elevSecondMomentFunction - elevMoments[0]^2.D)/elevCumulants[1]

	;;; Compute the bispectrum estimators
	elevBispectrum /= DOUBLE(N_r)
	elevComponentPower /= DOUBLE(N_r)
	elevSumPower /= DOUBLE(N_r)

	;;; Compute the bicoherence
	FOR j=0L,NN4 DO BEGIN
		FOR i=j,NN2-j DO BEGIN
			IF (SQRT(elevComponentPower[i,j])*SQRT(elevSumPower[i,j]) GT 10.D^(-12.D)) THEN BEGIN
				elevBicoherence[i,j] = ABS(elevBispectrum[i,j])/(SQRT(elevComponentPower[i,j])*SQRT(elevSumPower[i,j]))
			ENDIF ELSE BEGIN
				elevBicoherence[i,j] = 0.D
			ENDELSE
		ENDFOR
	ENDFOR

	;;; Fill the rest of the bispectrum and bicoherence array
	bispectrumSymmetry,elevBispectrum,NN
	bicoherenceSymmetry,elevBicoherence,NN

	;;; Compute the elevation third moment function
	elevThirdMomentFunction  = DBLARR(NN,NN)
	elevThirdMomentFunction =  FFT(elevBispectrum,/DOUBLE,/INVERSE)
	elevThirdMomentFunction /= elevThirdMomentFunction[0,0]

	;;; Compute the elevation third cumulant function
	elevThirdCumulantFunction  = DBLARR(NN,NN)
	FOR i=0L,NN/2L DO BEGIN
		FOR j=0L,i DO BEGIN
			elevThirdCumulantFunction[i,j] = (elevMoments[2]*elevThirdMomentFunction[i,j] $
				- elevMoments[0]*elevMoments[1]* $
				(elevSecondMomentFunction[i] + $
				elevSecondMomentFunction[j] + $
				elevSecondMomentFunction[ABS(j-i)]) $
				+ 2.D*elevMoments[0]^3.D)/elevCumulants[2]

			elevThirdCumulantFunction[j,i] = elevThirdCumulantFunction[i,j]
		ENDFOR
	ENDFOR

	;;; Fill the rest of the third cumulant function array
	biCovarianceSymmetry,elevThirdCumulantFunction,NN

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;;    Compute the slope estimators   ;;;
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	;;; Compute the estimators of the IDL moments
	slopeMean /= DOUBLE(N_r)
	slopeVariance /= DOUBLE(N_r)
	slopeSkewness /= DOUBLE(N_r)

	;;; Compute the average elevation moments
	slopeMoments /= DOUBLE(N_r)

	;;; Compute the average slope Cumulants
	slopeCumulants = DBLARR(numMoments)
	cumulantsFromMoments,slopeMoments,slopeCumulants

	;;; Compute the second order slope moment function and spectra
	totalSlopeAvgPower /= DOUBLE(N_r)
	slopeSecondMomentFunction = DCOMPLEXARR(N)
	slopeSecondMomentFunction = FFT(totalSlopeAvgPower,/DOUBLE,/INVERSE)
	slopeSecondMomentFunction /= slopeSecondMomentFunction[0]
	
	;;; Compute the second order slope cumulant function
	slopeSecondCumulantFunction  = DBLARR(N)
	slopeSecondCumulantFunction = $
		(slopeMoments[1]*slopeSecondMomentFunction - slopeMoments[0]^2.D)/slopeCumulants[1]
	
	;;; Compute the bispectrum estimators
	slopeBispectrum /= DOUBLE(N_r)
	slopeComponentPower /= DOUBLE(N_r)
	slopeSumPower /= DOUBLE(N_r)

	;;; Compute the bicoherence
	FOR j=0L,NN4 DO BEGIN
		FOR i=j,NN2-j DO BEGIN
			IF (SQRT(slopeComponentPower[i,j])*SQRT(slopeSumPower[i,j]) GT 10.D^(-12.D)) THEN BEGIN
				slopeBicoherence[i,j] = ABS(slopeBispectrum[i,j])/(SQRT(slopeComponentPower[i,j])*SQRT(slopeSumPower[i,j]))
			ENDIF ELSE BEGIN
				slopeBicoherence[i,j] = 0.D
			ENDELSE
		ENDFOR
	ENDFOR

	;;; Fill the rest of the bispectrum and bicoherence array
	bispectrumSymmetry,slopeBispectrum,NN
	bicoherenceSymmetry,slopeBicoherence,NN

	;;; Compute the slope third moment functions
	slopeThirdMomentFunction = DCOMPLEXARR(NN,NN)
	slopeThirdMomentFunction = FFT(slopeBispectrum,/DOUBLE,/INVERSE)
	slopeThirdMomentFunction /= slopeThirdMomentFunction[0,0]

	;;; Compute the slope third cumulant function
	slopeThirdCumulantFunction  = DBLARR(NN,NN)
	FOR i=0L,NN/2L DO BEGIN
		FOR j=0L,i DO BEGIN
			slopeThirdCumulantFunction[i,j] = (slopeMoments[2]*slopeThirdMomentFunction[i,j] $
				- slopeMoments[0]*slopeMoments[1]* $
				(slopeSecondMomentFunction[i] + $
				slopeSecondMomentFunction[j] + $
				slopeSecondMomentFunction[ABS(j-i)]) $
				+ 2.D*slopeMoments[0]^3.D)/slopeCumulants[2]

			slopeThirdCumulantFunction[j,i] = slopeThirdCumulantFunction[i,j]
		ENDFOR
	ENDFOR

	;;; Fill the rest of the third cumulant function array
	biCovarianceSymmetry,slopeThirdCumulantFunction,NN
	
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;;    Compute the glint estimators    ;;;
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	;;; Compute the estimators of the IDL moments
	glintMean /= DOUBLE(N_r)
	glintVariance /= DOUBLE(N_r)
	glintSkewness /= DOUBLE(N_r)

	;;; Compute the average elevation moments
	glintFirstMoments /= DOUBLE(angleRuns)

	;;; Compute the average glint Cumulants
	glintCumulants = DBLARR(numMoments,GEOM.N_angles)
	glintCumulantsFromMoments,glintFirstMoments,glintCumulants

	;;; Compute the second order glint moment function and spectra
	FOR angle=0L, GEOM.N_angles-1L DO BEGIN
		totalGlintAvgPower[*,angle] /= DOUBLE(angleRunsCum[angle])
	ENDFOR
	
	glintSecondMomentFunction = DCOMPLEXARR(N,GEOM.N_angles)
	FOR angle=0L, GEOM.N_angles-1L DO BEGIN
		glintSecondMomentFunction[*,angle] = FFT(totalGlintAvgPower[*,angle],/DOUBLE,/INVERSE)
		glintSecondMomentFunction[*,angle] /= glintSecondMomentFunction[0,angle]
	ENDFOR

	;;; Compute the second order glint cumulant functions
	glintSecondCumulantFunction  = DBLARR(N,GEOM.N_angles)
	FOR angle=0L, GEOM.N_angles-1L DO BEGIN
		glintSecondCumulantFunction[*,angle] = $
			(glintFirstMoments[angle]*glintSecondMomentFunction[*,angle] $
			- glintFirstMoments[angle]^2.D)/glintCumulants[1,angle]
	ENDFOR
	
	;;; Compute the bispectrum estimators
	FOR angle=0L, GEOM.N_angles-1L DO BEGIN
		glintBispectrum[*,*,angle] /= DOUBLE(angleRunsCum[angle])
		glintComponentPower[*,*,angle] /= DOUBLE(angleRunsCum[angle])
		glintSumPower[*,*,angle] /= DOUBLE(angleRunsCum[angle])
	ENDFOR

	;;; Compute the bicoherence
	FOR angle=0L, GEOM.N_angles-1L DO BEGIN
		FOR j=0L,NN4 DO BEGIN
			FOR i=j,NN2-j DO BEGIN
				IF (SQRT(glintComponentPower[i,j,angle])*SQRT(glintSumPower[i,j,angle]) GT 10.D^(-12.D)) THEN BEGIN
					glintBicoherence[i,j,angle] = $
						ABS(glintBispectrum[i,j,angle])/(SQRT(glintComponentPower[i,j,angle])*SQRT(glintSumPower[i,j,angle]))
				ENDIF ELSE BEGIN
					glintBicoherence[i,j,angle] = 0.D
				ENDELSE
			ENDFOR
		ENDFOR
	ENDFOR

	;;; Fill the rest of the bispectrum and bicoherence array
	FOR angle=0,GEOM.N_angles-1L DO BEGIN
		tempBispectrum = glintBispectrum[*,*,angle]
		bispectrumSymmetry,tempBispectrum,NN
		glintBispectrum[*,*,angle] = tempBispectrum 

		tempBicoherence = glintBicoherence[*,*,angle]
		bicoherenceSymmetry,tempBicoherence,NN
		glintBicoherence[*,*,angle] = tempBicoherence 
	ENDFOR
	
	;;; Compute the glint third moment functions
	glintThirdMomentFunction = DCOMPLEXARR(NN,NN,GEOM.N_angles)
	FOR angle=0L, GEOM.N_angles-1L DO BEGIN
		glintThirdMomentFunction[*,*,angle] = FFT(glintBispectrum[*,*,angle],/DOUBLE,/INVERSE)
		glintThirdMomentFunction[*,*,angle] /= glintThirdMomentFunction[0,0,angle]
	ENDFOR

	;;; Compute the glint third cumulant functions
	glintThirdCumulantFunction  = DBLARR(NN,NN,GEOM.N_angles)
	FOR angle=0,GEOM.N_angles-1L DO BEGIN
		FOR i=0L,NN/2L DO BEGIN
			FOR j=0L,i DO BEGIN
				glintThirdCumulantFunction[i,j,angle] = $
					(glintFirstMoments[angle]*glintThirdMomentFunction[i,j,angle] $
					- (glintFirstMoments[angle]^2.D)* $
						(glintSecondMomentFunction[i,angle] + $
						glintSecondMomentFunction[j,angle] + $
						glintSecondMomentFunction[ABS(j-i),angle]) $
					+ 2.D*glintFirstMoments[angle]^3.D)/glintCumulants[2,angle]

				glintThirdCumulantFunction[j,i,angle] = glintThirdCumulantFunction[i,j,angle]
			ENDFOR
		ENDFOR
	ENDFOR

	;;; Fill the rest of the third cumulant function array
	FOR angle=0,GEOM.N_angles-1L DO BEGIN
		tempThirdCumulantFunction = glintThirdCumulantFunction[*,*,angle]
		biCovarianceSymmetry,tempThirdCumulantFunction,NN
		glintThirdCumulantFunction[*,*,angle] = tempThirdCumulantFunction
	ENDFOR

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;; The elevation and slope summary results                   ;;;
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	print ''
	print 'IDL Elevation Moments:   ',elevMean,elevVariance,elevSkewness*(elevVariance^1.5D),FORMAT='(A,3E16.7)'
	print '    Elevation Moments:   ',elevMoments[0],elevMoments[1],elevMoments[2],FORMAT='(A,3E16.7)'
	print '    Elevation Cumulants: ',elevCumulants[0],elevCumulants[1],elevCumulants[2],FORMAT='(A,3E16.7)'
	print ''
	print 'IDL Slope Moments:   ',slopeMean,slopeVariance,slopeSkewness*(slopeVariance^1.5D),FORMAT='(A,3E16.7)'
	print '    Slope Moments:   ',slopeMoments[0],slopeMoments[1],slopeMoments[2],FORMAT='(A,3E16.7)'
	print '    Slope Cumulants: ',slopeCumulants[0],slopeCumulants[1],slopeCumulants[2],FORMAT='(A,3E16.7)'
	print ''
	;print 'IDL Glint Moments:   ',glintMean,glintVariance,glintSkewness*(glintVariance^1.5D),FORMAT='(A,3E16.7)'
	print 'Glint First Moments:   ',glintFirstMoments,FORMAT='(A/,'+STRING(GEOM.N_angles)+'E16.7/)'
	print 'Glint Cumulants:       ',TRANSPOSE(glintCumulants),FORMAT='(A/,3('+STRING(GEOM.N_angles)+'E16.7/))'

	print ''
	print "Elevation third moment from bicovariance: ",DOUBLE(elevThirdMomentFunction[0L,0L]),FORMAT='(A,E16.7)'
	print "  Elevation third moment from bispectrum: ",TOTAL(DOUBLE(elevBispectrum)),FORMAT='(A,E16.7)'
	print ''
	print "    Slope third moment from bicovariance: ",DOUBLE(slopeThirdMomentFunction[0L,0L]),FORMAT='(A,E16.7)'
	print "      Slope third moment from bispectrum: ",TOTAL(DOUBLE(slopeBispectrum)),FORMAT='(A,E16.7)'
	print ''
	print "    glint third moment from bicovariance: ",DOUBLE(glintThirdMomentFunction[0L,0L,*]),FORMAT='(A,5E16.7)'
	print "      glint third moment from bispectrum: "
	FOR angle=0L, GEOM.N_angles-1L DO BEGIN
		print  TOTAL(DOUBLE(glintBispectrum[*,*,angle])),FORMAT='(E16.7)'
	ENDFOR
	print ''

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;;   Open the output HDF file   ;;;
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	fileName = 'GlintSim_'+ $
				STRING(N)+'_'+ $
				STRING(NN)+'_'+ $
				STRING(FIX(100.*delta_x))+'_'+ $
				STRING(N_runs)+'_'+ $
				STRING(spectrumType)+'_'+ $
				STRING(FIX(specExp))+'_'+ $
				STRING(powWindType)+'_'+ $
				STRING(bispWindType)

	IF (nlSwitch EQ 1L) THEN BEGIN
		fileName += '_nl.hdf'
	ENDIF ELSE BEGIN
		fileName += '.hdf'
	ENDELSE

	fileName = STRCOMPRESS(fileName,/REMOVE_ALL)

	print "Open output filename: ",filename

	outFile = fileinfo(fileName)
	help, fileinfo(fileName), /structure
	print  outFile

	IF (NOT outFile.EXIST) THEN BEGIN
		;;; Create and open file using SD interface
		fileID = HDF_SD_START(fileName, /CREATE)
		;fileID = HDF_OPEN(fileName, /CREATE,/WRITE)
		print  'Created new HDF file: ',fileName
	ENDIF ELSE BEGIN
		;;; Create and open file using SD interface
		print  'HDF file ',fileName,' exists, opening...'
		fileID = HDF_SD_START(fileName, /RdWr)
		;fileID = HDF_OPEN(fileName, /WRITE)
		print  'Opened HDF file ',fileName,' for reading and writing'
	ENDELSE

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;;   Add some attributes to the file   ;;;
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	
	IF (NOT outFile.EXIST) THEN BEGIN
		HDF_SD_ATTRSET, fileID, 'DATE', SYSTIME()
		HDF_SD_ATTRSET, fileID, 'EXPERIMENT', 'cumulantFunctionSimulate.pro'
		HDF_SD_ATTRSET, fileID, 'NAME', 'Geoff Cureton'
		HDF_SD_ATTRSET, fileID, 'EMAIL ADDRESS', 'geoff.cureton@physics.org'
	ENDIF ELSE BEGIN
		HDF_SD_ATTRSET, fileID, 'DATE', SYSTIME()+" Zulu"
	ENDELSE

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;;   Add some datasets to the file   ;;;
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	   
	IF (NOT outFile.EXIST) THEN BEGIN

		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		;;;   Save the geometry information to global attributes, and variables   ;;;
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

		print  'Writing geometry angles and slopes...'

		sourceAngleID  = HDF_SD_CREATE(fileID, "Solar Zenith Angles", [GEOM.N_angles], /FLOAT)
		HDF_SD_ADDDATA, sourceAngleID, GEOM.source_angle
		numAnglesDimID = HDF_SD_DIMGETID(sourceAngleID, 0)
		HDF_SD_DIMSET, numAnglesDimID, LABEL='Number of Angles', NAME='N_angles', UNIT='radians'
		HDF_SD_ENDACCESS, sourceAngleID 

		detectorAngleID  = HDF_SD_CREATE(fileID, "Detector Zenith Angles", [GEOM.N_angles], /FLOAT)
		HDF_SD_ADDDATA, detectorAngleID, GEOM.detector_angle
		numAnglesDimID = HDF_SD_DIMGETID(detectorAngleID, 0)
		HDF_SD_DIMSET, numAnglesDimID, LABEL='Number of Angles', NAME='N_angles', UNIT='radians'
		HDF_SD_ENDACCESS, detectorAngleID 

		specularSlopeID  = HDF_SD_CREATE(fileID, "Specular Slopes", [GEOM.N_angles], /FLOAT)
		HDF_SD_ADDDATA, specularSlopeID, GEOM.xi_0
		numAnglesDimID = HDF_SD_DIMGETID(specularSlopeID, 0)
		HDF_SD_DIMSET, numAnglesDimID, LABEL='Number of Angles', NAME='N_angles', UNIT='radians'
		HDF_SD_ENDACCESS, specularSlopeID 

		minSlopeID  = HDF_SD_CREATE(fileID, "Min Slopes", [GEOM.N_angles], /FLOAT)
		HDF_SD_ADDDATA, minSlopeID, GEOM.xi_min
		numAnglesDimID = HDF_SD_DIMGETID(minSlopeID, 0)
		HDF_SD_DIMSET, numAnglesDimID, LABEL='Number of Angles', NAME='N_angles', UNIT='radians'
		HDF_SD_ENDACCESS, minSlopeID 

		maxSlopeID  = HDF_SD_CREATE(fileID, "Max Slopes", [GEOM.N_angles], /FLOAT)
		HDF_SD_ADDDATA, maxSlopeID, GEOM.xi_max
		numAnglesDimID = HDF_SD_DIMGETID(maxSlopeID, 0)
		HDF_SD_DIMSET, numAnglesDimID, LABEL='Number of Angles', NAME='N_angles', UNIT='radians'
		HDF_SD_ENDACCESS, maxSlopeID 

		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		;;;   Save the elevation, slope and glint moments   ;;;
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		
		print  'Writing the elevation, slope and glint moments...'

		elevMomentID  = HDF_SD_CREATE(fileID, 'Elevation Moments', [numMoments],    /FLOAT)
		HDF_SD_ADDDATA, elevMomentID , elevMoments
		numMomentsDimID = HDF_SD_DIMGETID(elevMomentID, 0)
		HDF_SD_DIMSET, numMomentsDimID, LABEL='Number of Moments', NAME='N_moments'
		HDF_SD_ENDACCESS, elevMomentID 
		
		slopeMomentID = HDF_SD_CREATE(fileID, 'Slope Moments',     [numMoments],    /FLOAT) 
		HDF_SD_ADDDATA, slopeMomentID, slopeMoments
		numMomentsDimID = HDF_SD_DIMGETID(slopeMomentID, 0)
		HDF_SD_DIMSET, numMomentsDimID, LABEL='Number of Moments', NAME='N_moments'
		HDF_SD_ENDACCESS, slopeMomentID
		
		glintMomentID = HDF_SD_CREATE(fileID, 'Glint Moments',     [GEOM.N_angles], /FLOAT)
		HDF_SD_ADDDATA, glintMomentID, glintFirstMoments
		numAnglesDimID = HDF_SD_DIMGETID(glintMomentID, 0)
		HDF_SD_DIMSET, numAnglesDimID, LABEL='Number of Angles', NAME='N_angles', UNIT='radians'
		HDF_SD_ENDACCESS, glintMomentID
		
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		;;;   Save the elevation, slope and glint cumulants ;;;
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		
		print  'Writing the elevation, slope and glint cumulants...'

		elevCumulantID  = HDF_SD_CREATE(fileID, 'Elevation Cumulants', [numMoments],    /FLOAT)
		HDF_SD_ADDDATA, elevCumulantID , elevCumulants
		numCumulantsDimID = HDF_SD_DIMGETID(elevCumulantID, 0)
		HDF_SD_DIMSET, numCumulantsDimID, LABEL='Number of Cumulants', NAME='N_moments'
		HDF_SD_ENDACCESS, elevCumulantID 
		
		slopeCumulantID = HDF_SD_CREATE(fileID, 'Slope Cumulants',     [numMoments],    /FLOAT) 
		HDF_SD_ADDDATA, slopeCumulantID, slopeCumulants
		numCumulantsDimID = HDF_SD_DIMGETID(slopeCumulantID, 0)
		HDF_SD_DIMSET, numCumulantsDimID, LABEL='Number of Cumulants', NAME='N_moments'
		HDF_SD_ENDACCESS, slopeCumulantID
		
		glintCumulantID = HDF_SD_CREATE(fileID, 'Glint Cumulants',     [numMoments,GEOM.N_angles], /FLOAT)
		HDF_SD_ADDDATA, glintCumulantID, glintCumulants
		numCumulantsDimID = HDF_SD_DIMGETID(glintCumulantID, 0)
		HDF_SD_DIMSET, numCumulantsDimID, LABEL='Number of Cumulants', NAME='N_moments'
		numAnglesDimID = HDF_SD_DIMGETID(glintCumulantID, 1)
		HDF_SD_DIMSET, numAnglesDimID, LABEL='Number of Angles', NAME='N_angles', UNIT='radians'
		HDF_SD_ENDACCESS, glintCumulantID
		
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		;;;   Save the indices of the interating modes   ;;;
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		
		print  'Writing indices of nonlinear components...'

		free1ID  = HDF_SD_CREATE(fileID, "Indices of first free wavenumber modes", [N_ELEMENTS(NLCOUPLING.free1)], /LONG)
		HDF_SD_ADDDATA, free1ID, NLCOUPLING.free1
		NLMODE_indexDimID = HDF_SD_DIMGETID(free1ID, 0)
		HDF_SD_DIMSET, NLMODE_indexDimID, LABEL='Number of interacting mode triplets', NAME='NumModeTriplets'
		HDF_SD_ENDACCESS, free1ID 
		
		free2ID  = HDF_SD_CREATE(fileID, "Indices of second free wavenumber modes", [N_ELEMENTS(NLCOUPLING.free2)], /LONG)
		HDF_SD_ADDDATA, free2ID, NLCOUPLING.free2
		NLMODE_indexDimID = HDF_SD_DIMGETID(free2ID, 0)
		HDF_SD_DIMSET, NLMODE_indexDimID, LABEL='Number of interacting mode triplets', NAME='NumModeTriplets'
		HDF_SD_ENDACCESS, free2ID 

		boundID  = HDF_SD_CREATE(fileID, "Indices of bound wavenumber modes", [N_ELEMENTS(NLCOUPLING.bound)], /LONG)
		HDF_SD_ADDDATA, boundID, NLCOUPLING.bound
		NLMODE_indexDimID = HDF_SD_DIMGETID(boundID, 0)
		HDF_SD_DIMSET, NLMODE_indexDimID, LABEL='Number of interacting mode triplets', NAME='NumModeTriplets'
		HDF_SD_ENDACCESS, boundID 
		
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		;;;   Set 1D wavenumber scale   ;;;
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		
		wavenumber = FINDGEN(N)*delta_k
		wavenumberID  = HDF_SD_CREATE(fileID, "Power Spectrum wavenumber scale", [N], /FLOAT)
		HDF_SD_ADDDATA, wavenumberID,  wavenumber
		HDF_SD_ATTRSET, wavenumberID,  'units', 'meters^{-1}'
		HDF_SD_ATTRSET, wavenumberID,  'increment', delta_k
		wavenumberDimID = HDF_SD_DIMGETID(wavenumberID, 0)
		HDF_SD_DIMSET, wavenumberDimID, LABEL='Data length', NAME='N', SCALE=wavenumber, UNIT='meters^{-1}'
		HDF_SD_ENDACCESS, wavenumberID 

		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		;;;   Save the average power spectra   ;;;
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

		print  'Writing average power spectra ...'

		;;; Elevation Power Spectrum
		
		elevPowID  = HDF_SD_CREATE(fileID, "Average Elevation Power Spectrum", [N], /FLOAT)
		HDF_SD_ADDDATA, elevPowID,  totalElevAvgPower
		HDF_SD_ATTRSET, elevPowID,  'long_name', $
			'Realisation averaged elevation power spectrum, power[0 ... N-1]'
		powerSpectrumDimID = HDF_SD_DIMGETID(elevPowID, 0)
		HDF_SD_DIMSET, powerSpectrumDimID, LABEL='Data length', NAME='N', SCALE=wavenumber, UNIT='meters^{-1}'
		HDF_SD_ENDACCESS, elevPowID 
		
		;;; Slope Power Spectrum
		
		slopePowID = HDF_SD_CREATE(fileID, "Average Slope Power Spectrum"    , [N], /FLOAT)
		HDF_SD_ADDDATA, slopePowID, totalSlopeAvgPower
		HDF_SD_ATTRSET, slopePowID, 'long_name', $
			'Realisation averaged slope power spectrum, power[0 ... N-1]'
		powerSpectrumDimID = HDF_SD_DIMGETID(slopePowID, 0)
		HDF_SD_DIMSET, powerSpectrumDimID, LABEL='Data length', NAME='N', SCALE=wavenumber, UNIT='meters^{-1}'
		HDF_SD_ENDACCESS, slopePowID
		
		;;; Glint Power Spectrum
		
		glintPowID = HDF_SD_CREATE(fileID, "Average Glint Power Spectrum"    , [N,GEOM.N_angles], /FLOAT)
		HDF_SD_ADDDATA, glintPowID, totalGlintAvgPower
		HDF_SD_ATTRSET, glintPowID, 'long_name', $
			'Realisation averaged glint power spectrum, power[0 ... N-1][0 ... N_angles-1]'
		powerSpectrumDimID = HDF_SD_DIMGETID(glintPowID, 0)
		HDF_SD_DIMSET, powerSpectrumDimID, LABEL='Data length', NAME='N', SCALE=wavenumber, UNIT='meters^{-1}'
		numAnglesDimID = HDF_SD_DIMGETID(glintPowID, 1)
		HDF_SD_DIMSET, numAnglesDimID, LABEL='Number of Angles', NAME='N_angles', UNIT='radians'
		HDF_SD_ENDACCESS, glintPowID

		;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		;;;   Set 1D lag scale   ;;;
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		
		length = FINDGEN(N)*delta_x
		lengthID  = HDF_SD_CREATE(fileID, "Cumulant Function 1D length scale", [N], /FLOAT)
		HDF_SD_ADDDATA, lengthID,  length
		HDF_SD_ATTRSET, lengthID,  'units', 'meters'
		HDF_SD_ATTRSET, lengthID,  'increment', delta_x
		lengthDimID = HDF_SD_DIMGETID(lengthID, 0)
		HDF_SD_DIMSET, lengthDimID, LABEL='Data length', NAME='N', SCALE=length, UNIT='meters'
		HDF_SD_ENDACCESS, lengthID 

		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		;;;   Save the Second Moment Functions   ;;;
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		
		print  'Writing average second moment functions ...'
		
		;;; Elevation Second Moment Function

		elevSecondMomentFunctionID  = HDF_SD_CREATE(fileID, "Average Elevation Second Moment Function", [N], /FLOAT)
		HDF_SD_ADDDATA, elevSecondMomentFunctionID , elevSecondMomentFunction
		HDF_SD_ATTRSET, elevSecondMomentFunctionID , 'long_name', $
			'Realisation averaged elevation Second Moment Function, secondMomentFunction[0 ... N-1]'
		secondMomentFunctionDimID = HDF_SD_DIMGETID(elevSecondMomentFunctionID, 0)
		HDF_SD_DIMSET, secondMomentFunctionDimID, LABEL='Second Moment Function Data length', NAME='N', SCALE=length, UNIT='meters'
		HDF_SD_ENDACCESS, elevSecondMomentFunctionID 
		
		;;; Slope Second Moment Function
		
		slopeSecondMomentFunctionID = HDF_SD_CREATE(fileID, "Average Slope Second Moment Function"    , [N], /FLOAT)
		HDF_SD_ADDDATA, slopeSecondMomentFunctionID, slopeSecondMomentFunction
		HDF_SD_ATTRSET, slopeSecondMomentFunctionID, 'long_name', $
			'Realisation averaged slope Second Moment Function, secondMomentFunction[0 ... N-1]'
		secondMomentFunctionDimID = HDF_SD_DIMGETID(slopeSecondMomentFunctionID, 0)
		HDF_SD_DIMSET, secondMomentFunctionDimID, LABEL='Second Moment Function Data length', NAME='N', SCALE=length, UNIT='meters'
		HDF_SD_ENDACCESS, slopeSecondMomentFunctionID
		
		;;; Glint Second Moment Function
		
		glintSecondMomentFunctionID = HDF_SD_CREATE(fileID, "Average Glint Second Moment Function"    , [N, GEOM.N_angles], /FLOAT)
		HDF_SD_ADDDATA, glintSecondMomentFunctionID, glintSecondMomentFunction
		HDF_SD_ATTRSET, glintSecondMomentFunctionID, 'long_name', $
			'Realisation averaged glint Second Moment Function, secondMomentFunction[0 ... N-1][0 ... N_angles-1]'
		secondMomentFunctionDimID = HDF_SD_DIMGETID(glintSecondMomentFunctionID, 0)
		HDF_SD_DIMSET, secondMomentFunctionDimID, LABEL='Second Moment Function Data length', NAME='N', SCALE=length, UNIT='meters'
		numAnglesDimID = HDF_SD_DIMGETID(glintSecondMomentFunctionID, 1)
		HDF_SD_DIMSET, numAnglesDimID, LABEL='Number of Angles', NAME='N_angles', UNIT='radians'
		HDF_SD_ENDACCESS, glintSecondMomentFunctionID

		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		;;;   Save the Second Cumulant Functions   ;;;
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		
		print  'Writing average second cumulant functions ...'
		
		;;; Elevation Second Cumulant Function

		elevSecondCumulantFunctionID  = HDF_SD_CREATE(fileID, "Average Elevation Second Cumulant Function", [N], /FLOAT)
		HDF_SD_ADDDATA, elevSecondCumulantFunctionID , elevSecondCumulantFunction
		HDF_SD_ATTRSET, elevSecondCumulantFunctionID , 'long_name', $
			'Realisation averaged elevation Second Cumulant Function, secondCumulantFunction[0 ... N-1]'
		secondCumulantFunctionDimID = HDF_SD_DIMGETID(elevSecondCumulantFunctionID, 0)
		HDF_SD_DIMSET, secondCumulantFunctionDimID, LABEL='Second Cumulant Function Data length', NAME='N', SCALE=length, UNIT='meters'
		HDF_SD_ENDACCESS, elevSecondCumulantFunctionID 
		
		;;; Slope Second Cumulant Function
		
		slopeSecondCumulantFunctionID = HDF_SD_CREATE(fileID, "Average Slope Second Cumulant Function"    , [N], /FLOAT)
		HDF_SD_ADDDATA, slopeSecondCumulantFunctionID, slopeSecondCumulantFunction
		HDF_SD_ATTRSET, slopeSecondCumulantFunctionID, 'long_name', $
			'Realisation averaged slope Second Cumulant Function, secondCumulantFunction[0 ... N-1]'
		secondCumulantFunctionDimID = HDF_SD_DIMGETID(slopeSecondCumulantFunctionID, 0)
		HDF_SD_DIMSET, secondCumulantFunctionDimID, LABEL='Second Cumulant Function Data length', NAME='N', SCALE=length, UNIT='meters'
		HDF_SD_ENDACCESS, slopeSecondCumulantFunctionID
		
		;;; Glint Second Cumulant Function
		
		glintSecondCumulantFunctionID = HDF_SD_CREATE(fileID, "Average Glint Second Cumulant Function"    , [N, GEOM.N_angles], /FLOAT)
		HDF_SD_ADDDATA, glintSecondCumulantFunctionID, glintSecondCumulantFunction
		HDF_SD_ATTRSET, glintSecondCumulantFunctionID, 'long_name', $
			'Realisation averaged glint Second Cumulant Function, secondCumulantFunction[0 ... N-1][0 ... N_angles-1]'
		secondCumulantFunctionDimID = HDF_SD_DIMGETID(glintSecondCumulantFunctionID, 0)
		HDF_SD_DIMSET, secondCumulantFunctionDimID, LABEL='Second Cumulant Function Data length', NAME='N', SCALE=length, UNIT='meters'
		numAnglesDimID = HDF_SD_DIMGETID(glintSecondCumulantFunctionID, 1)
		HDF_SD_DIMSET, numAnglesDimID, LABEL='Number of Angles', NAME='N_angles', UNIT='radians'
		HDF_SD_ENDACCESS, glintSecondCumulantFunctionID

		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		;;;   Set 2D wavenumber scale   ;;;
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		
		wavenumber2 = FINDGEN(NN)*delta_k
		wavenumber2ID  = HDF_SD_CREATE(fileID, "Bispectrum wavenumber scale", [NN], /FLOAT)
		HDF_SD_ADDDATA, wavenumber2ID,  wavenumber2
		HDF_SD_ATTRSET, wavenumber2ID,  'units', 'meters^{-1}'
		HDF_SD_ATTRSET, wavenumber2ID,  'increment', delta_k
		wavenumber2DimID = HDF_SD_DIMGETID(wavenumber2ID, 0)
		HDF_SD_DIMSET, wavenumber2DimID, LABEL='Data length', NAME='NN', SCALE=wavenumber2, UNIT='meters^{-1}'
		HDF_SD_ENDACCESS, wavenumber2ID
		
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		;;;   Save the Average Bispectra   ;;;
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		
		print  'Writing average bispectra ...'
		tempArr = DBLARR(NN,NN,2)
		tempGlintArr = DBLARR(NN,NN,2,GEOM.N_angles)
		
		;;; Elevation Bispectrum
		tempArr[*,*,0] = DOUBLE(elevBispectrum)
		tempArr[*,*,1] = IMAGINARY(elevBispectrum)

		elevBispectrumID  = HDF_SD_CREATE(fileID, "Average Elevation Bispectrum", [NN,NN,2], /FLOAT)
		HDF_SD_ADDDATA, elevBispectrumID , tempArr
		HDF_SD_ATTRSET, elevBispectrumID , 'long_name', $
			'Realisation averaged elevation bispectrum, bispectrum[0 ... NN-1][0 ... NN-1]'
		bispectrumDimID = HDF_SD_DIMGETID(elevBispectrumID, 0)
		HDF_SD_DIMSET, bispectrumDimID, LABEL='Bispectrum Data length', NAME='NN' ,SCALE=wavenumber2,  UNIT='meters^{-1}'
		bispectrumDimID = HDF_SD_DIMGETID(elevBispectrumID, 1)
		HDF_SD_DIMSET, bispectrumDimID, LABEL='Bispectrum Data length', NAME='NN' ,SCALE=wavenumber2,  UNIT='meters^{-1}'
		HDF_SD_ENDACCESS, elevBispectrumID 
		
		;;; Slope Bispectrum
		tempArr[*,*,0] = DOUBLE(slopeBispectrum)
		tempArr[*,*,1] = IMAGINARY(slopeBispectrum)
		
		slopeBispectrumID = HDF_SD_CREATE(fileID, "Average Slope Bispectrum"    , [NN,NN,2], /FLOAT)
		HDF_SD_ADDDATA, slopeBispectrumID, tempArr
		HDF_SD_ATTRSET, slopeBispectrumID, 'long_name', $
			'Realisation averaged slope bispectrum, bispectrum[0 ... NN-1][0 ... NN-1]'
		bispectrumDimID = HDF_SD_DIMGETID(slopeBispectrumID, 0)
		HDF_SD_DIMSET, bispectrumDimID, LABEL='Bispectrum Data length', NAME='NN' , SCALE=wavenumber2, UNIT='meters^{-1}'
		bispectrumDimID = HDF_SD_DIMGETID(slopeBispectrumID, 1)
		HDF_SD_DIMSET, bispectrumDimID, LABEL='Bispectrum Data length', NAME='NN' , SCALE=wavenumber2, UNIT='meters^{-1}'
		HDF_SD_ENDACCESS, slopeBispectrumID
		
		;;; Real Glint Bispectrum
		help,DOUBLE(glintBispectrum)

		glintBispectrumID = HDF_SD_CREATE(fileID, "Average Glint Real Bispectrum"    , [NN,NN,GEOM.N_angles], /FLOAT)
		HDF_SD_ADDDATA, glintBispectrumID, DOUBLE(glintBispectrum)
		HDF_SD_ATTRSET, glintBispectrumID, 'long_name', $
			'Realisation averaged real glint bispectrum, bispectrum[0 ... NN-1][0 ... NN-1][0 ... N_angles-1]'
		bispectrumDimID = HDF_SD_DIMGETID(glintBispectrumID, 0)
		HDF_SD_DIMSET, bispectrumDimID, LABEL='Bispectrum Data length', NAME='NN' , SCALE=wavenumber2, UNIT='meters^{-1}'
		bispectrumDimID = HDF_SD_DIMGETID(glintBispectrumID, 1)
		HDF_SD_DIMSET, bispectrumDimID, LABEL='Bispectrum Data length', NAME='NN' , SCALE=wavenumber2, UNIT='meters^{-1}'
		numAnglesDimID = HDF_SD_DIMGETID(glintBispectrumID, 2)
		HDF_SD_DIMSET, numAnglesDimID, LABEL='Number of Angles', NAME='N_angles', UNIT='radians'
		HDF_SD_ENDACCESS, glintBispectrumID

		;;; Imaginary Glint Bispectrum
		help,IMAGINARY(glintBispectrum)

		glintBispectrumID = HDF_SD_CREATE(fileID, "Average Glint Imaginary Bispectrum"    , [NN,NN,GEOM.N_angles], /FLOAT)
		HDF_SD_ADDDATA, glintBispectrumID, IMAGINARY(glintBispectrum)
		HDF_SD_ATTRSET, glintBispectrumID, 'long_name', $
			'Realisation averaged glint imaginary bispectrum, bispectrum[0 ... NN-1][0 ... NN-1][0 ... N_angles-1]'
		bispectrumDimID = HDF_SD_DIMGETID(glintBispectrumID, 0)
		HDF_SD_DIMSET, bispectrumDimID, LABEL='Bispectrum Data length', NAME='NN' , SCALE=wavenumber2, UNIT='meters^{-1}'
		bispectrumDimID = HDF_SD_DIMGETID(glintBispectrumID, 1)
		HDF_SD_DIMSET, bispectrumDimID, LABEL='Bispectrum Data length', NAME='NN' , SCALE=wavenumber2, UNIT='meters^{-1}'
		numAnglesDimID = HDF_SD_DIMGETID(glintBispectrumID, 2)
		HDF_SD_DIMSET, numAnglesDimID, LABEL='Number of Angles', NAME='N_angles', UNIT='radians'
		HDF_SD_ENDACCESS, glintBispectrumID

		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		;;;   Save the Average Bicoherence ;;;
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		
		print  'Writing average bicoherence ...'
		
		;;; Elevation Bicoherence
		
		elevBicoherenceID  = HDF_SD_CREATE(fileID, "Average Elevation Bicoherence", [NN,NN], /FLOAT)
		HDF_SD_ADDDATA, elevBicoherenceID , elevBicoherence
		HDF_SD_ATTRSET, elevBicoherenceID , 'long_name', $
			'Realisation averaged elevation Bicoherence, Bicoherence[0 ... NN-1][0 ... NN-1]'
		BicoherenceDimID = HDF_SD_DIMGETID(elevBicoherenceID, 0)
		HDF_SD_DIMSET, BicoherenceDimID, LABEL='Bicoherence Data length', NAME='NN' ,SCALE=wavenumber2,  UNIT='meters^{-1}'
		BicoherenceDimID = HDF_SD_DIMGETID(elevBicoherenceID, 1)
		HDF_SD_DIMSET, BicoherenceDimID, LABEL='Bicoherence Data length', NAME='NN' ,SCALE=wavenumber2,  UNIT='meters^{-1}'
		HDF_SD_ENDACCESS, elevBicoherenceID 
		
		;;; Slope Bicoherence
		
		slopeBicoherenceID = HDF_SD_CREATE(fileID, "Average Slope Bicoherence"    , [NN,NN], /FLOAT)
		HDF_SD_ADDDATA, slopeBicoherenceID, slopeBicoherence
		HDF_SD_ATTRSET, slopeBicoherenceID, 'long_name', $
			'Realisation averaged slope Bicoherence, Bicoherence[0 ... NN-1][0 ... NN-1]'
		BicoherenceDimID = HDF_SD_DIMGETID(slopeBicoherenceID, 0)
		HDF_SD_DIMSET, BicoherenceDimID, LABEL='Bicoherence Data length', NAME='NN' , SCALE=wavenumber2, UNIT='meters^{-1}'
		BicoherenceDimID = HDF_SD_DIMGETID(slopeBicoherenceID, 1)
		HDF_SD_DIMSET, BicoherenceDimID, LABEL='Bicoherence Data length', NAME='NN' , SCALE=wavenumber2, UNIT='meters^{-1}'
		HDF_SD_ENDACCESS, slopeBicoherenceID
		
		;;; Glint Bicoherence

		glintBicoherenceID = HDF_SD_CREATE(fileID, "Average Glint Bicoherence"    , [NN,NN,GEOM.N_angles], /FLOAT)
		HDF_SD_ADDDATA, glintBicoherenceID, glintBicoherence
		HDF_SD_ATTRSET, glintBicoherenceID, 'long_name', $
			'Realisation averaged glint Bicoherence, Bicoherence[0 ... NN-1][0 ... NN-1][0 ... N_angles-1]'
		BicoherenceDimID = HDF_SD_DIMGETID(glintBicoherenceID, 0)
		HDF_SD_DIMSET, BicoherenceDimID, LABEL='Bicoherence Data length', NAME='NN' , SCALE=wavenumber2, UNIT='meters^{-1}'
		BicoherenceDimID = HDF_SD_DIMGETID(glintBicoherenceID, 1)
		HDF_SD_DIMSET, BicoherenceDimID, LABEL='Bicoherence Data length', NAME='NN' , SCALE=wavenumber2, UNIT='meters^{-1}'
		numAnglesDimID = HDF_SD_DIMGETID(glintBicoherenceID, 2)
		HDF_SD_DIMSET, numAnglesDimID, LABEL='Number of Angles', NAME='N_angles', UNIT='radians'
		HDF_SD_ENDACCESS, glintBicoherenceID

		;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		;;;   Set 2D lag scale   ;;;
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		
		length2 = FINDGEN(NN)*delta_x
		length2ID  = HDF_SD_CREATE(fileID, "Cumulant Function 2D length scale", [NN], /FLOAT)
		HDF_SD_ADDDATA, length2ID,  length2
		HDF_SD_ATTRSET, length2ID,  'units', 'meters'
		HDF_SD_ATTRSET, length2ID,  'increment', delta_x
		length2DimID = HDF_SD_DIMGETID(length2ID, 0)
		HDF_SD_DIMSET, length2DimID, LABEL='Data length', NAME='NN', SCALE=length2, UNIT='meters'
		HDF_SD_ENDACCESS, length2ID 

		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		;;;   Save the Third Moment Functions   ;;;
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		
		print  'Writing average third moment functions ...'
		
		;;; Elevation Third Moment Function

		elevThirdMomentFunctionID  = HDF_SD_CREATE(fileID, "Average Elevation Third Moment Function", [NN,NN], /FLOAT)
		HDF_SD_ADDDATA, elevThirdMomentFunctionID , elevThirdMomentFunction
		HDF_SD_ATTRSET, elevThirdMomentFunctionID , 'long_name', $
			'Realisation averaged elevation Third Moment Function, ThirdMomentFunction[0 ... NN-1][0 ... NN-1]'
		thirdMomentFunctionDimID = HDF_SD_DIMGETID(elevThirdMomentFunctionID, 0)
		HDF_SD_DIMSET, thirdMomentFunctionDimID, LABEL='Third Moment Function Data length', NAME='NN', SCALE=length2, UNIT='meters'
		thirdMomentFunctionDimID = HDF_SD_DIMGETID(elevThirdMomentFunctionID, 1)
		HDF_SD_DIMSET, thirdMomentFunctionDimID, LABEL='Third Moment Function Data length', NAME='NN', SCALE=length2, UNIT='meters'
		HDF_SD_ENDACCESS, elevThirdMomentFunctionID 
		
		;;; Slope Third Moment Function
		
		slopeThirdMomentFunctionID = HDF_SD_CREATE(fileID, "Average Slope Third Moment Function"    , [NN,NN], /FLOAT)
		HDF_SD_ADDDATA, slopeThirdMomentFunctionID, slopeThirdMomentFunction
		HDF_SD_ATTRSET, slopeThirdMomentFunctionID, 'long_name', $
			'Realisation averaged slope Third Moment Function, ThirdMomentFunction[0 ... NN-1][0 ... NN-1]'
		thirdMomentFunctionDimID = HDF_SD_DIMGETID(slopeThirdMomentFunctionID, 0)
		HDF_SD_DIMSET, thirdMomentFunctionDimID, LABEL='Third Moment Function Data length', NAME='NN', SCALE=length2, UNIT='meters'
		thirdMomentFunctionDimID = HDF_SD_DIMGETID(slopeThirdMomentFunctionID, 1)
		HDF_SD_DIMSET, thirdMomentFunctionDimID, LABEL='Third Moment Function Data length', NAME='NN', SCALE=length2, UNIT='meters'
		HDF_SD_ENDACCESS, slopeThirdMomentFunctionID
		
		;;; Glint Third Moment Function
		
		glintThirdMomentFunctionID = HDF_SD_CREATE(fileID, "Average Glint Third Moment Function"    , [NN,NN, GEOM.N_angles], /FLOAT)
		HDF_SD_ADDDATA, glintThirdMomentFunctionID, glintThirdMomentFunction
		HDF_SD_ATTRSET, glintThirdMomentFunctionID, 'long_name', $
			'Realisation averaged glint Third Moment Function, ThirdMomentFunction[0 ... NN-1][0 ... NN-1][0 ... N_angles-1]'
		thirdMomentFunctionDimID = HDF_SD_DIMGETID(glintThirdMomentFunctionID, 0)
		HDF_SD_DIMSET, thirdMomentFunctionDimID, LABEL='Third Moment Function Data length', NAME='NN', SCALE=length2, UNIT='meters'
		thirdMomentFunctionDimID = HDF_SD_DIMGETID(glintThirdMomentFunctionID, 1)
		HDF_SD_DIMSET, thirdMomentFunctionDimID, LABEL='Third Moment Function Data length', NAME='NN', SCALE=length2, UNIT='meters'
		numAnglesDimID = HDF_SD_DIMGETID(glintThirdMomentFunctionID, 2)
		HDF_SD_DIMSET, numAnglesDimID, LABEL='Number of Angles', NAME='N_angles', UNIT='radians'
		HDF_SD_ENDACCESS, glintThirdMomentFunctionID

		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		;;;   Save the Third Cumulant Functions   ;;;
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		
		print  'Writing average third Cumulant functions ...'
		
		;;; Elevation Third Cumulant Function

		elevThirdCumulantFunctionID  = HDF_SD_CREATE(fileID, "Average Elevation Third Cumulant Function", [NN,NN], /FLOAT)
		HDF_SD_ADDDATA, elevThirdCumulantFunctionID , elevThirdCumulantFunction
		HDF_SD_ATTRSET, elevThirdCumulantFunctionID , 'long_name', $
			'Realisation averaged elevation Third Cumulant Function, ThirdCumulantFunction[0 ... NN-1][0 ... NN-1]'
		thirdCumulantFunctionDimID = HDF_SD_DIMGETID(elevThirdCumulantFunctionID, 0)
		HDF_SD_DIMSET, thirdCumulantFunctionDimID, LABEL='Third Cumulant Function Data length', NAME='NN', SCALE=length2, UNIT='meters'
		thirdCumulantFunctionDimID = HDF_SD_DIMGETID(elevThirdCumulantFunctionID, 1)
		HDF_SD_DIMSET, thirdCumulantFunctionDimID, LABEL='Third Cumulant Function Data length', NAME='NN', SCALE=length2, UNIT='meters'
		HDF_SD_ENDACCESS, elevThirdCumulantFunctionID 
		
		;;; Slope Third Cumulant Function
		
		slopeThirdCumulantFunctionID = HDF_SD_CREATE(fileID, "Average Slope Third Cumulant Function"    , [NN,NN], /FLOAT)
		HDF_SD_ADDDATA, slopeThirdCumulantFunctionID, slopeThirdCumulantFunction
		HDF_SD_ATTRSET, slopeThirdCumulantFunctionID, 'long_name', $
			'Realisation averaged slope Third Cumulant Function, ThirdCumulantFunction[0 ... NN-1][0 ... NN-1]'
		thirdCumulantFunctionDimID = HDF_SD_DIMGETID(slopeThirdCumulantFunctionID, 0)
		HDF_SD_DIMSET, thirdCumulantFunctionDimID, LABEL='Third Cumulant Function Data length', NAME='NN', SCALE=length2, UNIT='meters'
		thirdCumulantFunctionDimID = HDF_SD_DIMGETID(slopeThirdCumulantFunctionID, 1)
		HDF_SD_DIMSET, thirdCumulantFunctionDimID, LABEL='Third Cumulant Function Data length', NAME='NN', SCALE=length2, UNIT='meters'
		HDF_SD_ENDACCESS, slopeThirdCumulantFunctionID
		
		;;; Glint Third Cumulant Function
		
		glintThirdCumulantFunctionID = HDF_SD_CREATE(fileID, "Average Glint Third Cumulant Function"    , [NN,NN, GEOM.N_angles], /FLOAT)
		HDF_SD_ADDDATA, glintThirdCumulantFunctionID, glintThirdCumulantFunction
		HDF_SD_ATTRSET, glintThirdCumulantFunctionID, 'long_name', $
			'Realisation averaged glint Third Cumulant Function, ThirdCumulantFunction[0 ... NN-1][0 ... NN-1][0 ... N_angles-1]'
		thirdCumulantFunctionDimID = HDF_SD_DIMGETID(glintThirdCumulantFunctionID, 0)
		HDF_SD_DIMSET, thirdCumulantFunctionDimID, LABEL='Third Cumulant Function Data length', NAME='NN', SCALE=length2, UNIT='meters'
		thirdCumulantFunctionDimID = HDF_SD_DIMGETID(glintThirdCumulantFunctionID, 1)
		HDF_SD_DIMSET, thirdCumulantFunctionDimID, LABEL='Third Cumulant Function Data length', NAME='NN', SCALE=length2, UNIT='meters'
		numAnglesDimID = HDF_SD_DIMGETID(glintThirdCumulantFunctionID, 2)
		HDF_SD_DIMSET, numAnglesDimID, LABEL='Number of Angles', NAME='N_angles', UNIT='radians'
		HDF_SD_ENDACCESS, glintThirdCumulantFunctionID

	ENDIF ELSE BEGIN
		
		;;; Update the existing variables

		print  'Updating the elevation, slope and glint moments...'

		sds_index  = hdf_sd_nametoindex(fileID, 'Elevation Moments')
		sds_id = hdf_sd_select(fileID,sds_index)
		HDF_SD_ADDDATA, sds_id , elevMoments
		HDF_SD_ENDACCESS, sds_id 
		
		sds_index = hdf_sd_nametoindex(fileID, 'Slope Moments') 
		sds_id = hdf_sd_select(fileID,sds_index)
		HDF_SD_ADDDATA, sds_id, slopeMoments
		HDF_SD_ENDACCESS, sds_id
		
		sds_index = hdf_sd_nametoindex(fileID, 'Glint Moments')
		sds_id = hdf_sd_select(fileID,sds_index)
		HDF_SD_ADDDATA, sds_id, glintFirstMoments
		HDF_SD_ENDACCESS, sds_id

		;;;;;;;;;
		
		print  'Updating the free and bound indices...'
		
		sds_index  = hdf_sd_nametoindex(fileID, "Indices of first free wavenumber modes")
		sds_id = hdf_sd_select(fileID,sds_index)
		HDF_SD_ADDDATA, sds_id, NLCOUPLING.free1
		HDF_SD_ENDACCESS, sds_id 
		
		sds_index  = hdf_sd_nametoindex(fileID, "Indices of second free wavenumber modes")
		sds_id = hdf_sd_select(fileID,sds_index)
		HDF_SD_ADDDATA, sds_id, NLCOUPLING.free2
		HDF_SD_ENDACCESS, sds_id 

		sds_index  = hdf_sd_nametoindex(fileID, "Indices of bound wavenumber modes")
		sds_id = hdf_sd_select(fileID,sds_index)
		HDF_SD_ADDDATA, sds_id, NLCOUPLING.bound
		HDF_SD_ENDACCESS, sds_id 

		;;;;;;;;;

		print  'Updating the elevation, slope and glint power spectra...'

		sds_index = hdf_sd_nametoindex(fileID,"Average Elevation Power Spectrum")
		sds_id = hdf_sd_select(fileID,sds_index)
		HDF_SD_ADDDATA, sds_id, totalElevAvgPower
		HDF_SD_ENDACCESS, sds_id

		sds_index = hdf_sd_nametoindex(fileID,"Average Slope Power Spectrum")
		sds_id = hdf_sd_select(fileID,sds_index)
		HDF_SD_ADDDATA, sds_id, totalSlopeAvgPower
		HDF_SD_ENDACCESS, sds_id 

		sds_index = hdf_sd_nametoindex(fileID,"Average Glint Power Spectrum")
		sds_id = hdf_sd_select(fileID,sds_index)
		HDF_SD_ADDDATA, sds_id, totalGlintAvgPower
		HDF_SD_ENDACCESS, sds_id

		;;;;;;;;;

		print  'Updating the elevation, slope and glint bispectra...'

		sds_index = hdf_sd_nametoindex(fileID,"Average Elevation Bispectrum")
		sds_id = hdf_sd_select(fileID,sds_index)
		HDF_SD_ADDDATA, sds_id, elevBispectrum
		HDF_SD_ENDACCESS, sds_id

		sds_index = hdf_sd_nametoindex(fileID,"Average Slope Bispectrum")
		sds_id = hdf_sd_select(fileID,sds_index)
		HDF_SD_ADDDATA, sds_id, slopeBispectrum
		HDF_SD_ENDACCESS, sds_id
                                                                              
		sds_index = hdf_sd_nametoindex(fileID,"Average Glint Bispectrum")
		sds_id = hdf_sd_select(fileID,sds_index)
		HDF_SD_ADDDATA, sds_id, glintBispectrum
		HDF_SD_ENDACCESS, sds_id

		;;;;;;;;;

		print  'Updating the elevation, slope and glint third moment functions...'

		sds_index = hdf_sd_nametoindex(fileID,"Average Elevation Third Moment Function")
		sds_id = hdf_sd_select(fileID,sds_index)
		HDF_SD_ADDDATA, sds_id, elevThirdMomentFunction
		HDF_SD_ENDACCESS, sds_id
                                                                              
		sds_index = hdf_sd_nametoindex(fileID,"Average Slope Third Moment Function")
		sds_id = hdf_sd_select(fileID,sds_index)
		HDF_SD_ADDDATA, sds_id, slopeThirdMomentFunction
		HDF_SD_ENDACCESS, sds_id
                                                                              
		sds_index = hdf_sd_nametoindex(fileID,"Average Glint Third Moment Function")
		sds_id = hdf_sd_select(fileID,sds_index)
		HDF_SD_ADDDATA, sds_id, glintThirdMomentFunction
		HDF_SD_ENDACCESS, sds_id

	ENDELSE

	

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;;   Close the output HDF file   ;;;
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	HDF_SD_END, fileID
	print  'Write Operation Completed'
	print  ''






	
END

"""


###################################################
#                  Main Function                  #
###################################################

def main():

    spectrumChoices=['phillips_3','phillips_4','gaussian']

    description = \
    '''
    This is a brief description of %prog
    '''
    usage = "usage: %prog [mandatory args] [options]"
    version = __version__
    parser = optparse.OptionParser(description=description,usage=usage,version=version)

    # Mandatory arguments
    mandatoryGroup = optparse.OptionGroup(parser, "Mandatory Arguments",
                        "At a minimum these arguments must be specified")

    mandatoryGroup.add_option('-n','--numdatapoints',
                      action="store",
                      dest="N" ,
                      #default='1024',
                      type="int",
                      help="Number of points of the glint dataset. Must be an integer power of 2")
    mandatoryGroup.add_option('-N','--num2Dpoints',
                      action="store",
                      dest="NN" ,
                      #default='64',
                      type="int",
                      help="Number of points of the bispectrum array side. Must be an integer power of 2")
    mandatoryGroup.add_option('-d','--deltax',
                      action="store",
                      dest="delta_x",
                      default='0.2',
                      type="float",
                      help="The spatial increment in meters. [default: %default]")
    mandatoryGroup.add_option('-r','--num_realisations',
                      action="store",
                      dest="N_r" ,
                      #default='100',
                      type="int",
                      help="Number of realisations")
    mandatoryGroup.add_option('-S','--spectrum_type',
                      action="store",
                      dest="spectrumType",
                      type="string",
                      help='''Form of the elevation power spectrum.\n\n
                                                   Possible values are...
                                                   %s
                                                   ''' % (spectrumChoices.__str__()[1:-1]))

    parser.add_option_group(mandatoryGroup)

    # Optional arguments
    optionalGroup = optparse.OptionGroup(parser, "Extra Options",
                        "These options may be used to customize plot characteristics.")

    optionalGroup.add_option('-l','--nonlinear_modes',
                      action="store_true",
                      dest="nlSwitch",
                      help="Switch on nonlinear modes.")
    optionalGroup.add_option('-o','--output_file',
                      action="store",
                      dest="outputFile",
                      default="outGrid.h5",
                      type="string",
                      help="The full path of the output HDF5 file. [default: %default]")


    parser.add_option_group(optionalGroup)

    # Parse the arguments from the command line
    (options, args) = parser.parse_args()

    # Check that all of the mandatory options are given. If one or more 
    # are missing, print error message and exit...
    mandatories = ['N', 'NN','delta_x','N_r', 'spectrumType','nlSwitch']
    mand_errors = ["Missing mandatory argument [-n N            | --numdatapoints=N]",
                   "Missing mandatory argument [-N NN           | --num2Dpoints=NN]",
                   "Missing mandatory argument [-d delta_x      | --deltax=delta_x]",
                   "Missing mandatory argument [-r N_r          | --num_realisations=N_r]",
                   "Missing mandatory argument [-S spectrumType | --spectrum_type=spectrumType]"
                  ]
    isMissingMand = False
    for m,m_err in zip(mandatories,mand_errors):
        if not options.__dict__[m]:
            isMissingMand = True
            print m_err
    if isMissingMand :
        parser.error("Incomplete mandatory arguments, aborting...")

    N            = options.N
    NN           = options.NN
    delta_x      = options.delta_x
    N_r          = options.N_r
    spectrumType = options.spectrumType
    specExp      = 4 if options.spectrumType=='phillips_4' else 3
    nlSwitch     = options.nlSwitch

    cumulantFunctionSimulate(N,NN,delta_x,N_r,spectrumType,specExp,nlSwitch)

    sys.exit(0)

if __name__ == '__main__':
    main()


