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

import os, sys, logging, traceback
from os import path,uname,environ

import string, copy
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

# This must come *after* the backend is specified.
#import matplotlib.pyplot as ppl

import tables as pytables
from tables import exceptions as pyEx

# every module should have a LOG object
sourcename= file_Id.split(" ")
#LOG = logging.getLogger(sourcename[1])
LOG = logging.getLogger(__file__)

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

    # Determine the various scale parameters and populate the 
    # Scale class 

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

    #   Populate the GEOM structure with angular quantities

    N_angles = 5L
    angleLo=10.
    angleHi=30.
    Geom = GeomStruct(N_angles,angleLo,angleHi)

    #   Populate the elevation power spectrum structures ElevPower, 
    #   and NLCoupling                                          

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

    # Initialise the slope power spectrum structure
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

    # Initialise the curvature power spectrum structure
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

    #   Compute the total elevation amplitude, phase and spectrum,
    #   and the second moment function 

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

    #    Compute the total slope amplitude, phase and spectrum

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

    #    Compute the total curvature amplitude, phase and spectrum

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

    #   Define the glint, glint spectrum and glint power

    glint = np.zeros(N,dtype=double)
    glintSpectrum = np.zeros(N,dtype=np.complex64)
    totalGlintAvgPower = np.zeros((Geom.N_angles,N),dtype=double) # DBLARR(N,GEOM.N_angles)

    #   Define the various point estimators for the elevation,
    #   slope and glint

    numMoments = 3

    elevStats = DataStatsStruct(numMoments)
    slopeStats = DataStatsStruct(numMoments)
    curvatureStats = DataStatsStruct(numMoments)
    glintStats = [DataStatsStruct(numMoments) for geoms in np.arange(Geom.N_angles) ]

    #   Loop through the surface realisations for the quadratically
    #   coupled oscillations

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

        #   Compute the elevation realisations from the elevation spectra
        #   and the synthesised phases
        
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

        #   Compute the slope realisations from the slope spectra
        #   and the synthesised phases
        
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

        #   Compute the curvature realisations from the curvature spectra
        #   and the synthesised phases
        
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

        # Compute the curvature moments

        curvatureStats.mean     += np.mean(totalCurvatureSurface.real)
        curvatureStats.variance += np.var(totalCurvatureSurface.real)
        curvatureStats.skewness += stats.skew(totalCurvatureSurface.real)

        curvatureStats.moments += [ np.sum(totalCurvatureSurface    )/double(N), \
                                    np.sum(totalCurvatureSurface**2.)/double(N), \
                                    np.sum(totalCurvatureSurface**3.)/double(N) ]

        # Compute the Fourier spectrum of the total surface
        totalCurvatureSpectrum = ifft(totalCurvatureSurface)

        # Loop through the geometries in the GEOM structure

        for angle in np.arange(Geom.N_angles) :

            # Check if we have finished processing for this
            # angle.

            if (angleRuns[angle] < N_r) :

                #print "\n\tProcessing angle ",angle," for run ",angleRuns[angle]+1, \
                #" --> attempt ",angleRunsCum[angle]+1

                #   Compute the glint realisation from the slope
                #   realisations

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

                    # Compute the average glint power spectrum

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

    #    Have a look at the scipy calculated stats

    ########################################
    #   Compute the elevation estimators   #
    ########################################

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

    #   Compute the moments and cumulants our way
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

    # compute the average elevation moment and cumulant functions.
    elevSecondMomentFunction =  fft(totalElevAvgPower)
    elevSecondMomentFunction /= elevSecondMomentFunction.real[0]
    #elevSecondCumulantFunction = (elevMoments[1]*elevSecondMomentFunction - elevMoments[0]^2.D)/elevCumulants[1]

    # compute the average slope moment and cumulant functions.
    slopeSecondMomentFunction =  fft(totalSlopeAvgPower)
    slopeSecondMomentFunction /= slopeSecondMomentFunction.real[0]

    # compute the average glint moment and cumulant functions.
    glintSecondMomentFunction = np.zeros((Geom.N_angles,N),dtype=double)
    for geoms in np.arange(Geom.N_angles) :
        glintSecondMomentFunction[geoms] = fft(totalGlintAvgPower[geoms]).real
        glintSecondMomentFunction[geoms] /= glintSecondMomentFunction[geoms][0]


def _argparse():

    import argparse as argparse

    spectrumChoices=['phillips_3','phillips_4','gaussian']

    defaults = {'N':1024,
                'NN':64,
                'delta_x':0.02,
                'N_r':100,
                'spectrumType':'phillips_3',
                'nlSwitch':False,
                'outputFile':"outSimulatedGlint.h5",
                }

    description = "This is a brief description of %prog"
    usage = "usage: %prog [mandatory args] [optional args]"
    version = __version__

    parser = argparse.ArgumentParser()

    # Positional arguments...
    #parser.add_argument("x", type=int, help="the base")
    #parser.add_argument("y", type=int, help="the exponent")

    # Mandatory arguments...
    parser.add_argument('-n','--numdatapoints',
                      action="store",
                      dest="N" ,
                      default=defaults['N'],
                      type=int,
                      help="Number of points of the glint dataset. Must be an integer power of 2. [default: {}]".format(defaults['N']),
                      metavar="N")

    parser.add_argument('-N','--num2Dpoints',
                      action="store",
                      dest="NN" ,
                      default=defaults['NN'],
                      type=int,
                      help="Number of points of the bispectrum array side. Must be an integer power of 2. [default: {}]".format(defaults['NN']),
                      metavar="NN")

    parser.add_argument('-d','--deltax',
                      action="store",
                      dest="delta_x",
                      default=defaults['delta_x'],
                      type=float,
                      help="The spatial increment in meters. [default: {}]".format(defaults['delta_x']),
                      metavar="DELTA_X")

    parser.add_argument('-r','--num_realisations',
                      action="store",
                      dest="N_r" ,
                      default=defaults['N_r'],
                      type=int,
                      help="Number of realisations. [default: {}]".format(defaults['N_r']),
                      metavar="NUMREALS")

    parser.add_argument('-S','--spectrum_type',
                      action="store",
                      dest="spectrumType",
                      default=defaults['spectrumType'],
                      type=str,
                      help='''Form of the elevation power spectrum.\n\n
                                                   Possible values are...
                                                   {}. [default: '{}']
                                                   '''.format(spectrumChoices.__str__()[1:-1],defaults['spectrumType']),
                      metavar="SPECTRUMTYPE")

    # Optional arguments
    parser.add_argument('-l','--nonlinear_modes',
                      action="store_true",
                      dest="nlSwitch",
                      default=defaults['nlSwitch'],
                      help="Switch on nonlinear modes. [default: {}]".format(defaults['nlSwitch']))

    parser.add_argument('-o','--output_file',
                      action="store",
                      dest="outputFile",
                      default=defaults['outputFile'],
                      type=str,
                      help="The full path of the output HDF5 file. [default: '{}']".format(defaults['outputFile']),
                      metavar="OUTFILE")

    parser.add_argument("-v", "--verbose",
                      dest='verbosity',
                      action="count", 
                      default=0,
                      help='each occurrence increases verbosity 1 level from ERROR: -v=WARNING -vv=INFO -vvv=DEBUG')

    args = parser.parse_args()

    print args

    # Check that all of the mandatory args are given. If one or more 
    # are missing, print error message and exit...
    #mandatories = ['N', 'NN','delta_x','N_r', 'spectrumType','nlSwitch']
    #mand_errors = ["Missing mandatory argument [-n N            | --numdatapoints=N]",
                   #"Missing mandatory argument [-N NN           | --num2Dpoints=NN]",
                   #"Missing mandatory argument [-d delta_x      | --deltax=delta_x]",
                   #"Missing mandatory argument [-r N_r          | --num_realisations=N_r]",
                   #"Missing mandatory argument [-S spectrumType | --spectrum_type=spectrumType]"
                  #]
    #isMissingMand = False
    #for m,m_err in zip(mandatories,mand_errors):
        #if not args.__dict__[m]:
            #isMissingMand = True
            #print m_err
    #if isMissingMand :
        #parser.error("Incomplete mandatory arguments, aborting...")


    # Set up the logging
    #console_logFormat = '%(asctime)s : %(name)-12s: %(levelname)-8s %(message)s'
    #console_logFormat = '%(levelname)s:%(name)s:%(msg)s') # [%(filename)s:%(lineno)d]'
    console_logFormat = '%(asctime)s : (%(levelname)s):%(filename)s:%(funcName)s:%(lineno)d:  %(message)s'
    dateFormat = '%Y-%m-%d %H:%M:%S'
    levels = [logging.ERROR, logging.WARN, logging.INFO, logging.DEBUG]
    logging.basicConfig(level = levels[args.verbosity], 
            format = console_logFormat, 
            datefmt = dateFormat)

    return args


###################################################
#                  Main Function                  #
###################################################

def main():

    args = _argparse()

    N            = args.N
    NN           = args.NN
    delta_x      = args.delta_x
    N_r          = args.N_r
    spectrumType = args.spectrumType
    specExp      = 4 if args.spectrumType=='phillips_4' else 3
    nlSwitch     = args.nlSwitch

    cumulantFunctionSimulate(N,NN,delta_x,N_r,spectrumType,specExp,nlSwitch)

    LOG.info("Exiting...")
    sys.exit(0)


if __name__ == '__main__':
    main()


