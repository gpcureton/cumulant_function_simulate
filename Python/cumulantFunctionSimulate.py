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
from utility import bispectrumSymmetry, bicoherenceSymmetry, biCovarianceSymmetry


# Initialise the scale structure
class ScaleStruct :
    def __init__(self, N, NN, delta_x) :
        self.N = N
        self.delta_x = delta_x
        self.x_max = double(N-1)*delta_x
        self.x = np.arange(N)*delta_x

        self.k_max = 2.*pi/delta_x
        self.delta_k = self.k_max/double(N-1)
        self.k_N = double(N/2)*self.delta_k # Nyquist "wavenumber"
        self.k = np.arange(N)*self.delta_k

        self.N2 = N/2

        self.NN = NN
        self.NN2 = NN/2
        self.NN4 = NN/4

class GeomStruct :
    def __init__(self, N_geoms, angleLo, angleHi) :
        self.N_geoms = N_geoms
        self.source_angle = np.arange(N_geoms)
        self.detector_angle = np.arange(N_geoms)
        self.xi_min = np.arange(N_geoms)
        self.xi_0 = np.arange(N_geoms)
        self.xi_max = np.arange(N_geoms)

        self.angleRange=(angleHi-angleLo)

        if ((N_geoms-1) == 0) :
            self.d_angle = 0
        else :
            self.d_angle = self.angleRange/(N_geoms-1)

        self.start_angle = angleLo
        d2r = pi/180.
        r2d = 180./pi
        beta = 0.68*d2r

        self.source_angle = ((self.start_angle + np.arange(N_geoms,dtype=double)*self.d_angle))*d2r
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

    Scale = ScaleStruct(N, NN, delta_x)

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

    N2 = Scale.N2

    NN2 = Scale.NN2
    NN4 = Scale.NN4

    #   Populate the GEOM structure with angular quantities

    N_geoms = 5L
    angleLo=10.
    angleHi=30.
    Geom = GeomStruct(N_geoms,angleLo,angleHi)

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
    elevBispectrum = np.zeros((NN,NN),dtype=np.complex)
    elevComponentPower = np.zeros((NN,NN),dtype=np.float)
    elevSumPower = np.zeros((NN,NN),dtype=np.float)

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

    slopeBispectrum = np.zeros((NN,NN),dtype=np.complex)
    slopeComponentPower = np.zeros((NN,NN),dtype=np.float)
    slopeSumPower = np.zeros((NN,NN),dtype=np.float)

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

    curvatureBispectrum = np.zeros((NN,NN),dtype=np.complex)
    curvatureComponentPower = np.zeros((NN,NN),dtype=np.float)
    curvatureSumPower = np.zeros((NN,NN),dtype=np.float)


    #   Define the glint, glint spectrum and glint power
    glint               = np.zeros(N,dtype=double)
    glintSpectrum       = np.zeros(N,dtype=np.complex)
    totalGlintAvgPower  = np.zeros((N_geoms,N),dtype=double)
    glintBispectrum     = np.zeros((N_geoms,NN,NN),dtype=np.complex)
    glintComponentPower = np.zeros((N_geoms,NN,NN),dtype=np.float)
    glintSumPower       = np.zeros((N_geoms,NN,NN),dtype=np.float)

    #   Define the various point estimators for the elevation,
    #   slope and glint

    numMoments = 3

    elevStats = DataStatsStruct(numMoments)
    slopeStats = DataStatsStruct(numMoments)
    curvatureStats = DataStatsStruct(numMoments)
    glintStats = [DataStatsStruct(numMoments) for geoms in np.arange(N_geoms) ]

    #   Loop through the surface realisations for the quadratically
    #   coupled oscillations

    seed = 30
    N_r_cum = 0
    geom_runs = np.zeros(N_geoms,np.long)
    geom_runsCum = np.zeros(N_geoms,np.long)

    #time.sleep(3.)
    t1 = time.time()

    while (geom_runs.sum() < N_r * N_geoms) :

        N_r_cum += 1

        t2 = time.time()
        #print "Elapsed time = ",t2-t1
        if ((t2-t1) > 0.5):
            print "\n>>>>>>>>>>>>>>>>>>>>>\n"
            print "Computing realisation: %d at time %f" % (N_r_cum,(t2-t1))
            t1 = time.time()


        ### Compute the independent phases for this realisation
        primaryElevPhase = np.random.rand(N)*2.*pi - pi
        nlElevPhase      = np.random.rand(N)*2.*pi - pi

        ### Apply the phase correlations between the free and bound wavenumbers for the nonlinear
        ### component, if (nlSwitch==1)
        if (nlSwitch==1) :
            nlElevPhase[NLCoupling.bound] = primaryElevPhase[NLCoupling.free1] + \
                    primaryElevPhase[NLCoupling.free2]

        ##################################################################
        # Compute the elevation realisations from the elevation spectra  #
        # and the synthesised phases                                     #
        ##################################################################
        
        ### Calculate the elevation spectrum for the free waves
        primaryElevSpectrum = primaryElevAmplitude*(cos(primaryElevPhase) + 1j*sin(primaryElevPhase))
        primaryElevSpectrum[N/2+1 :] = np.conjugate(primaryElevSpectrum[1 : N/2][::-1])

        ### Calculate the elevation spectrum for the bound waves
        nlElevSpectrum = nlElevAmplitude*(cos(nlElevPhase) + 1j*sin(nlElevPhase))
        nlElevSpectrum[N/2+1 :] = np.conjugate(nlElevSpectrum[1 : N/2][::-1])

        ### Compute specific realisation of the free and bound waves. Nonlinear elevation
        ### (totalElevSurface) is sum of free and bound waves.
        primaryElevSurface = fft(primaryElevSpectrum)                ### Free waves
        nlElevSurface      = fft(nlElevSpectrum)                     ### Bound waves
        totalElevSurface   = primaryElevSurface + nlElevSurface      ### Total surface
 
        ### Compute the average power spectrum for free, bound and total elevation waves
        primaryElevAvgPower += abs(ifft(primaryElevSurface))**2.
        nlElevAvgPower      += abs(ifft(nlElevSurface))**2.
        totalElevAvgPower   += abs(ifft(totalElevSurface))**2.

        ### Compute the elevation estimators for this realisation
        elevStats.mean     += np.mean(totalElevSurface.real)
        elevStats.variance += np.var(totalElevSurface.real)
        elevStats.skewness += stats.skew(totalElevSurface.real)

        elevStats.moments += [ np.sum(totalElevSurface.real    )/double(N), \
                               np.sum(totalElevSurface.real**2.)/double(N), \
                               np.sum(totalElevSurface.real**3.)/double(N)]

        # Compute the Fourier spectrum of the total surface
        totalElevSpectrum = ifft(totalElevSurface)

        # Calculate the elevation bispectrum (reduced domain) for this realisation
        for j in np.arange(NN4+1):
            for i in np.arange(j,NN2-j+1):
                elevBispectrum[i,j] += totalElevSpectrum[i]*totalElevSpectrum[j] \
                        * np.conjugate(totalElevSpectrum[i+j])

        # Calculate the elevation component power spectrum (reduced domain) for this realisation
        for j in np.arange(NN4+1):
            for i in np.arange(j,NN2-j+1):
                elevComponentPower[i,j] += (abs(totalElevSpectrum[i]*totalElevSpectrum[j]))**2.

        # Calculate the elevation sum power spectrum (reduced domain) for this realisation
        for j in np.arange(NN4+1):
            for i in np.arange(j,NN2-j+1):
                elevSumPower[i,j] += (abs(totalElevSpectrum[i+j]))**2.

        #########################################################
        # Compute the slope realisations from the slope spectra #
        # and the synthesised phases                            #
        #########################################################
        
        ### Calculate the slope spectrum for the free waves
        primarySlopeSpectrum = primarySlopeAmplitude*(sin(primaryElevPhase) - 1j*cos(primaryElevPhase))
        primarySlopeSpectrum[N/2+1 :] = np.conjugate(primarySlopeSpectrum[1 : N/2][::-1])

        ### Calculate the slope spectrum for the bound waves
        nlSlopeSpectrum = nlSlopeAmplitude*(sin(nlElevPhase) - 1j*cos(nlElevPhase))
        nlSlopeSpectrum[N/2+1 :] = np.conjugate(nlSlopeSpectrum[1 : N/2][::-1])

        ### Compute specific realisation of the free and bound waves. Nonlinear slope
        ### (totalSlopeSurface) is sum of free and bound waves.
        primarySlopeSurface = fft(primarySlopeSpectrum)             ### Free waves
        nlSlopeSurface = fft(nlSlopeSpectrum)                       ### Bound waves
        totalSlopeSurface = primarySlopeSurface + nlSlopeSurface    ### Total surface

        ### Compute the average power spectrum for free, bound and total elevation waves
        primarySlopeAvgPower += abs(ifft(primarySlopeSurface))**2.
        nlSlopeAvgPower += abs(ifft(nlSlopeSurface))**2.
        totalSlopeAvgPower += abs(ifft(totalSlopeSurface))**2.

        ### Compute the slope estimators
        slopeStats.mean     += np.mean(totalSlopeSurface.real)
        slopeStats.variance += np.var(totalSlopeSurface.real)
        slopeStats.skewness += stats.skew(totalSlopeSurface.real)

        slopeStats.moments += [ np.sum(totalSlopeSurface    )/double(N), \
                                np.sum(totalSlopeSurface**2.)/double(N), \
                                np.sum(totalSlopeSurface**3.)/double(N) ]

        # Compute the Fourier spectrum of the total surface
        totalSlopeSpectrum = ifft(totalSlopeSurface)

        # Calculate the slope bispectrum (reduced domain) for this realisation
        for j in np.arange(NN4+1):
            for i in np.arange(j,NN2-j+1):
                slopeBispectrum[i,j] += totalSlopeSpectrum[i]*totalSlopeSpectrum[j] \
                        * np.conjugate(totalSlopeSpectrum[i+j])

        # Calculate the slope component power spectrum (reduced domain) for this realisation
        for j in np.arange(NN4+1):
            for i in np.arange(j,NN2-j+1):
                slopeComponentPower[i,j] += (abs(totalSlopeSpectrum[i]*totalSlopeSpectrum[j]))**2.

        # Calculate the slope sum power spectrum (reduced domain) for this realisation
        for j in np.arange(NN4+1):
            for i in np.arange(j,NN2-j+1):
                slopeSumPower[i,j] += (abs(totalSlopeSpectrum[i+j]))**2.

        #################################################################
        # Compute the curvature realisations from the curvature spectra #
        # and the synthesised phases                                    #
        #################################################################
        
        ### Calculate the curvature spectrum for the free waves
        primaryCurvatureSpectrum = primaryCurvatureAmplitude*(-cos(primaryElevPhase) - 1j*sin(primaryElevPhase))
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

        # Compute the curvature estimators
        curvatureStats.mean     += np.mean(totalCurvatureSurface.real)
        curvatureStats.variance += np.var(totalCurvatureSurface.real)
        curvatureStats.skewness += stats.skew(totalCurvatureSurface.real)

        curvatureStats.moments += [ np.sum(totalCurvatureSurface    )/double(N), \
                                    np.sum(totalCurvatureSurface**2.)/double(N), \
                                    np.sum(totalCurvatureSurface**3.)/double(N) ]

        # Compute the Fourier spectrum of the total surface
        totalCurvatureSpectrum = ifft(totalCurvatureSurface)

        # Calculate the curvature bispectrum (reduced domain) for this realisation
        for j in np.arange(NN4+1):
            for i in np.arange(j,NN2-j+1):
                curvatureBispectrum[i,j] += totalCurvatureSpectrum[i]*totalCurvatureSpectrum[j] \
                        * np.conjugate(totalCurvatureSpectrum[i+j])

        # Calculate the curvature component power spectrum (reduced domain) for this realisation
        for j in np.arange(NN4+1):
            for i in np.arange(j,NN2-j+1):
                curvatureComponentPower[i,j] += (abs(totalCurvatureSpectrum[i]*totalCurvatureSpectrum[j]))**2.

        # Calculate the curvature sum power spectrum (reduced domain) for this realisation
        for j in np.arange(NN4+1):
            for i in np.arange(j,NN2-j+1):
                curvatureSumPower[i,j] += (abs(totalCurvatureSpectrum[i+j]))**2.



        #####################################################
        # Loop through the geometries in the GEOM structure #
        #####################################################

        for geom in np.arange(N_geoms) :

            # Check if we have finished processing for this
            # geom.

            if (geom_runs[geom] < N_r) :

                #print "\n\tProcessing geom ",geom," for run ",geom_runs[geom]+1, \
                #" --> attempt ",geom_runsCum[geom]+1

                #   Compute the glint realisation from the slope
                #   realisations

                slopeMin = Geom.xi_min[geom]
                slopeMax = Geom.xi_max[geom]

                glint = np.double(totalSlopeSurface.real > slopeMin) * np.double(totalSlopeSurface.real < slopeMax)

                ### Check if all glint elements vanish
                result = np.where(glint)

                if (np.shape(np.squeeze(np.where(glint))) == (0,)) :

                    #print "\tZero-glint realisation geom ",geom, \
                    #" for run ",geom_runs[geom]+1, \
                    #" --> attempt ",geom_runsCum[geom]+1

                    ### There are no glints, add to attempts count
                    ### If this geom fails and there are no glints, then steeper 
                    ### geoms will fail also, so break out of the geom loop and 
                    ### proceed to the next realisation...

                    geom_runsCum[geom:N_geoms] += 1
                    break

                else :

                    #print "\tSuccessful realisation geom ",geom, \
                    #" for run ",geom_runs[geom] + 1, \
                    #" --> attempt ",geom_runsCum[geom]+1

                    geom_runs[geom] += 1
                    geom_runsCum[geom] += 1

                    ### Compute the glint moments

                    glintStats[geom].mean     += np.mean(glint.real)
                    glintStats[geom].variance += np.var( glint.real)
                    glintStats[geom].skewness += stats.skew(glint.real)

                    glintStats[geom].moments += [ np.sum(glint.real    )/double(N), \
                                                   np.sum(glint.real**2.)/double(N), \
                                                   np.sum(glint.real**3.)/double(N) ]

                    ### Compute the Fourier spectrum of this glint realisation (using ifft 
                    ### which has same normalisation as IDL FFT routine).
                    glintSpectrum = ifft(glint)

                    # Compute the average glint power spectrum
                    # TODO: incorporate windowing to minimise aliasing
                    totalGlintAvgPower[geom] += abs(glintSpectrum)**2.

                    # Calculate the glint bispectrum (reduced domain) for this realisation
                    for j in np.arange(NN4+1):
                        for i in np.arange(j,NN2-j+1):
                            glintBispectrum[geom,i,j] += glintSpectrum[i]*glintSpectrum[j] \
                                    * np.conjugate(glintSpectrum[i+j])

                    # Calculate the glint component power spectrum (reduced domain) for this realisation
                    for j in np.arange(NN4+1):
                        for i in np.arange(j,NN2-j+1):
                            glintComponentPower[geom,i,j] += (abs(glintSpectrum[i]*glintSpectrum[j]))**2.

                    # Calculate the glint sum power spectrum (reduced domain) for this realisation
                    for j in np.arange(NN4+1):
                        for i in np.arange(j,NN2-j+1):
                            glintSumPower[geom,i,j] += (abs(glintSpectrum[i+j]))**2.


                ### End checking for zero-glint of this geom

            ### End checking for completion of this geom

        ### End geom loop

    ### End realisation loop
    
    print ""
    print "geom_runs:    ",geom_runs," ... for total of ", \
        int(np.sum(geom_runs))," / ",N_r * N_geoms
    print "geom_runsCum: ",geom_runsCum
    
    N_runs = N_r
    N_r = N_r_cum
    print  "Final N_r_cum = ",N_r_cum

    #    Have a look at the scipy calculated stats

    ########################################
    #   Compute the elevation estimators   #
    ########################################

    # Normalise the estimators
    elevStats.mean     /= double(N_runs)
    elevStats.variance /= double(N_runs)
    elevStats.skewness /= double(N_runs)
    elevStats.moments /= double(N_runs)

    # Compute the cumulants
    elevStats.cumulantsFromMoments()

    # Compute the average elevation moment and cumulant functions.
    elevSecondMomentFunction =  fft(totalElevAvgPower)
    elevSecondMomentFunction /= elevSecondMomentFunction.real[0]

    # Compute the second order elev cumulant function
    elevSecondCumulantFunction = (elevStats.moments[1]*elevSecondMomentFunction - 
            elevStats.moments[0]**2.)/elevStats.cumulants[1]

	# Compute the bispectrum estimators
    elevBispectrum /= float(N_r)
    elevComponentPower /= float(N_r)
    elevSumPower /= float(N_r)

	# Compute the bicoherence
    elevBicoherence = np.zeros((NN,NN),dtype=np.float)
    for j in np.arange(NN4+1):
        for i in np.arange(j,NN2-j+1):
            if (sqrt(elevComponentPower[i,j])*sqrt(elevSumPower[i,j]) > 10.**(-12.)):
				elevBicoherence[i,j] = abs(elevBispectrum[i,j])/ \
                        (sqrt(elevComponentPower[i,j])*sqrt(elevSumPower[i,j]))
            else:
				elevBicoherence[i,j] = 0.

	# Fill the rest of the bispectrum and bicoherence array
	bispectrumSymmetry(elevBispectrum,NN)
	bicoherenceSymmetry(elevBicoherence,NN)

	# Compute the elevation third moment function
    elevThirdMomentFunction = np.zeros((NN,NN),dtype=np.float)
    elevThirdMomentFunction =  fft(elevBispectrum).real
    elevThirdMomentFunction /= elevThirdMomentFunction[0,0]

	# Compute the elevation third cumulant function
    elevThirdCumulantFunction = np.zeros((NN,NN),dtype=np.float)
    for i in np.arange(0,NN/2+1):
        for j in np.arange(0,i+1):
			elevThirdCumulantFunction[i,j] = (elevStats.moments[2]*elevThirdMomentFunction[i,j] \
				- elevStats.moments[0]*elevStats.moments[1] * \
				(elevSecondMomentFunction[i] + \
				 elevSecondMomentFunction[j] + \
				 elevSecondMomentFunction[abs(j-i)]) \
				+ 2.*elevStats.moments[0]**3.)/elevStats.cumulants[2]

			elevThirdCumulantFunction[j,i] = elevThirdCumulantFunction[i,j]


    ########################################
    #     Compute the slope estimators     #
    ########################################

    # Normalise the estimators
    slopeStats.mean     /= double(N_runs)
    slopeStats.variance /= double(N_runs)
    slopeStats.skewness /= double(N_runs)
    slopeStats.moments /= double(N_runs)

    # Compute the cumulants
    slopeStats.cumulantsFromMoments()

    # compute the average slope moment and cumulant functions.
    slopeSecondMomentFunction =  fft(totalSlopeAvgPower)
    slopeSecondMomentFunction /= slopeSecondMomentFunction.real[0]

	# Compute the bispectrum estimators
    slopeBispectrum /= float(N_r)
    slopeComponentPower /= float(N_r)
    slopeSumPower /= float(N_r)

	# Compute the bicoherence
    slopeBicoherence = np.zeros((NN,NN),dtype=np.float)
    for j in np.arange(NN4+1):
        for i in np.arange(j,NN2-j+1):
            if (sqrt(slopeComponentPower[i,j])*sqrt(slopeSumPower[i,j]) > 10.**(-12.)):
				slopeBicoherence[i,j] = abs(slopeBispectrum[i,j])/ \
                        (sqrt(slopeComponentPower[i,j])*sqrt(slopeSumPower[i,j]))
            else:
				slopeBicoherence[i,j] = 0.

	# Fill the rest of the bispectrum and bicoherence array
	bispectrumSymmetry(slopeBispectrum,NN)
	bicoherenceSymmetry(slopeBicoherence,NN)

	# Compute the slope third moment function
    slopeThirdMomentFunction = np.zeros((NN,NN),dtype=np.float)
    slopeThirdMomentFunction =  fft(slopeBispectrum).real
    slopeThirdMomentFunction /= slopeThirdMomentFunction[0,0]

	# Compute the slope third cumulant function
    slopeThirdCumulantFunction = np.zeros((NN,NN),dtype=np.float)
    for i in np.arange(0,NN/2+1):
        for j in np.arange(0,i+1):
			slopeThirdCumulantFunction[i,j] = (slopeStats.moments[2]*slopeThirdMomentFunction[i,j] \
				- slopeStats.moments[0]*slopeStats.moments[1] * \
				(slopeSecondMomentFunction[i] + \
				 slopeSecondMomentFunction[j] + \
				 slopeSecondMomentFunction[abs(j-i)]) \
				+ 2.*slopeStats.moments[0]**3.)/slopeStats.cumulants[2]

			slopeThirdCumulantFunction[j,i] = slopeThirdCumulantFunction[i,j]


    ########################################
    #   Compute the curvature estimators   #
    ########################################

    # Normalise the estimators
    curvatureStats.mean     /= double(N_runs)
    curvatureStats.variance /= double(N_runs)
    curvatureStats.skewness /= double(N_runs)
    curvatureStats.moments /= double(N_runs)

    # Compute the cumulants
    curvatureStats.cumulantsFromMoments()

    # compute the average curvature moment and cumulant functions.
    curvatureSecondMomentFunction =  fft(totalCurvatureAvgPower)
    curvatureSecondMomentFunction /= curvatureSecondMomentFunction.real[0]

	# Compute the bispectrum estimators
    curvatureBispectrum /= float(N_r)
    curvatureComponentPower /= float(N_r)
    curvatureSumPower /= float(N_r)

	# Compute the bicoherence
    curvatureBicoherence = np.zeros((NN,NN),dtype=np.float)
    for j in np.arange(NN4+1):
        for i in np.arange(j,NN2-j+1):
            if (sqrt(curvatureComponentPower[i,j])*sqrt(curvatureSumPower[i,j]) > 10.**(-12.)):
				curvatureBicoherence[i,j] = abs(curvatureBispectrum[i,j])/ \
                        (sqrt(curvatureComponentPower[i,j])*sqrt(curvatureSumPower[i,j]))
            else:
				curvatureBicoherence[i,j] = 0.

	# Fill the rest of the bispectrum and bicoherence array
	bispectrumSymmetry(curvatureBispectrum,NN)
	bicoherenceSymmetry(curvatureBicoherence,NN)

	# Compute the curvature third moment function
    curvatureThirdMomentFunction = np.zeros((NN,NN),dtype=np.float)
    curvatureThirdMomentFunction =  fft(curvatureBispectrum).real
    curvatureThirdMomentFunction /= curvatureThirdMomentFunction[0,0]

	# Compute the curvature third cumulant function
    curvatureThirdCumulantFunction = np.zeros((NN,NN),dtype=np.float)
    for i in np.arange(0,NN/2+1):
        for j in np.arange(0,i+1):
			curvatureThirdCumulantFunction[i,j] = (curvatureStats.moments[2]*curvatureThirdMomentFunction[i,j] \
				- curvatureStats.moments[0]*curvatureStats.moments[1] * \
				(curvatureSecondMomentFunction[i] + \
				 curvatureSecondMomentFunction[j] + \
				 curvatureSecondMomentFunction[abs(j-i)]) \
				+ 2.*curvatureStats.moments[0]**3.)/curvatureStats.cumulants[2]

			curvatureThirdCumulantFunction[j,i] = curvatureThirdCumulantFunction[i,j]


    ########################################
    #     Compute the glint estimators     #
    ########################################

    # Normalise the estimators
    for geom in np.arange(N_geoms) :
        glintStats[geom].mean     /= double(geom_runsCum[geom])
        glintStats[geom].variance /= double(geom_runsCum[geom])
        glintStats[geom].skewness /= double(geom_runsCum[geom])
        glintStats[geom].moments  /= double(geom_runsCum[geom])

    # Compute the cumulants
    for geom in np.arange(N_geoms) :
        glintStats[geom].cumulantsFromMoments()

    # compute the average glint moment and cumulant functions.
    glintSecondMomentFunction = np.zeros((N_geoms,N),dtype=double)
    for geom in np.arange(N_geoms) :
        glintSecondMomentFunction[geom] = fft(totalGlintAvgPower[geom]).real
        glintSecondMomentFunction[geom] /= glintSecondMomentFunction[geom][0]

	# Compute the bispectrum estimators
    glintBispectrum /= float(N_r)
    glintComponentPower /= float(N_r)
    glintSumPower /= float(N_r)

	# Compute the bicoherence
    glintBicoherence = np.zeros((N_geoms,NN,NN),dtype=np.float)
    for geom in np.arange(N_geoms) :
        for j in np.arange(NN4+1):
            for i in np.arange(j,NN2-j+1):
                if (sqrt(glintComponentPower[geom,i,j])*sqrt(glintSumPower[geom,i,j]) > 10.**(-12.)):
                    glintBicoherence[geom,i,j] = \
                            abs(glintBispectrum[geom,i,j]) / (sqrt(glintComponentPower[geom,i,j]) * sqrt(glintSumPower[geom,i,j]))
                else:
                    glintBicoherence[geom,i,j] = 0.

	# Fill the rest of the bispectrum and bicoherence array
    for geom in np.arange(N_geoms) :
        bispectrumSymmetry(glintBispectrum[geom],NN)
        bicoherenceSymmetry(glintBicoherence[geom],NN)

	# Compute the glint third moment function
    glintThirdMomentFunction = np.zeros((N_geoms,NN,NN),dtype=np.float)
    for geom in np.arange(N_geoms) :
        glintThirdMomentFunction[geom] =  fft(glintBispectrum[geom]).real
        glintThirdMomentFunction[geom] /= glintThirdMomentFunction[geom,0,0]

	# Compute the glint third cumulant function
    glintThirdCumulantFunction = np.zeros((N_geoms,NN,NN),dtype=np.float)
    for geom in np.arange(N_geoms) :
        for i in np.arange(0,NN/2+1):
            for j in np.arange(0,i+1):
                glintThirdCumulantFunction[geom,i,j] = (glintStats[geom].moments[2]*glintThirdMomentFunction[geom,i,j] \
                    - glintStats[geom].moments[0]*glintStats[geom].moments[1] * \
                    (glintSecondMomentFunction[geom,i] + \
                     glintSecondMomentFunction[geom,j] + \
                     glintSecondMomentFunction[geom,abs(j-i)]) \
                    + 2.*glintStats[geom].moments[0]**3.)/glintStats[geom].cumulants[2]

                glintThirdCumulantFunction[geom,j,i] = glintThirdCumulantFunction[geom,i,j]


    #################################################
    ###  The elevation and slope summary results  ###
    #################################################

    #strFormat = "{:30}: {:10.6e}"
    strFormat = "{:30}: {:15.12f}"

    print '\n####################################'
    print '            Elevation'
    print '####################################'

    print ""
    print strFormat.format("Python Elevation mean",elevStats.mean)
    #print strFormat.format("Python Elevation stdev",sqrt(elevStats.variance))
    print strFormat.format("Python Elevation variance",elevStats.variance)
    print strFormat.format("Python Elevation skewness",elevStats.skewness)

    print ""
    print strFormat.format("Elevation first moment",elevStats.moments[0])
    print strFormat.format("Elevation second moment",elevStats.moments[1])
    print strFormat.format("Elevation third moment",elevStats.moments[2])

    print ""
    print strFormat.format("Elevation first cumulant",elevStats.cumulants[0])
    print strFormat.format("Elevation second cumulant",elevStats.cumulants[1])
    print strFormat.format("Elevation third cumulant",elevStats.cumulants[2])

    print ""
    print strFormat.format("Elevation mean",elevStats.mean)
    #print strFormat.format("Elevation stdev",sqrt(elevStats.variance))
    print strFormat.format("Elevation variance",elevStats.variance)
    print strFormat.format("Elevation skewness",elevStats.skewness)

    print '\n####################################'
    print '               Slope'
    print '####################################'

    print ""
    print strFormat.format("Python Slope mean",slopeStats.mean)
    #print strFormat.format("Python Slope stdev",sqrt(slopeStats.variance))
    print strFormat.format("Python Slope variance",slopeStats.variance)
    print strFormat.format("Python Slope skewness",slopeStats.skewness)

    print ""
    print strFormat.format("Slope first moment",slopeStats.moments[0])
    print strFormat.format("Slope second moment",slopeStats.moments[1])
    print strFormat.format("Slope third moment",slopeStats.moments[2])

    print ""
    print strFormat.format("Slope first cumulant",slopeStats.cumulants[0])
    print strFormat.format("Slope second cumulant",slopeStats.cumulants[1])
    print strFormat.format("Slope third cumulant",slopeStats.cumulants[2])

    print ""
    print strFormat.format("Slope mean",slopeStats.mean)
    #print strFormat.format("Slope stdev",sqrt(slopeStats.variance))
    print strFormat.format("Slope variance",slopeStats.variance)
    print strFormat.format("Slope skewness",slopeStats.skewness)

    print '\n####################################'
    print '              Curvature'
    print '####################################'

    print ""
    print strFormat.format("Python Curvature mean",curvatureStats.mean)
    #print strFormat.format("Python Curvature stdev",sqrt(curvatureStats.variance))
    print strFormat.format("Python Curvature variance",curvatureStats.variance)
    print strFormat.format("Python Curvature skewness",curvatureStats.skewness)

    print ""
    print strFormat.format("Curvature first moment",curvatureStats.moments[0])
    print strFormat.format("Curvature second moment",curvatureStats.moments[1])
    print strFormat.format("Curvature third moment",curvatureStats.moments[2])

    print ""
    print strFormat.format("Curvature first cumulant",curvatureStats.cumulants[0])
    print strFormat.format("Curvature second cumulant",curvatureStats.cumulants[1])
    print strFormat.format("Curvature third cumulant",curvatureStats.cumulants[2])

    print ""
    print strFormat.format("Curvature mean",curvatureStats.mean)
    #print strFormat.format("Curvature stdev",sqrt(curvatureStats.variance))
    print strFormat.format("Curvature variance",curvatureStats.variance)
    print strFormat.format("Curvature skewness",curvatureStats.skewness)

    print '\n####################################'
    print '              Glint'
    print '####################################'

    print "\nPython Glint mean, variance and skewness ...\n"
    for geom in np.arange(Geom.N_geoms) :
        print "\tGeometry %1d:\t\t%10.6e\t%10.6e\t%10.6e" % (geom,\
            glintStats[geom].mean,\
            glintStats[geom].variance,\
            glintStats[geom].skewness)

    print "Glint moments ...\n"
    for geom in np.arange(Geom.N_geoms) :
        print "\tGeometry %1d:\t\t%10.6e\t%10.6e\t%10.6e" % (geom,\
            glintStats[geom].moments[0],\
            glintStats[geom].moments[1],\
            glintStats[geom].moments[2])

    print "\nGlint cumulants ...\n"
    for geom in np.arange(Geom.N_geoms) :
        print "\tGeometry %1d:\t\t%10.6e\t%10.6e\t%10.6e" % (geom,\
            glintStats[geom].cumulants[0],\
            glintStats[geom].cumulants[1],\
            glintStats[geom].cumulants[2])

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


