#!/usr/bin/env python
# encoding: utf-8
"""
utility.py

Created by Geoff Cureton <geoff.cureton@ssec.wisc.edu> on 2011-09-30.
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

__author__ = 'Geoff Cureton <geoff.cureton@ssec.wisc.edu>'
__version__ = '$Id$'
__docformat__ = 'Epytext'


import numpy as np

"""
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Routine to use the symmetry properties of the bispectrum of a real 
;;; sequence to populate an entire NxN array from the primary octant. Takes
;;; as input an NxN complex array, and the array size N, and returns the 
;;; fully populated array in the input array
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO bispectrumSymmetry,bispectrum,N

    FOR j=0L,N/4L DO BEGIN
        FOR i=j,(N/2L-j) DO BEGIN
            bispectrum[(j LT 0L) ? N+(j) : j , (i LT 0L) ? N+(i) : i] = bispectrum[i,j]
            bispectrum[(j LT 0L) ? N+(j) : j , (-i-j LT 0L) ? N+(-i-j) : -i-j] = bispectrum[i,j]
            bispectrum[(-i-j LT 0L) ? N+(-i-j) : -i-j , (j LT 0L) ? N+(j) : j] = bispectrum[i,j]
            bispectrum[(-i-j LT 0L) ? N+(-i-j) : -i-j , (i LT 0L) ? N+(i) : i] = bispectrum[i,j]
            bispectrum[(i LT 0L) ? N+(i) : i , (-i-j LT 0L) ? N+(-i-j) : -i-j] = bispectrum[i,j]

            bispectrum[(-i LT 0L) ? N+(-i) : -i , (-j LT 0L) ? N+(-j) : -j   ] = CONJ(bispectrum[i,j])
            bispectrum[(-j LT 0L) ? N+(-j) : -j , (-i LT 0L) ? N+(-i) : -i   ] = CONJ(bispectrum[i,j])
            bispectrum[(-j LT 0L) ? N+(-j) : -j , (i+j LT 0L) ? N+(i+j) : i+j] = CONJ(bispectrum[i,j])
            bispectrum[(i+j LT 0L) ? N+(i+j) : i+j , (-j LT 0L) ? N+(-j) : -j] = CONJ(bispectrum[i,j])
            bispectrum[(i+j LT 0L) ? N+(i+j) : i+j , (-i LT 0L) ? N+(-i) : -i] = CONJ(bispectrum[i,j])
            bispectrum[(-i LT 0L) ? N+(-i) : -i , (i+j LT 0L) ? N+(i+j) : i+j] = CONJ(bispectrum[i,j])
        ENDFOR
    ENDFOR
END
"""
def bispectrumSymmetry(bispectrum,N):
    for j in range(N/4+1):
        for i in range(N/2-j+1):
            try :

                bispectrum[ N+(j)    if (j < 0L)    else  j    , N+(i)    if  (i < 0L)    else  i    ]  = bispectrum[i,j]
                bispectrum[ N+(j)    if (j < 0L)    else  j    , N+(-i-j) if  (-i-j < 0L) else  -i-j ]  = bispectrum[i,j]
                bispectrum[ N+(-i-j) if (-i-j < 0L) else  -i-j , N+(j)    if  (j < 0L)    else  j    ]  = bispectrum[i,j]
                bispectrum[ N+(-i-j) if (-i-j < 0L) else  -i-j , N+(i)    if  (i < 0L)    else  i    ]  = bispectrum[i,j]
                bispectrum[ N+(i)    if (i < 0L)    else  i    , N+(-i-j) if  (-i-j < 0L) else  -i-j ]  = bispectrum[i,j]

                bispectrum[ N+(-i)   if (-i < 0L)   else  -i   , N+(-j)   if  (-j < 0L)   else   -j   ] = np.conjugate(bispectrum[i,j])
                bispectrum[ N+(-j)   if (-j < 0L)   else  -j   , N+(-i)   if  (-i < 0L)   else   -i   ] = np.conjugate(bispectrum[i,j])
                bispectrum[ N+(-j)   if (-j < 0L)   else  -j   , N+(i+j)  if  (i+j < 0L)  else   i+j  ] = np.conjugate(bispectrum[i,j])
                bispectrum[ N+(i+j)  if (i+j < 0L)  else  i+j  , N+(-j)   if  (-j < 0L)   else   -j   ] = np.conjugate(bispectrum[i,j])
                bispectrum[ N+(i+j)  if (i+j < 0L)  else  i+j  , N+(-i)   if  (-i < 0L)   else   -i   ] = np.conjugate(bispectrum[i,j])
                bispectrum[ N+(-i)   if (-i < 0L)   else  -i   , N+(i+j)  if  (i+j < 0L)  else   i+j  ] = np.conjugate(bispectrum[i,j])

            except IndexError :
                print "(i,j) = ",i,j
                #print "Property ",symProp
                sys.exit(0)

"""
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Routine to use the symmetry properties of the bispectrum of a real 
;;; sequence to populate an entire NxN bicoherence array from the primary 
;;; octant. Takes as input an NxN real array, and the array size N, and 
;;; returns the fully populated array in the input array
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO bicoherenceSymmetry,bicoherence,N

    FOR j=0L,N/4L DO BEGIN
        FOR i=j,(N/2L-j) DO BEGIN
            bicoherence[ (j LT 0L)    ? N+(j)    : j    , (i LT 0L)    ? N+(i)    : i   ] = bicoherence[i,j]
            bicoherence[ (j LT 0L)    ? N+(j)    : j    , (-i-j LT 0L) ? N+(-i-j) : -i-j] = bicoherence[i,j]
            bicoherence[ (-i-j LT 0L) ? N+(-i-j) : -i-j , (j LT 0L)    ? N+(j)    : j   ] = bicoherence[i,j]
            bicoherence[ (-i-j LT 0L) ? N+(-i-j) : -i-j , (i LT 0L)    ? N+(i)    : i   ] = bicoherence[i,j]
            bicoherence[ (i LT 0L)    ? N+(i)    : i    , (-i-j LT 0L) ? N+(-i-j) : -i-j] = bicoherence[i,j]

            bicoherence[ (-i LT 0L)   ? N+(-i)   : -i   , (-j LT 0L)  ? N+(-j)  : -j    ] = bicoherence[i,j]
            bicoherence[ (-j LT 0L)   ? N+(-j)   : -j   , (-i LT 0L)  ? N+(-i)  : -i    ] = bicoherence[i,j]
            bicoherence[ (-j LT 0L)   ? N+(-j)   : -j   , (i+j LT 0L) ? N+(i+j) : i+j   ] = bicoherence[i,j]
            bicoherence[ (i+j LT 0L)  ? N+(i+j)  : i+j  , (-j LT 0L)  ? N+(-j)  : -j    ] = bicoherence[i,j]
            bicoherence[ (i+j LT 0L)  ? N+(i+j)  : i+j  , (-i LT 0L)  ? N+(-i)  : -i    ] = bicoherence[i,j]
            bicoherence[ (-i LT 0L)   ? N+(-i)   : -i   , (i+j LT 0L) ? N+(i+j) : i+j   ] = bicoherence[i,j]
        ENDFOR
    ENDFOR
END
"""
def bicoherenceSymmetry(bicoherence,N):
    for j in range(N/4+1):
        for i in range(N/2-j+1):
            try :

                bicoherence[ N+(j)    if (j < 0L)    else  j    , N+(i)    if  (i < 0L)    else  i    ]  = bicoherence[i,j]
                bicoherence[ N+(j)    if (j < 0L)    else  j    , N+(-i-j) if  (-i-j < 0L) else  -i-j ]  = bicoherence[i,j]
                bicoherence[ N+(-i-j) if (-i-j < 0L) else  -i-j , N+(j)    if  (j < 0L)    else  j    ]  = bicoherence[i,j]
                bicoherence[ N+(-i-j) if (-i-j < 0L) else  -i-j , N+(i)    if  (i < 0L)    else  i    ]  = bicoherence[i,j]
                bicoherence[ N+(i)    if (i < 0L)    else  i    , N+(-i-j) if  (-i-j < 0L) else  -i-j ]  = bicoherence[i,j]

                bicoherence[ N+(-i)   if (-i < 0L)   else  -i   , N+(-j)   if  (-j < 0L)   else   -j   ] = bicoherence[i,j]
                bicoherence[ N+(-j)   if (-j < 0L)   else  -j   , N+(-i)   if  (-i < 0L)   else   -i   ] = bicoherence[i,j]
                bicoherence[ N+(-j)   if (-j < 0L)   else  -j   , N+(i+j)  if  (i+j < 0L)  else   i+j  ] = bicoherence[i,j]
                bicoherence[ N+(i+j)  if (i+j < 0L)  else  i+j  , N+(-j)   if  (-j < 0L)   else   -j   ] = bicoherence[i,j]
                bicoherence[ N+(i+j)  if (i+j < 0L)  else  i+j  , N+(-i)   if  (-i < 0L)   else   -i   ] = bicoherence[i,j]
                bicoherence[ N+(-i)   if (-i < 0L)   else  -i   , N+(i+j)  if  (i+j < 0L)  else   i+j  ] = bicoherence[i,j]

            except IndexError :
                print "(i,j) = ",i,j
                #print "Property ",symProp
                sys.exit(0)


"""
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Routine to use the symmetry properties of the bicovariance of a real 
;;; sequence to populate an entire NxN array from the primary sextant. Takes
;;; as input an NxN complex array, and the array size N, and returns the 
;;; fully populated array in the input array
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO biCovarianceSymmetry,biCovariance,N

    FOR i=0L,N/2L DO BEGIN
        FOR j=0L,i DO BEGIN
            biCovariance[(j LT 0L) ? N+(j) : j , (i LT 0L) ? N+(i) : i] = biCovariance[i,j]
            biCovariance[(-j LT 0L) ? N+(-j) : -j , (i-j LT 0L) ? N+(i-j) : i-j] = biCovariance[i,j]
            biCovariance[(i-j LT 0L) ? N+(i-j) : i-j , (-j LT 0L) ? N+(-j) : -j] = biCovariance[i,j]
            biCovariance[(j-i LT 0L) ? N+(j-i) : j-i , (-i LT 0L) ? N+(-i) : -i] = biCovariance[i,j]
            biCovariance[(-i LT 0L) ? N+(-i) : -i , (j-i LT 0L) ? N+(j-i) : j-i] = biCovariance[i,j]
        ENDFOR
    ENDFOR
END

"""
def biCovarianceSymmetry(biCovariance,NN):
    for i in range(NN/2+1):
        for j in range(i+1):
            try :

                symProp=5 ; biCovariance[ NN+(j)   if (  j < 0L)  else   j , NN+(i)    if (  i < 0L) else   i ] = biCovariance[i,j]
                symProp=4 ; biCovariance[ NN+(-j)  if ( -j < 0L)  else  -j , NN+(i-j)  if (i-j < 0L) else i-j ] = biCovariance[i,j]
                symProp=3 ; biCovariance[ NN+(i-j) if (i-j < 0L)  else i-j , NN+(-j)   if ( -j < 0L) else  -j ] = biCovariance[i,j]
                symProp=2 ; biCovariance[ NN+(j-i) if (j-i < 0L)  else j-i , NN+(-i)   if ( -i < 0L) else  -i ] = biCovariance[i,j]
                symProp=1 ; biCovariance[ NN+(-i)  if ( -i < 0L)  else  -i , NN+(j-i)  if (j-i < 0L) else j-i ] = biCovariance[i,j]

            except IndexError :
                print "(i,j) = ",i,j
                print "Property ",symProp
                sys.exit(0)


def plotInstance() :
    """
        Plot a single instance of the elevation, slope, curvature and glint, for
        comparison purposes.
    """

    #left  = 0.125  # the left side of the subplots of the figure
    #right = 0.9    # the right side of the subplots of the figure
    #bottom = 0.1   # the bottom of the subplots of the figure
    #top = 0.9      # the top of the subplots of the figure
    #wspace = 0.2   # the amount of width reserved for blank space between subplots
    #hspace = 0.2   # the amount of height reserved for white space between subplots

    pl.figure(figsize=(12,10))
    pl.subplots_adjust(left=0.05,right=0.95,bottom=0.05,top=0.95,wspace=0.15,hspace=0.15)

    #pl.subplot(4,1,1)
    #pl.plot(Scale.x,totalElevSurface.real,label="Elevation")
    #pl.xlim(11.,12.)
    #pl.ylim(-0.0005,0.005)
    #pl.legend()
    #pl.show()

    pl.subplot(4,1,1)
    testAng = 0
    slopeMin = Geom.xi_min[testAng]
    slopeMax = Geom.xi_max[testAng]
    pl.plot(Scale.x,totalSlopeSurface.real,label="Slope")
    pl.plot(Scale.x,np.ones(N)*slopeMin,label="Slope Min")
    pl.plot(Scale.x,np.ones(N)*slopeMax,label="Slope Max")
    #pl.plot(Scale.x,np.gradient(totalElevSurface.real,delta_x),label="Slope (gradient)")
    #pl.plot(Scale.x[:-1],np.diff(totalElevSurface.real,n=1)/delta_x,label="Slope (diff)")
    #pl.plot(Scale.x[:-1],np.diff(totalElevSurface.real,n=1)/delta_x,label="Slope (diff)")
    pl.xlim(20.,30.)
    #pl.ylim(-0.0005,0.005)
    pl.legend()
    pl.grid(b=True)
    pl.show()

    glint = np.double(totalSlopeSurface.real > slopeMin) * np.double(totalSlopeSurface.real < slopeMax)

    pl.subplot(4,1,2)
    pl.plot(Scale.x,glint,label="Glint")
    pl.xlim(20.,30.)
    pl.ylim(-0.2,1.2)
    pl.grid(b=True)
    pl.legend()
    pl.show()

    filterLen=64
    filterStdev = 0.6
    convX = np.linspace(-5.,5.,filterLen)
    convFilter = exp(-0.5*(convX**2.)/(filterStdev**2.))
    glintConv = np.convolve(glint,convFilter,mode='same')
    glintConv = (glintConv<=1.)*glintConv + (glintConv>1.)
    SLmag,SLstdev = 0.5,1./sqrt(200.)
    Skylight = SLmag*exp(-0.5*((totalSlopeSurface.real - Geom.xi_0[testAng])**2.)/(SLstdev**2.))
    combGlint = (glintConv > Skylight)*glintConv + (glintConv < Skylight)*Skylight

    pl.subplot(4,1,3)
    pl.plot(Scale.x,combGlint,label="Combined Glint")
    pl.xlim(20.,30.)
    pl.ylim(-0.2,1.2)
    pl.grid(b=True)
    pl.legend()
    pl.show()

    threshold = 0.45
    thresholdGlint = (combGlint >= threshold)*combGlint

    pl.subplot(4,1,4)
    pl.plot(Scale.x,thresholdGlint,label="Thresholded Glint")
    pl.xlim(20.,30.)
    pl.ylim(-0.2,1.2)
    pl.legend()
    pl.grid(b=True)
    pl.show()

    nbins = 100
    bins = np.linspace(0.,1.,nbins)
    pl.figure()
    glintHistogram = np.histogram(glint,bins,normed=True)
    pl.plot(bins[:-1],glintHistogram[0],label="Glint")
    glintHistogram = np.histogram(combGlint,bins,normed=True)
    pl.plot(bins[:-1],glintHistogram[0],label="Combined Glint")
    glintHistogram = np.histogram(thresholdGlint,bins,normed=True)
    pl.plot(bins[:-1],glintHistogram[0],label="Thresholded Glint")
    pl.xlim(-0.1,1.1)
    pl.ylim(-1.,10.)
    pl.legend()
    pl.show()

    #pl.subplot(4,1,4)
    #pl.plot(Scale.x,totalCurvatureSurface,label=r"$\eta^{''}(x)$")
    #pl.plot(Scale.x,np.gradient(np.gradient(totalElevSurface.real))/(delta_x*delta_x),label=r"$\eta^{''}(x)$ (gradient)")
    #pl.plot(Scale.x[1:-1],np.diff(totalElevSurface.real,n=2)/(delta_x**2.),label="$\eta^{''}(x)$ (diff)")
    #pl.plot(Scale.x,totalCurvatureSurface/(sqrt(1. + totalSlopeSurface**2.)**3.),label=r"$\kappa(x)$")
    #pl.plot(Scale.x[1:-1],(np.diff(totalElevSurface.real,n=2)/(delta_x**2.),label="Curvature (diff)")
    #pl.xlim(11.,12.)
    #pl.ylim(-0.0005,0.005)
    #pl.legend()
    #pl.grid(b=True)
    #pl.show()


    #return

    #left  = 0.125  # the left side of the subplots of the figure
    #right = 0.9    # the right side of the subplots of the figure
    #bottom = 0.1   # the bottom of the subplots of the figure
    #top = 0.9      # the top of the subplots of the figure
    #wspace = 0.2   # the amount of width reserved for blank space between subplots
    #hspace = 0.2   # the amount of height reserved for white space between subplots

    #pl.figure(figsize=(15,12))
    #pl.subplot(2,3,1)
    #pl.subplots_adjust(left=0.05,right=0.95,bottom=0.05,top=0.95,wspace=0.15,hspace=0.15)
    #pl.plot(Scale.k,ElevPower.primaryPower,label="primary")
    #pl.plot(Scale.k,ElevPower.nlPower,label="nonLinear")
    #pl.plot(Scale.k,ElevPower.totalPower,label="total")
    #pl.xlim(0.,5.)
    #pl.ylim(-0.0005,0.005)
    #pl.legend()
    #pl.title("Elevation Power Spectrum")
    #pl.show()

    #pl.subplot(2,3,2)
    #pl.plot(Scale.k,SlopePower.primaryPower,label="primary")
    #pl.plot(Scale.k,SlopePower.nlPower,label="nonLinear")
    #pl.plot(Scale.k,SlopePower.totalPower,label="total")
    #pl.xlim(0.,5.)
    #pl.ylim(-0.0005,0.005)
    #pl.legend()
    #pl.title("Slope Power Spectrum")
    #pl.show()

    #pl.subplot(2,3,3)
    #pl.plot(Scale.k,CurvaturePower.primaryPower,label="primary")
    #pl.plot(Scale.k,CurvaturePower.nlPower,label="nonLinear")
    #pl.plot(Scale.k,CurvaturePower.totalPower,label="total")
    #pl.xlim(0.,5.)
    #pl.ylim(-0.0005,0.005)
    #pl.legend()
    #pl.title("Curvature Power Spectrum")
    #pl.show()

    #pl.subplot(2,3,4)
    #pl.plot(Scale.k,2.*totalElevAvgPower/(delta_k*N_r),label="total")
    #pl.xlim(0.,5.)
    #pl.ylim(-0.0005,0.005)
    #pl.legend()
    #pl.title("Average Elevation Power Spectrum")
    #pl.show()

    #pl.subplot(2,3,5)
    #pl.plot(Scale.k,2.*totalSlopeAvgPower/(delta_k*N_r),label="total")
    #pl.xlim(0.,5.)
    #pl.ylim(-0.0005,0.005)
    #pl.legend()
    #pl.title("Average Slope Power Spectrum")
    #pl.show()

    #pl.subplot(2,3,6)
    #pl.plot(Scale.k,2.*totalCurvatureAvgPower/(delta_k*N_r),label="total")
    #pl.xlim(0.,5.)
    #pl.ylim(-0.0005,0.005)
    #pl.legend()
    #pl.title("Average Curvature Power Spectrum")
    #pl.show()

