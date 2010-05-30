import sys, time, string, getopt, copy
import numpy as np
from numpy import pi,sin,cos,tan,sqrt,abs
from numpy.fft import fft,ifft
from numpy import float64 as double

import matplotlib as mpl
from matplotlib import pylab as pl

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

    
    #left  = 0.125  # the left side of the subplots of the figure
    #right = 0.9    # the right side of the subplots of the figure
    #bottom = 0.1   # the bottom of the subplots of the figure
    #top = 0.9      # the top of the subplots of the figure
    #wspace = 0.2   # the amount of width reserved for blank space between subplots
    #hspace = 0.2   # the amount of height reserved for white space between subplots

    #pl.figure(figsize=(12,8))
    #pl.subplot(2,2,1)
    #pl.subplots_adjust(left=0.05,right=0.95,bottom=0.05,top=0.95,wspace=0.15,hspace=0.15)
    #pl.plot(Scale.k,ElevPower.primaryPower,label="primary")
    #pl.plot(Scale.k,ElevPower.nlPower,label="nonLinear")
    #pl.plot(Scale.k,ElevPower.totalPower,label="total")
    #pl.xlim(0.,3.)
    #pl.ylim(-0.0005,0.005)
    #pl.legend()
    #pl.title("Elevation Power Spectrum")
    #pl.show()

    print '\nFirst component indicies for free waves ...'
    print NLCoupling.free1
    print Scale.k[NLCoupling.free1]/Scale.k_N

    print 'Second component indicies for free waves ...'
    print NLCoupling.free2
    print Scale.k[NLCoupling.free2]/Scale.k_N

    print 'Indicies for bound waves...'
    print NLCoupling.bound
    print Scale.k[NLCoupling.bound]/Scale.k_N

    totalElevPower = ElevPower.totalPower
    primaryElevPower = ElevPower.primaryPower
    nlElevPower = ElevPower.nlPower

    print "\nElevation stdev from power vector: %10.6f meters " % \
        (sqrt(np.sum(totalElevPower)*delta_k))
    print "Elevation variance from power vector: %10.6f meters^{2} " % \
        (np.sum(totalElevPower)*delta_k)

    print "\nTotal elevation power at the bound wavenumbers..."
    print totalElevPower[NLCoupling.bound]
    print "Free elevation power at the bound wavenumbers..."
    print ElevPower.primaryPower[NLCoupling.bound]
    print "Bound elevation power at the bound wavenumbers..."
    print ElevPower.nlPower[NLCoupling.bound]
    print "Ratio of bound to free elevation power at the bound wavenumbers..."
    print ElevPower.nlPower[NLCoupling.bound]/totalElevPower[NLCoupling.bound]

    """
        Initialise the slope power spectrum structure
    """
    SlopePower = copy.deepcopy(ElevPower)
    SlopePower.primaryPower = k*k*ElevPower.primaryPower
    SlopePower.nlPower = k*k*ElevPower.nlPower
    SlopePower.totalPower = k*k*ElevPower.totalPower

    #pl.subplot(2,2,2)
    #pl.plot(Scale.k,SlopePower.primaryPower,label="primary")
    #pl.plot(Scale.k,SlopePower.nlPower,label="nonLinear")
    #pl.plot(Scale.k,SlopePower.totalPower,label="total")
    #pl.xlim(0.,3.)
    #pl.ylim(-0.0005,0.005)
    #pl.legend()
    #pl.title("Slope Power Spectrum")
    #pl.show()

    totalSlopePower = SlopePower.totalPower
    primarySlopePower = SlopePower.primaryPower
    nlSlopePower = SlopePower.nlPower

    print "\nSlope stdev from power vector: %10.6f meters " % \
        (sqrt(np.sum(totalSlopePower)*delta_k))
    print "Slope variance from power vector: %10.6f meters^{2} " % \
        (np.sum(totalSlopePower)*delta_k)

    print "\nTotal slope power at the bound wavenumbers..."
    print totalSlopePower[NLCoupling.bound]
    print "Free slope power at the bound wavenumbers..."
    print SlopePower.primaryPower[NLCoupling.bound]
    print "Bound slope power at the bound wavenumbers..."
    print SlopePower.nlPower[NLCoupling.bound]
    print "Ratio of bound to free slope power at the bound wavenumbers..."
    print SlopePower.nlPower[NLCoupling.bound]/totalSlopePower[NLCoupling.bound]

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

    print "\nSlope stdev from amplitude vector: %10.6f meters " % \
            (sqrt(np.sum(totalSlopeAmplitude**2.)))
    print "Slope variance from amplitude vector: %10.6f meters^{2}" % \
            (np.sum(totalSlopeAmplitude**2.))

    testSlopePhase = np.random.rand(N)*2.*pi - pi
    totalSlopeSpectrum = totalSlopeAmplitude*(-sin(testSlopePhase) + 1j*cos(testSlopePhase))
    totalSlopeSpectrum[N/2+1 :] = np.conjugate(totalSlopeSpectrum[1L : N/2][::-1])
    totalSlopeSurface = fft(totalSlopeSpectrum)

    print "\nSlope mean from surface:    %10.6f meters " % \
            np.mean(totalSlopeSurface.real)
    print "Slope stdev from surface:    %10.6f meters " % \
            np.std(totalSlopeSurface.real)
    print "Slope variance from surface: %10.6f meters^{2} " % \
            np.var(totalSlopeSurface.real)

    totalSlopeAvgPower = np.zeros(N,dtype=double)
    primarySlopeAvgPower = np.zeros(N,dtype=double)
    nlSlopeAvgPower = np.zeros(N,dtype=double)

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
    glintStats = [DataStatsStruct(numMoments) for geoms in np.arange(Geom.N_angles) ]

    """
        Loop through the surface realisations for the quadratically
        coupled oscillations
    """

    seed = 30
    N_r_cum = 0
    angleRuns = np.zeros(Geom.N_angles,np.long)
    angleRunsCum = np.zeros(Geom.N_angles,np.long)

    #while (N_r_cum < N_r) :
    while (angleRuns.sum() < N_r*Geom.N_angles) :

        N_r_cum += 1

        print "\n>>>>>>>>>>>>>>>>>>>>>\n"
        print "Computing realisation: %d" % (N_r_cum)

		### Compute the independent phases for this realisation
        primaryElevPhase = np.random.rand(N)*2.*pi - pi # RANDOMU(seed,N)*2.D*!DPI - !DPI
        nlElevPhase      = np.random.rand(N)*2.*pi - pi # RANDOMU(seed,N)*2.D*!DPI - !DPI

        ### Apply the phase correlations between the free and bound wavenumbers for the nonlinear
        ### component, if (nlSwitch==1)
        if (nlSwitch) :
            nlElevPhase[NLCoupling.bound] = primaryElevPhase[NLCoupling.free1] + \
                    primaryElevPhase[NLCoupling.free2]

        """
		    Compute the elevation realisations from the elevation spectra
		    and the synthesised phases
        """
		
		### Calculate the elevation spectrum for the free waves
        primaryElevSpectrum = primaryElevAmplitude*(cos(primaryElevPhase) + 1j*sin(primaryElevPhase))
        #primaryElevSpectrum += 0.00001*MAX(totalElevAmplitude)*RANDOMN(seed,N)
        primaryElevSpectrum[N/2+1:] = np.conjugate(primaryElevSpectrum[1 : N/2][::-1])

        ### Calculate the elevation spectrum for the bound waves
        nlElevSpectrum = nlElevAmplitude*(cos(nlElevPhase) + 1j*sin(nlElevPhase))
        nlElevSpectrum[N/2+1:] = np.conjugate(nlElevSpectrum[1:N/2][::-1])

		### Compute specific realisation of the free and bound waves. Nonlinear elevation
		### (totalElevSurface) is sum of free and bound waves.
        primaryElevSurface = fft(primaryElevSpectrum)                      ### Free waves
        nlElevSurface = fft(nlElevSpectrum)                                ### Bound waves
        totalElevSurface = primaryElevSurface + nlElevSurface              ### Total surface
 
        #print "\n\tElevation mean from surface:     %10.6e meters" % \
                #np.mean(totalElevSurface.real)
        #print "\tElevation stdev from surface:    %10.6e meters" % \
                #np.std(totalElevSurface.real)
        #print "\tElevation variance from surface: %10.6e meters^{2}\n" % \
                #np.var(totalElevSurface.real)

		### Compute the average power spectrum for free, bound and total elevation waves
        primaryElevAvgPower += abs(ifft(primaryElevSurface))**2.
        nlElevAvgPower += abs(ifft(nlElevSurface))**2.
        totalElevAvgPower += abs(ifft(totalElevSurface))**2.

        #print "\tElevation stdev from power vector:    %10.6e meters" % \
            #(sqrt(np.sum(abs(ifft(totalElevSurface.real))**2.)))
        #print "\tElevation variance from power vector: %10.6e meters^{2}\n" % \
            #(np.sum(abs(ifft(totalElevSurface.real))**2.))

        ### Compute the elevation moments

        #print "\tElevation mean from elevStats.moments:     %10.6e" % double(np.sum(totalElevSurface.real)/double(N))
        #print "\tElevation stdev from elevStats.moments:    %10.6e" % double(sqrt(np.sum(totalElevSurface.real**2.)/double(N)))
        #print "\tElevation variance from elevStats.moments: %10.6e" % double(np.sum(totalElevSurface.real**2.)/double(N))
        #print "\tElevation skewness from elevStats.moments: %10.6e" % double(np.sum(totalElevSurface.real**3.)/double(N))

        elevStats.moments += [ np.sum(totalElevSurface.real    )/double(N), \
                               np.sum(totalElevSurface.real**2.)/double(N), \
                               np.sum(totalElevSurface.real**3.)/double(N)]

        """
		    Compute the slope realisations from the slope spectra
		    and the synthesisised phases
        """
		
        ### Calculate the slope spectrum for the free waves
        primarySlopeSpectrum = primarySlopeAmplitude*(-sin(primaryElevPhase) + 1j*cos(primaryElevPhase))
        #primarySlopeSpectrum += 0.00001*MAX(totalSlopeAmplitude)*RANDOMN(seed,N)
        primarySlopeSpectrum[N/2+1:] = np.conjugate(primarySlopeSpectrum[1 : N/2][::-1])

        ### Calculate the slope spectrum for the bound waves
        nlSlopeSpectrum = nlSlopeAmplitude*(-sin(nlElevPhase) + 1j*cos(nlElevPhase))
        nlSlopeSpectrum[N/2+1:] = np.conjugate(nlSlopeSpectrum[1:N/2][::-1])

        ### Compute specific realisation of the free and bound waves. Nonlinear slope
        ### (totalSlopeSurface) is sum of free and bound waves.
        primarySlopeSurface = fft(primarySlopeSpectrum)                    ### Free waves
        nlSlopeSurface = fft(nlSlopeSpectrum)                              ### Bound waves
        totalSlopeSurface = primarySlopeSurface + nlSlopeSurface           ### Total surface

        #print "\n\tSlope mean from surface:     %10.6e meters " % \
                #np.mean(totalSlopeSurface.real)
        #print "\tSlope stdev from surface:    %10.6e meters " % \
                #np.std(totalSlopeSurface.real)
        #print "\tSlope variance from surface: %10.6e meters^{2} " % \
                #np.var(totalSlopeSurface.real)

        ### Compute the average power spectrum for free, bound and total elevation waves
        primarySlopeAvgPower += abs(ifft(primarySlopeSurface))**2.
        nlSlopeAvgPower += abs(ifft(nlSlopeSurface))**2.
        totalSlopeAvgPower += abs(ifft(totalSlopeSurface))**2.

        #print "\n\tSlope stdev from power vector:    %10.6e meters" % \
            #(sqrt(np.sum(abs(ifft(totalSlopeSurface.real))**2.)))
        #print "\tSlope variance from power vector: %10.6e meters^{2}\n" % \
            #(np.sum(abs(ifft(totalSlopeSurface.real))**2.))

        ### Compute the slope moments

        #print "\tSlope mean from elevStats.moments:     %10.6e" % double(np.sum(totalSlopeSurface.real)/double(N))
        #print "\tSlope stdev from elevStats.moments:    %10.6e" % double(sqrt(np.sum(totalSlopeSurface.real**2.)/double(N)))
        #print "\tSlope variance from elevStats.moments: %10.6e" % double(np.sum(totalSlopeSurface.real**2.)/double(N))
        #print "\tSlope skewness from elevStats.moments: %10.6e" % double(np.sum(totalSlopeSurface.real**3.)/double(N))

        slopeStats.moments += [ np.sum(totalSlopeSurface    )/double(N), \
                                np.sum(totalSlopeSurface**2.)/double(N), \
                                np.sum(totalSlopeSurface**3.)/double(N) ]

		### Compute the Fourier spectra of the total surfaces
        totalElevSpectrum = ifft(totalElevSurface)
        totalSlopeSpectrum = ifft(totalSlopeSurface)

        """
		    Loop through the geometries in the GEOM structure
        """

        for angle in np.arange(Geom.N_angles) :

            """
            Check if we have finished processing for this
            angle.
            """

            if (angleRuns[angle] < N_r) :

                print "\n\tProcessing angle ",angle," for run ",angleRuns[angle]+1, \
                " --> attempt ",angleRunsCum[angle]+1

                """
                    Compute the glint realisation from the slope
                    realisations
                """

                slopeMin = Geom.xi_min[angle]
                slopeMax = Geom.xi_max[angle]

                #print "\t(slopeMin,slopeMax) = (%f,%f)" % (slopeMin,slopeMax)

                #glint = np.array(((totalSlopeSurface.real > slopeMin) and (totalSlopeSurface.real < slopeMax)))
                glint = np.double(totalSlopeSurface.real > slopeMin) * np.double(totalSlopeSurface.real < slopeMax)

                ### Check if all glint elements vanish
                result = np.where(glint)

                #if ((result.squeeze()).shape == (0,)) :
                #if (glint.sum() == 0.) :
                if (np.shape(np.squeeze(np.where(glint))) == (0,)) :

                    print "\tZero-glint realisation angle ",angle, \
                    " for run ",angleRuns[angle]+1, \
                    " --> attempt ",angleRunsCum[angle]+1

                    ### There are no glints, add to attempts count

                    ### If this angle fails and there are no glints, then steeper 
                    ### angles will fail also, so break out of the angle loop and 
                    ### proceed to the next realisation...
                    angleRunsCum[angle:Geom.N_angles] += 1
                    break

                else :

                    print "\tSuccessful realisation angle ",angle, \
                    " for run ",angleRuns[angle] + 1, \
                    " --> attempt ",angleRunsCum[angle]+1

                    angleRuns[angle] += 1
                    angleRunsCum[angle] += 1

                    ### Compute the glint moments

                    glintStats[angle].moments += [ np.sum(glint.real    )/double(N), \
                                                   np.sum(glint.real**2.)/double(N), \
                                                   np.sum(glint.real**3.)/double(N) ]

                    ### Compute the Fourier spectrum of this glint realisation
                    glintSpectrum = ifft(glint)

                    """
                    Compute the average glint power spectrum
                    """

                    totalGlintAvgPower[angle] += abs(glintSpectrum)**2.

                ### End checking for zero-glint of this angle

            ### End checking for completion of this angle

        ### End angle loop

    ### End realisation loop

    #pl.subplot(2,2,3)
    #pl.plot(Scale.k,2.*totalElevAvgPower/(delta_k*N_r),label="total")
    #pl.xlim(0.,3.)
    #pl.ylim(-0.0005,0.005)
    #pl.legend()
    #pl.title("Average Elevation Power Spectrum")
    #pl.show()

    #pl.subplot(2,2,4)
    #pl.plot(Scale.k,2.*totalSlopeAvgPower/(delta_k*N_r),label="total")
    #pl.xlim(0.,3.)
    #pl.ylim(-0.0005,0.005)
    #pl.legend()
    #pl.title("Average Slope Power Spectrum")
    #pl.show()

    print ""
    print "AngleRuns:    ",angleRuns," ... for total of ", \
        int(np.sum(angleRuns))," / ",N_r*Geom.N_angles
    print "AngleRunsCum: ",angleRunsCum
	
    N_runs = N_r
    N_r = N_r_cum
    print  "Final N_r_cum = ",N_r_cum

    elevStats.moments /= double(N_runs)
    elevStats.cumulantsFromMoments()
    slopeStats.moments /= double(N_runs)
    slopeStats.cumulantsFromMoments()
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
    print "Slope skewness %10.6e" % slopeStats.skewness

    for geoms in np.arange(Geom.N_angles) :
        print "\nAngle ",geoms,"..."
        print "\n\tGlint first moment  %10.6e" % glintStats[geoms].moments[0]
        print "\tGlint second moment %10.6e" % glintStats[geoms].moments[1]
        print "\tGlint third moment  %10.6e" % glintStats[geoms].moments[2]

    for geoms in np.arange(Geom.N_angles) :
        print "\nAngle ",geoms,"..."
        print "\n\tGlint first cumulant  %10.6e" % glintStats[geoms].cumulants[0]
        print "\tGlint second cumulant %10.6e" % glintStats[geoms].cumulants[1]
        print "\tGlint third cumulant  %10.6e" % glintStats[geoms].cumulants[2]

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


if __name__ == '__main__':
    N            = int(sys.argv[1])
    NN           = int(sys.argv[2])   
    delta_x      = double(sys.argv[3]) 
    N_r          = int(sys.argv[4])  
    spectrumType = str(sys.argv[5])    
    specExp      = double(sys.argv[6])  
    nlSwitch     = int(sys.argv[7])  

    cumulantFunctionSimulate(N,NN,delta_x,N_r,spectrumType,specExp,nlSwitch)
