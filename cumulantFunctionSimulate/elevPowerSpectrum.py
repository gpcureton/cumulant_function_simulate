import sys, string
import numpy as np

######################################################################
###   Function to calculate a Phillips elevation power spectrum,   ###
###   and the desired phase relationships                          ###
######################################################################

def phillips_elev_spectrum(Scale,Power,NLCoupling,powExponent):
	
    print '\nCalculating the Phillips power spectrum...'

    N = Scale.N
    dk = Scale.delta_k
    k_max = Scale.k_max
    k = Scale.k

    print 'N     = %15d' % (N)
    print 'k_max = %15.6f meters^{-1}' % (k_max)
    print 'dk    = %15.6f meters^{-1}' % (dk)
    print 'k_N   = %15.6f meters^{-1}' % (np.double(N/2)*dk)

    B      =  0.005	    # Dimensionless constant
    g      =  9.81      # Acceleration due to gravity in m * s^{-2}
    v      =  5.        # Wind velocity in m * s^{-1}
    gammaT =  0.073     # Water surface tension in N * m^{-1}
    rho    =  1000.     # Water density in kg * m^{3}

    kc0 = g/(v*v)	                                # Lower wavenumber cutoff for wind waves
    k_gamma = np.min(np.sqrt(rho*g/gammaT),k_max)	# Upper wavenumber cutoff for 
                                                    #the capillary wave region

    nlCutoff   = 1.5	### Upper limit of source waves is nlCutoff*kc0
    elevStDev  = np.sqrt(0.25*B*(kc0**(-2.)-k_gamma**(-2.)))
    elevVar    = 0.25*B*(kc0**(-2.)-k_gamma**(-2.))
    slopeStDev = np.sqrt(0.5*B*np.log(k_gamma/kc0))
    slopeVar   = 0.5*B*np.log(k_gamma/kc0)

    print ''
    print 'Theoretical elevation stdev:    ',elevStDev,' meters'
    print 'Theoretical elevation variance: ',elevVar,' meters^{2}'
    print 'Theoretical slope stdev:        ',slopeStDev
    print 'Theoretical slope variance:     ',slopeVar
    print ''

    print 'kc0 cutoff is:     ',kc0,' meters^{-1}'
    print 'nlCutoff is:       ',nlCutoff
    print '(nlCutoff*kc0) is: ',nlCutoff*kc0,' meters^{-1}'
    print 'k_gamma cutoff is: ',k_gamma,' meters^{-1}'

    ### Determine the indices of components in the wavenumber range for the free waves
    sourceIndex = np.squeeze(np.where((k > kc0)*(k < nlCutoff*kc0)))

    print "Indices between kc0 and nlCutoff*kc0 (free1 and free2, the primary power indicies)..."
    print "sourceindex...",sourceIndex," at wavenumber ",k[sourceIndex]

    ### Define the structure containing the phase relationships between the free
    ### and bound waves
    NLCoupling.Nbound = np.shape(sourceIndex)[0]    # Number of coupling relationships
    NLCoupling.bound = NLCoupling.free1 = NLCoupling.free2 \
        = np.zeros(NLCoupling.Nbound)

    print 'There are %3d bound-free pairs' % (NLCoupling.Nbound)

    ### Determine the indices of the bound waves
    coupleIndex = 2*sourceIndex

    print "Indices between 2*kc0 and 2*nlCutoff*kc0 (bound, the coupled power indicies)..."
    print "coupleIndex...",coupleIndex," at wavenumber ",k[coupleIndex]

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
	print "Calculating the Elevation and Slope variances..."

	B      =  0.005	    # Dimensionless constant
	g      =  9.81      # Acceleration due to gravity in m * s^{-2}
	v      =  5.        # Wind velocity in m * s^{-1}
	gammaT =  0.073     # Water surface tension in N * m^{-1}
	rho    =  1000.     # Water density in kg * m^{3}

	k_max = k[-1]
	kc0 = g/(v*v)	                                # Lower wavenumber cutoff for wind waves
	k_gamma = np.min(np.sqrt(rho*g/gammaT),k_max)	# Upper wavenumber cutoff for the 
                                                    # capillary wave region

	nlCutoff   = 1.5	### Upper limit of source waves is nlCutoff*kc0
	elevStDev  = np.sqrt(0.25*B*(kc0**(-2.)-k_gamma**(-2.)))
	elevVar    = 0.25*B*(kc0**(-2.)-k_gamma**(-2.))
	slopeStDev = np.sqrt(0.5*B*np.log(k_gamma/kc0))
	slopeVar   = 0.5*B*np.log(k_gamma/kc0)

	return elevVar,slopeVar


"""

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   Function to calculate the gaussian elevation power spectrum   ;;;
;;;   and the desired phase relationships                           ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO gaussian_elev_spectrum,SCALE,POWER,NLCOUPLING
	COMMON DEBUGGING,windowIndex
	print,'Calculating the Gaussian power spectrum...', Format='(/A)'


	elevStDev = 0.13D ; meters
	elevVar = elevStDev*elevStDev
	slopeStDev = 0.2121D
	slopeVar = slopeStDev*slopeStDev
	elevCorr = SQRT(2.0D)*elevStDev/slopeStDev ; meters

	PRINT,''
	PRINT,'Theoretical elevation stdev:   ',elevStDev,' meters',Format='(A,F12.6,A)'
	PRINT,'Theoretical elevation variance: ',elevVar,' meters^{2}',Format='(A,F12.6,A)'
	PRINT,'Theoretical slope stdev:       ',slopeStDev,Format='(A,F12.6,A)'
	PRINT,'Theoretical slope variance:       ',slopeVar,Format='(A,F12.6,A)'
;	PRINT,'Elevation correlation length:  ',elevCorr,' meters',Format='(A,F12.6,A)'
	PRINT,''

	N = SCALE.N
	dk = SCALE.delta_k
	k_max = SCALE.k_max
	k = DBLARR(LONG(N))
    k = DINDGEN(LONG(N))*dk;


	;; Compute the total one sided power S(k)
	totalPower = DBLARR(N)
	exp_arg = -(elevCorr*elevCorr*k*k)/4.0D;
	;;; Causes arithmetic underflow! OK though
	totalPower = elevCorr*elevVar*exp(exp_arg)/SQRT(!DPI)

	;;; Determine the indices of components in the wavenumber range for the free waves
	sourceIndex=WHERE((totalPower GT 0.5D*totalPower[0L])*(totalPower LT 0.8D*totalPower[0L]))

	;;; Define the structure containing the phase relationships between the free
	;;; and bound waves
	NLCOUPLING = {Nbound:N_ELEMENTS(sourceIndex),bound:LONARR(N_ELEMENTS(sourceIndex)),$
		free1:LONARR(N_ELEMENTS(sourceIndex)),free2:LONARR(N_ELEMENTS(sourceIndex))}

	;;; Determine the indices of the bound waves
	coupleIndex = 2L*sourceIndex

	totalPower[0L] = 0.D

	NLCOUPLING.free1 = sourceIndex
	NLCOUPLING.free2 = sourceIndex
	NLCOUPLING.bound = coupleIndex

	;;; Set the bound power at the bound wavenumbers to 35% of the total power
	;;; at those wavenumbers
	nlPower = DBLARR(N)
	nlPower[0L] = 0.D
	nlPower[coupleIndex] = 0.35D*totalPower[coupleIndex]
	POWER.nlPower = nlPower

	;;; Define the primary power spectrum
	primaryPower = DBLARR(N)
	primaryPower = totalPower - nlPower
	POWER.primaryPower = primaryPower

	totalPower[0L] = 0.D

	PRINT,'Elevation variance from one-sided power: ',TOTAL(totalPower*dk),' meters^{2}',$
		Format='(A,F12.6,A)'
	PRINT,'Elevation variance from one-sided power: ',TOTAL(totalPower)*dk,' meters^{2}',$
		Format='(A,F12.6,A)'

	xwinsize=800
	ywinsize=450
 	!P.MULTI=0
	WINDOW,0,xsize = xwinsize,ysize = ywinsize,title='Elevation Power Spectrum (sub)',RETAIN=2
	PLOT,k,totalPower,xrange=[0.01D,10.0],xtitle='k',/xstyle,$
		ytitle='power',linestyle=0,/ystyle,charsize=chsize
	OPLOT,k,totalPower,PSYM=2
	OPLOT,k,POWER.primaryPower,linestyle=1
	OPLOT,k,POWER.nlPower,linestyle=2
	OPLOT,k[sourceIndex],totalPower[sourceIndex],PSYM=6

END

"""
