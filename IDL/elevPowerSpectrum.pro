;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   Function to calculate a Phillips elevation power spectrum,   ;;;
;;;   and the desired phase relationships                          ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO phillips_elev_spectrum,SCALE,POWER,NLCOUPLING,powExponent
	
	COMMON DEBUGGING,windowIndex
	
	PRINT,'Calculating the Phillips power spectrum...', Format='(/A)'

	N = SCALE.N
	dk = SCALE.delta_k
	k_max = SCALE.k_max

;	PRINT,'k_max = ',k_max,' meters^{-1}', Format='(A25,F20.10,A)'
;	PRINT,'dk = ',dk,' meters^{-1}', Format='(A25,F20.10,A)'
;	PRINT,'Nyquist Wavenumber = ',DOUBLE(N/2L)*dk,' meters^{-1}',Format='(A25,F20.12,A)'

	B      =  0.005D	;;; Dimensionless constant
	g      =  9.81D		;;; Acceleration due to gravity in m * s^{-2}
	v      =  5.D		;;; Wind velocity in m * s^{-1}
	gammaT =  0.073D    ;;; Water surface tension in N * m^{-1}
	rho    =  1000.D    ;;; Water density in kg * m^{3}

	kc0 = g/(v*v)							;;; Lower wavenumber cutoff for wind waves
	k_gamma = MIN(SQRT(rho*g/gammaT),k_max)	;;; Upper wavenumber cutoff for the capillary wave region

	nlCutoff   = 1.5D	;;; Upper limit of source waves is nlCutoff*kc0
	elevStDev  = SQRT(0.25D*B*(kc0^(-2.D)-k_gamma^(-2.D)))
	elevVar    = 0.25D*B*(kc0^(-2.D)-k_gamma^(-2.D))
	slopeStDev = SQRT(0.5D*B*ALOG(k_gamma/kc0))
	slopeVar   = 0.5D*B*ALOG(k_gamma/kc0)

	k = DBLARR(LONG(N))
    k = DINDGEN(LONG(N))*dk

	PRINT,''
	PRINT,'Theoretical elevation stdev:    ',elevStDev,' meters',Format='(A,F12.6,A)'
	PRINT,'Theoretical elevation variance: ',elevVar,' meters^{2}',Format='(A,F12.6,A)'
	PRINT,'Theoretical slope stdev:        ',slopeStDev,Format='(A,F12.6,A)'
	PRINT,'Theoretical slope variance:     ',slopeVar,Format='(A,F12.6,A)'
	PRINT,''

	PRINT,'kc0 cutoff is:     ',kc0,' meters^{-1}',Format='(A,F12.6,A)'
	PRINT,'nlCutoff is:       ',nlCutoff,Format='(A,F12.3)'
	PRINT,'(nlCutoff*kc0) is: ',nlCutoff*kc0,' meters^{-1}',Format='(A,F12.6,A)'
	PRINT,'k_gamma cutoff is: ',k_gamma,' meters^{-1}',Format='(A,F12.6,A/)'

	;;; Determine the indices of components in the wavenumber range for the free waves
	sourceIndex = WHERE((k[0:N-1] GT kc0)*(k[0:N-1] LT nlCutoff*kc0))

	PRINT,"Indices between kc0 and nlCutoff*kc0 (free1 and free2, the primary power indicies)..."
	PRINT,"sourceindex...",sourceIndex," at wavenumber ",k[sourceIndex]

	;;; Define the structure containing the phase relationships between the free
	;;; and bound waves
	NLCOUPLING = {Nbound:N_ELEMENTS(sourceIndex),bound:LONARR(N_ELEMENTS(sourceIndex)),$
		free1:LONARR(N_ELEMENTS(sourceIndex)),free2:LONARR(N_ELEMENTS(sourceIndex))}

	PRINT,"There are ",NLCOUPLING.Nbound," bound-free pairs",Format='(/A,I3,A/)'

	;;; Determine the indices of the bound waves
	coupleIndex = 2L*sourceIndex
	;coupleIndex = 2L*WHERE((k[0:N-1] GT kc0)*(k[0:N-1] LT 1.5D*kc0))

	PRINT,"Indices between 2*kc0 and 2*nlCutoff*kc0 (bound, the coupled power indicies)..."
	PRINT,"coupleIndex...",coupleIndex," at wavenumber ",k[coupleIndex]

	;;; Assign the free and bound wave indicies to the NLCOUPLING structure
	NLCOUPLING.free1 = sourceIndex
	NLCOUPLING.free2 = sourceIndex
	NLCOUPLING.bound = coupleIndex

	powExponent = DOUBLE(powExponent)   ;-- power law spectra are proportional to k^(-powExponent)
	nlProp      = 0.65D	;-- proportion of the total wave power due to bound waves

	;;; Compute the total power S(k)
	totalPower = DBLARR(N)
	totalPower[0] = 0.D
	totalPower[1:N-1] = (k[1:N-1] GT kc0)*(k[1:N-1] LT k_gamma)*B/(2.D*(k[1:N-1]^powExponent))

	;;; Set the bound power at the bound wavenumbers to 35% of the total power
	;;; at those wavenumbers
	nlPower = DBLARR(N)
	nlPower[0L] = 0.D
	nlPower[coupleIndex] = nlProp*B/(2.D*(k[coupleIndex]^powExponent))
	POWER.nlPower = nlPower

	;;; Define the primary power spectrum
	primaryPower = DBLARR(N)
	primaryPower = totalPower - nlPower
	POWER.primaryPower = primaryPower

	;xwinsize=800
	;ywinsize=450
     ;!P.MULTI=0
	;WINDOW,0,xsize = xwinsize,ysize = ywinsize,title='Elevation Power Spectrum (sub)',RETAIN=2
	;PLOT,k,totalPower,xrange=[0.0D,4.0],xtitle='k',/xstyle,$
		;ytitle='power',linestyle=0,/ystyle,charsize=chsize
	;OPLOT,k,totalPower,PSYM=2
	;OPLOT,k,POWER.primaryPower,linestyle=1
	;OPLOT,k,POWER.nlPower,linestyle=2
	;OPLOT,k[sourceIndex],totalPower[sourceIndex],PSYM=6
END

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

