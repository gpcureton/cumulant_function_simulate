PRO cumulantFunctionSimulate,N,NN,delta_x,N_r,spectrumType,specExp,powWindType,bispWindType,nlSwitch

	;+
	; NAME:
	;    cumulantFunctionSimulate
	;
	; PURPOSE:
	;    This program simulates multiple elevation and slope realisations with
    ;    given spectra, introducing quadratic phase coupling at selected wavenumbers.
    ;    Glint realisations are computed from the slope realisations. The average 
    ;    elevation, slope and glint moments and cumulants, and moment and cumulant 
    ;    functions are computed. The results are written to a HDF4 file.
	;
	; CATEGORY:
	;    Main routine
	;
	; CALLING SEQUENCE:
	;
	; INPUTS:
	;    N            : Data length
	;    NN           : Bispectrum data length
	;    delta_x      : Spatial increment in meters
	;    N_r          : Number of realisations
	;    spectrumType : Form of the elevation power spectrum
	;    powWindType  : Form of the glint power spectrum windowing
	;    bispWindType : Form of the glint bispectrum windowing
	;    specExp      : Elevation power spectrum is proportional to
	;                   k^{-specExp}
	;    nlSwitch     : Elevation phase coupling on/off
	;
	; OPTIONAL INPUTS:
	;    None.
	;
	; KEYWORD PARAMETERS:
	;
	; OUTPUTS:
	;
	; OPTIONAL OUTPUTS:
	;    None
	;
	; COMMON BLOCKS:
	;    None
	;
	; SIDE EFFECTS:
	;    None.
	;
	; RESTRICTIONS:
	;
	; EXAMPLE:
	;
	;
	; Created by Geoff Cureton <geoff.cureton@ssec.wisc.edu> on 2009-04-15.
	; Copyright (c) 2009-2013 Geoff Cureton. All rights reserved.
	; 
	; file_Date = '$Date$'
	; file_Revision = '$Revision$'
	; file_Author = '$Author$'
	; file_HeadURL = '$HeadURL$'
	; file_Id = '$Id$'
	;
	;
	;     This program is free software: you can redistribute it and/or modify
	;     it under the terms of the GNU General Public License as published by
	;     the Free Software Foundation, either version 3 of the License, or
	;     (at your option) any later version.
	; 
	;     This program is distributed in the hope that it will be useful,
	;     but WITHOUT ANY WARRANTY; without even the implied warranty of
	;     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	;     GNU General Public License for more details.
	; 
	;     You should have received a copy of the GNU General Public License
	;     along with this program.  If not, see <http://www.gnu.org/licenses/>.
	;-


	;;; Turn on error reporting
	!EXCEPT=1
	;COMMON DEBUGGING,windowIndex
	windowIndex = -1

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;; Determine the various scale parameters and populate the  ;;;
	;;; SCALE structure                                          ;;;
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	x_max = FLOAT(N-1L)*delta_x ; meters
	x = FINDGEN(N)*delta_x
	k_max = 2.D*!DPI/delta_x
	delta_k = k_max/DOUBLE(N-1L)
	k_N = DOUBLE(N/2L)*delta_k ; Nyquist "wavenumber"
	k = FINDGEN(N)*delta_k

	PRINT,'N = ',N,Format='(A25,F6.0)'
	PRINT,'delta_x = ',delta_x,' meters', Format='(A25,F20.10,A)'
	PRINT,'x_max = ',x_max,' meters',Format='(A25,F20.12,A)'
	PRINT,'k_max = ',k_max,' meters^{-1}', Format='(A25,F20.10,A)'
	PRINT,'delta_k = ',delta_k,' meters^{-1}', Format='(A25,F20.10,A)'
	PRINT,'Nyquist Wavenumber = ',k_N,' meters^{-1}',Format='(A25,F20.12,A)'

	;;; Initialise the scale structure
	SCALE = {N:N,delta_x:delta_x,x_max:x_max,k_max:k_max,delta_k:delta_k}

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;;   Populate the GEOM structure with angular quantities    ;;;
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	N_angles = 5L
	GEOM = {N_angles:0L, $
			source_angle:DINDGEN(N_angles), $
			detector_angle:DINDGEN(N_angles), $
			xi_min:DINDGEN(N_angles), $
			xi_0:DINDGEN(N_angles), $
			xi_max:DINDGEN(N_angles) $
	}

	angleLo=10.D
	angleHi=30.D
	angleRange=(angleHi-angleLo)
	
	IF ((N_angles-1L) EQ 0) THEN BEGIN
		d_angle = 0L
	ENDIF ELSE BEGIN
		d_angle = angleRange/(N_angles-1L)
	ENDELSE
	
	start_angle = angleLo
	d2r = !DPI/180.D
	r2d = 180.D/!DPI
	beta = 0.68D*d2r

	GEOM.N_angles = N_angles
	GEOM.source_angle = ((start_angle + DINDGEN(N_angles)*d_angle))*d2r
	GEOM.detector_angle = 0.0D*d2r
	gamma = (GEOM.source_angle - GEOM.detector_angle)/2.D
	GEOM.xi_0 = tan(gamma)
	GEOM.xi_min = GEOM.xi_0 - (1.0D + GEOM.xi_0^2.D)*(beta/4.D)
	GEOM.xi_max = GEOM.xi_0 + (1.0D + GEOM.xi_0^2.D)*(beta/4.D)
	dxi = GEOM.xi_max - GEOM.xi_min

	angleRuns = LONARR(GEOM.N_angles)
	angleRunsCum = LONARR(GEOM.N_angles)

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;;   Populate the elevation power spectrum structure POWER,  ;;;
	;;;   and NLCOUPLING                                          ;;;
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	POWER = {spectrumType:spectrumType,primaryPower:DBLARR(N),nlPower:DBLARR(N)}	
	;POWER = {spectrumType:'phillips',primaryPower:DBLARR(N),nlPower:DBLARR(N)}
	;POWER = {spectrumType:'gaussian',primaryPower:DBLARR(N),nlPower:DBLARR(N)}
	NLCOUPLING = {Nbound:0L,bound:LONARR(N),free1:LONARR(N),free2:LONARR(N)}

	CASE POWER.spectrumType OF
	'gaussian': BEGIN
			PRINT, 'The spectrum is gaussian'
			gaussian_elev_spectrum,SCALE,POWER,NLCOUPLING
		END
	'phillips': BEGIN
			PRINT, 'The spectrum is phillips'
			phillips_elev_spectrum,SCALE,POWER,NLCOUPLING,specExp
		END
	ENDCASE

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;; Check some of the data from the POWER structure			 ;;;
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	PRINT,"First component indicies for free waves ...",FORMAT='(/A)'
	PRINT,NLCOUPLING.free1
	PRINT,k[NLCOUPLING.free1]/k_N
	
	PRINT,"Second component indicies for free waves ...",FORMAT='(/A)'
	PRINT,NLCOUPLING.free2
	PRINT,k[NLCOUPLING.free2]/k_N
	
	PRINT,"Indicies for bound waves...",FORMAT='(/A)'
	PRINT,NLCOUPLING.bound

	totalElevPower = POWER.primaryPower + POWER.nlPower
	PRINT,"Elevation stdev from power vector:    ",SQRT(TOTAL(totalElevPower)*delta_k)," meters",FORMAT='(/A,F10.6,A)'
	PRINT,"Elevation variance from power vector: ",TOTAL(totalElevPower)*delta_k," meters^{2}",FORMAT='(A,F10.6,A)'

	totalSlopePower = DBLARR(N)
	totalSlopePower = k*k*totalElevPower
	primarySlopePower = DBLARR(N)
	primarySlopePower = k*k*POWER.primaryPower
	nlSlopePower = DBLARR(N)
	nlSlopePower = k*k*POWER.nlPower

	PRINT,"Slope stdev from power vector: ",SQRT(TOTAL(totalSlopePower)*delta_k),FORMAT='(A,F10.6)'
	PRINT,"Slope variance from power vector: ",TOTAL(totalSlopePower)*delta_k,FORMAT='(A,F10.6/)'

	PRINT,"Total elevation power at the bound wavenumbers...",FORMAT='(/A)'
	PRINT,totalElevPower[NLCOUPLING.bound]
	PRINT,"Free elevation power at the bound wavenumbers..."
	PRINT,POWER.primaryPower[NLCOUPLING.bound]
	PRINT,"Bound elevation power at the bound wavenumbers..."
	PRINT,POWER.nlPower[NLCOUPLING.bound]
	PRINT,"Ratio of bound to free elevation power at the bound wavenumbers..."
	PRINT,POWER.nlPower[NLCOUPLING.bound]/totalElevPower[NLCOUPLING.bound]

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;; Compute the total elevation amplitude, phase and spectrum, ;;;
	;;; and the second moment function                             ;;;
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	totalElevAmplitude = DBLARR(N)
	totalElevAmplitude = SQRT(0.5D*totalElevPower*delta_k)
	totalElevAmplitude[N/2L+1L:N-1L] = REVERSE(totalElevAmplitude[1L:N/2L-1L])
	totalElevSpectrum = DCOMPLEXARR(N)

	primaryElevAmplitude = DBLARR(N)
	primaryElevAmplitude = SQRT(0.5D*POWER.primaryPower*delta_k)
	primaryElevAmplitude[N/2L+1L:N-1L] = REVERSE(primaryElevAmplitude[1L:N/2L-1L])

	nlElevAmplitude = DBLARR(N)
	nlElevAmplitude = SQRT(0.5D*POWER.nlPower*delta_k)
	nlElevAmplitude[N/2L+1L:N-1L] = REVERSE(nlElevAmplitude[1L:N/2L-1L])

	PRINT,"Elevation stdev from amplitude vector:    ",SQRT(TOTAL(totalElevAmplitude^2.D)),FORMAT='(/A,F10.6)'
	PRINT,"Elevation variance from amplitude vector: ",TOTAL(totalElevAmplitude^2.D),FORMAT='(A,F10.6/)'

	totalElevAvgPower = DBLARR(N)
	primaryElevAvgPower = DBLARR(N)
	nlElevAvgPower = DBLARR(N)

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;; Compute the total slope amplitude, phase and spectrum     ;;;
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	totalSlopeAmplitude = DBLARR(N)
	totalSlopeAmplitude = SQRT(0.5D*totalSlopePower*delta_k)
	totalSlopeAmplitude[N/2L+1L:N-1L] = REVERSE(totalSlopeAmplitude[1L:N/2L-1L])
	totalSlopeSpectrum = DCOMPLEXARR(N)
	totalSlopeSurface = DCOMPLEXARR(N)

	primarySlopeAmplitude = DBLARR(N)
	primarySlopeAmplitude = SQRT(0.5D*primarySlopePower*delta_k)
	primarySlopeAmplitude[N/2L+1L:N-1L] = REVERSE(primarySlopeAmplitude[1L:N/2L-1L])
	primarySlopeSpectrum = DCOMPLEXARR(N)
	primarySlopeSurface = DCOMPLEXARR(N)

	nlSlopeAmplitude = DBLARR(N)
	nlSlopeAmplitude = SQRT(0.5D*nlSlopePower*delta_k)
	nlSlopeAmplitude[N/2L+1L:N-1L] = REVERSE(nlSlopeAmplitude[1L:N/2L-1L])
	nlSlopeSpectrum = DCOMPLEXARR(N)
	nlSlopeSurface = DCOMPLEXARR(N)

	PRINT,"Slope stdev from amplitude vector: ",SQRT(TOTAL(totalSlopeAmplitude^2.D)),FORMAT='(A,F10.6)'
	PRINT,"Slope variance from amplitude vector: ",TOTAL(totalSlopeAmplitude^2.D),FORMAT='(A,F10.6/)'

	totalSlopeAvgPower = DBLARR(N)
	primarySlopeAvgPower = DBLARR(N)
	nlSlopeAvgPower = DBLARR(N)

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;;     Define the glint, glint spectrum and glint power      ;;;
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	glint = DBLARR(N)
	glintSpectrum = DCOMPLEXARR(N)
	totalGlintAvgPower = DBLARR(N,GEOM.N_angles)

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;; Define the various point estimators for the elevation,    ;;;
	;;; slope and glint                                           ;;;
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	numMoments = 3L

	elevMoments       = DBLARR(numMoments)
	slopeMoments      = DBLARR(numMoments)
	glintFirstMoments = DBLARR(GEOM.N_angles)
	
	;;; The IDL centered moments

	elevMean      = 0.D
	slopeMean     = 0.D
	glintMean     = DBLARR(GEOM.N_angles)
	elevVariance  = 0.D
	slopeVariance = 0.D
	glintVariance = DBLARR(GEOM.N_angles)
	elevSkewness  = 0.D
	slopeSkewness = 0.D
	glintSkewness = DBLARR(GEOM.N_angles)

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;; Define the various quantities for the calculation of the  ;;;
	;;; bispectrum and the component power spectra                ;;;
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	;NN = 512L				;-- Size of the "entire" bispectrum domain
	NN2 = LONG(FIX(NN/2L))	;-- The Nyquist number of the bispectrum domain
	NN4 = LONG(FIX(NN/4L))
	PRINT,"NN = ",NN,FORMAT='(/A,I8)'

	;;; The bispectrum domain abcissa
	kx = DINDGEN(NN)*delta_k
	ky = DINDGEN(NN)*delta_k
	
	;;; The bicovariance domain abcissa
	tau_x = DINDGEN(NN)*delta_x
	tau_y = DINDGEN(NN)*delta_x

	elevBispectrum =     DCOMPLEXARR(NN,NN)
	elevComponentPower = DBLARR(NN,NN)
	elevSumPower =       DBLARR(NN,NN)
	elevBicoherence =    DBLARR(NN,NN)

	slopeBispectrum =     DCOMPLEXARR(NN,NN)
	slopeComponentPower = DBLARR(NN,NN)
	slopeSumPower =       DBLARR(NN,NN)
	slopeBicoherence =    DBLARR(NN,NN)

	glintBispectrum =     DCOMPLEXARR(NN,NN,GEOM.N_angles)
	glintComponentPower = DBLARR(NN,NN,GEOM.N_angles)
	glintSumPower =       DBLARR(NN,NN,GEOM.N_angles)
	glintBicoherence =    DBLARR(NN,NN,GEOM.N_angles)
	
	;;; The Hanning window function used to ensure that the glint power
	;;; spectrum vanishes outside of the Nyquist wavenumbers.

	hanWindow_N = DBLARR(N)
	hanWindow_N = SHIFT(HANNING(N,ALPHA=0.5,/DOUBLE),N/2)

	hanWindow_small = SHIFT(HANNING(NN,ALPHA=0.5,/DOUBLE),NN/2)
	hanWindow = DBLARR(N)
	hanWindow[0L:NN/2L] = hanWindow_small[0L:NN/2L]
	hanWindow[N/2L+1L:N-1L] = REVERSE(hanWindow[1L:N/2L-1L])

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;;   Loop through the surface realisations for the quadratically   ;;;
	;;;   coupled oscillations											;;;
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	seed = 30L
	N_r_cum = 0L
	angleRuns[*] = 0L
	angleRunsCum[*] = 0L

	WHILE (TOTAL(angleRuns) LT N_r*GEOM.N_angles) DO BEGIN

		N_r_cum++

		PRINT,">>>>>>>>>>>>>>>>>>>>>",FORMAT='(/A/)'
		PRINT,"Computing realisation: ",N_r_cum,FORMAT='(A,I4/)'

		;;; Compute the independent phases for this realisation
		primaryElevPhase = RANDOMU(seed,N)*2.D*!DPI - !DPI
		nlElevPhase = RANDOMU(seed,N)*2.D*!DPI - !DPI

		;;; Apply the phase correlations between the free and bound wavenumbers for the nonlinear
		;;; component, if (nlSwitch==1)
		IF (nlSwitch EQ 1) THEN $
			nlElevPhase[NLCOUPLING.bound] = primaryElevPhase[NLCOUPLING.free1] + primaryElevPhase[NLCOUPLING.free2]

		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		;;;   Compute the elevation realisations from the elevation spectra ;;;
		;;;   and the synthesisised phases                                  ;;;
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		
		;;; Calculate the elevation spectrum for the free waves
		primaryElevSpectrum = primaryElevAmplitude*DCOMPLEX(COS(primaryElevPhase),SIN(primaryElevPhase))
;		primaryElevSpectrum += 0.00001*MAX(totalElevAmplitude)*RANDOMN(seed,N)
		primaryElevSpectrum[N/2L+1L:N-1L] = CONJ(REVERSE(primaryElevSpectrum[1L:N/2L-1L]))

		;;; Calculate the elevation spectrum for the bound waves
		nlElevSpectrum = nlElevAmplitude*DCOMPLEX(COS(nlElevPhase),SIN(nlElevPhase))
		nlElevSpectrum[N/2L+1L:N-1L] = CONJ(REVERSE(nlElevSpectrum[1L:N/2L-1L]))

		;;; Compute specific realisation of the free and bound waves. Nonlinear elevation
		;;; (totalElevSurface) is sum of free and bound waves.
		primaryElevSurface = FFT(primaryElevSpectrum,/INVERSE,/DOUBLE)		;;; Free waves
		nlElevSurface = FFT(nlElevSpectrum,/INVERSE,/DOUBLE)				;;; Bound waves
		totalElevSurface = primaryElevSurface + nlElevSurface				;;; Total surface

		;;; Compute the average power spectrum for free, bound and total elevation waves
		primaryElevAvgPower += ABS(FFT(primaryElevSurface,/DOUBLE))^2.D
		nlElevAvgPower += ABS(FFT(nlElevSurface,/DOUBLE))^2.D
		totalElevAvgPower += ABS(FFT(totalElevSurface,/DOUBLE))^2.D

		;;; Compute the elevation moments
		elevMean     += MEAN(DOUBLE(totalElevSurface),/DOUBLE)
		elevVariance += VARIANCE(DOUBLE(totalElevSurface),/DOUBLE)
		elevSkewness += SKEWNESS(DOUBLE(totalElevSurface),/DOUBLE)

		elevMoments += [	TOTAL(DOUBLE(totalElevSurface))/DOUBLE(N), $
							TOTAL(DOUBLE(totalElevSurface)^2.D)/DOUBLE(N), $ 
							TOTAL(DOUBLE(totalElevSurface)^3.D)/DOUBLE(N) ]

		;;; Compute the Fourier spectra of the total surfaces
		totalElevSpectrum = FFT(totalElevSurface,/DOUBLE)

		;;; Calculate the average bispectra (for the reduced domain)
		FOR j=0L,NN4 DO BEGIN
			FOR i=j,NN2-j DO BEGIN
				elevBispectrum[i,j] += totalElevSpectrum[i]*totalElevSpectrum[j]*CONJ(totalElevSpectrum[i+j])
			ENDFOR
		ENDFOR

		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		;;;   Compute the slope realisations from the slope spectra ;;;
		;;;   and the synthesisised phases                          ;;;
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		
		;;; Calculate the slope spectrum for the free waves
		primarySlopeSpectrum = primarySlopeAmplitude*DCOMPLEX(-SIN(primaryElevPhase),COS(primaryElevPhase))
		primarySlopeSpectrum[N/2L+1L:N-1L] = CONJ(REVERSE(primarySlopeSpectrum[1L:N/2L-1L]))

		;;; Calculate the slope spectrum for the bound waves
		nlSlopeSpectrum = nlSlopeAmplitude*DCOMPLEX(-SIN(nlElevPhase),COS(nlElevPhase))
		nlSlopeSpectrum[N/2L+1L:N-1L] = CONJ(REVERSE(nlSlopeSpectrum[1L:N/2L-1L]))

		;;; Compute specific realisation of the free and bound waves. Nonlinear slope
		;;; (totalSlopeSurface) is sum of free and bound waves.
		primarySlopeSurface = FFT(primarySlopeSpectrum,/INVERSE,/DOUBLE)    ;;; Free waves
		nlSlopeSurface = FFT(nlSlopeSpectrum,/INVERSE,/DOUBLE)              ;;; Bound waves
		totalSlopeSurface = primarySlopeSurface + nlSlopeSurface            ;;; Total surface

		;;; Compute the average power spectrum for free, bound and total elevation waves
		primarySlopeAvgPower += ABS(FFT(primarySlopeSurface,/DOUBLE))^2.D
		nlSlopeAvgPower += ABS(FFT(nlSlopeSurface,/DOUBLE))^2.D
		totalSlopeAvgPower += ABS(FFT(totalSlopeSurface,/DOUBLE))^2.D

		;;; Compute the slope moments
		slopeMean += MEAN(DOUBLE(totalSlopeSurface),/DOUBLE)
		slopeVariance += VARIANCE(DOUBLE(totalSlopeSurface),/DOUBLE)
		slopeSkewness += SKEWNESS(DOUBLE(totalSlopeSurface),/DOUBLE)

		slopeMoments += [	TOTAL(DOUBLE(totalSlopeSurface))/DOUBLE(N), $
							TOTAL(DOUBLE(totalSlopeSurface)^2.D)/DOUBLE(N), $ 
							TOTAL(DOUBLE(totalSlopeSurface)^3.D)/DOUBLE(N) ]

		;;; Compute the Fourier spectrum of the total surface
		totalSlopeSpectrum = FFT(totalSlopeSurface,/DOUBLE)

		;;; Calculate the average bispectra (for the reduced domain)
		FOR j=0L,NN4 DO BEGIN
			FOR i=j,NN2-j DO BEGIN
				slopeBispectrum[i,j] += totalSlopeSpectrum[i]*totalSlopeSpectrum[j]*CONJ(totalSlopeSpectrum[i+j])
			ENDFOR
		ENDFOR

		;;; Calculate the average component power spectra (for the reduced domain)
		FOR j=0L,NN4 DO BEGIN
			FOR i=j,NN2-j DO BEGIN
				elevComponentPower[i,j] += (ABS(totalElevSpectrum[i]*totalElevSpectrum[j]))^2.D
				slopeComponentPower[i,j] += (ABS(totalSlopeSpectrum[i]*totalSlopeSpectrum[j]))^2.D
			ENDFOR
		ENDFOR

		;;; Calculate the average sum power spectra (for the reduced domain)
		FOR j=0L,NN4 DO BEGIN
			FOR i=j,NN2-j DO BEGIN
				elevSumPower[i,j] += (ABS(totalElevSpectrum[i+j]))^2.D
				slopeSumPower[i,j] += (ABS(totalSlopeSpectrum[i+j]))^2.D
			ENDFOR
		ENDFOR

		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		;;;   Loop through the geometries in the GEOM structure             ;;;
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

		FOR angle=0L, GEOM.N_angles-1L DO BEGIN

			;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
			;;;   Check if we have finished processing for this   ;;;
			;;;   angle.                                          ;;;
			;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
			
			IF (angleRuns[angle] LT N_r) THEN BEGIN

				PRINT,"Processing angle ",angle," for run ",angleRuns[angle] +1L, $
					" --> attempt ",angleRunsCum[angle]+1L,FORMAT='(A,I3,2(A,I7))'

				;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
				;;;   Compute the glint realisation from the slope    ;;;
				;;;   realisations                                    ;;;
				;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

				slopeMin = GEOM.xi_min[angle]
				slopeMax = GEOM.xi_max[angle]
				glint = (REAL_PART(totalSlopeSurface) GT slopeMin)*$
					(REAL_PART(totalSlopeSurface) LT slopeMax)
				;glint = glint - MEAN(glint)

				;;; Check if all glint elements vanish
				result=WHERE(glint,NCOMPLEMENT=vanishCount)
				
				IF(vanishCount EQ N) THEN BEGIN
					
					PRINT,"Zero-glint realisation angle ",angle, $
						" for run ",angleRuns[angle] + 1L, $
						" --> attempt ",angleRunsCum[angle]+1L,FORMAT='(A,I3,2(A,I7))'

					;;; There are no glints, add to attempts count
					angleRunsCum[angle:GEOM.N_angles-1L]++
						
					;;; If this angle fails, then steeper angles will also,
					;;; so break out of the angle loop and proceed to the 
					;;; next realisation...
					break
					
				ENDIF ELSE BEGIN

					PRINT,"Successful realisation angle ",angle, $
						" for run ",angleRuns[angle] + 1L, $
						" --> attempt ",angleRunsCum[angle]+1L,FORMAT='(A,I3,2(A,I4))'
					
					angleRuns[angle]++
					angleRunsCum[angle]++

					;;; Compute the glint moments

					glintFirstMoments[angle] += TOTAL(DOUBLE(glint))/DOUBLE(N)

					glintMean[angle]     += MEAN(DOUBLE(glint),/DOUBLE)
					glintVariance[angle] += VARIANCE(DOUBLE(glint),/DOUBLE)
					glintSkewness[angle] += SKEWNESS(DOUBLE(glint),/DOUBLE)

					;;; Compute the Fourier spectrum of this glint realisation
					glintSpectrum = FFT(glint,/DOUBLE)

					;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
					;;;   Compute the average glint power spectrum   ;;;
					;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

					;-- Windowing types are...
					;
					;  0 : No windowing 
					;  1 : Windowed glint Fourier spectrum
					;  2 : Windowed glint Amplitude
					;  3 : Windowed glint power spectrum
					;
					CASE powWindType OF
						0: BEGIN
							totalGlintAvgPower[*,angle] += ABS(glintSpectrum)^2.D
						END
						1: BEGIN
							totalGlintAvgPower[*,angle] += ABS(hanWindow_N*glintSpectrum)^2.D
						END
						2: BEGIN
							totalGlintAvgPower[*,angle] += (hanWindow_N*ABS(glintSpectrum))^2.D
						END
						3: BEGIN
							totalGlintAvgPower[*,angle] += hanWindow_N*(ABS(glintSpectrum)^2.D)
						END
					ENDCASE

					;;; Apply the Hanning window function to the glint Fourier 
					;;; spectrum, so that we can calculate the glint power spectrum
					;;; and bispectrum, while minimising aliasing.

					;-- Windowing types are...
					;
					;  0 : No windowing 
					;  1 : Windowed glint Fourier spectrum
					;  2 : Windowed glint Amplitude
					;  3 : Windowed glint power spectrum
					;
					CASE bispWindType OF
						0: BEGIN
							glintSpectrum = glintSpectrum
						END
						1: BEGIN
							glintSpectrum = hanWindow*glintSpectrum
						END
					ENDCASE

					;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
					;;;   Compute the average glint bispectrum   ;;;
					;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

					;;; Calculate the average bispectra (for the reduced domain)
					FOR j=0L,NN4 DO BEGIN
						FOR i=j,NN2-j DO BEGIN
							glintBispectrum[i,j,angle] += glintSpectrum[i]*glintSpectrum[j]*CONJ(glintSpectrum[i+j])
						ENDFOR
					ENDFOR

					;;; Calculate the average component power spectra (for the reduced domain)
					FOR j=0L,NN4 DO BEGIN
						FOR i=j,NN2-j DO BEGIN
							glintComponentPower[i,j,angle] += (ABS(glintSpectrum[i]*glintSpectrum[j]))^2.D
						ENDFOR
					ENDFOR

					;;; Calculate the average sum power spactra (for the reduced domain)
					FOR j=0L,NN4 DO BEGIN
						FOR i=j,NN2-j DO BEGIN
							glintSumPower[i,j,angle] += (ABS(glintSpectrum[i+j]))^2.D
						ENDFOR
					ENDFOR

				ENDELSE		;-- End checking for zero-glint of this angle

			ENDIF		;-- End checking for completion of this angle

		ENDFOR		;-- End angle loop

	ENDWHILE	;-- End realisation loop

	PRINT,""
	PRINT,"AngleRuns:    ",angleRuns," ... for total of ", $
		FIX(TOTAL(angleRuns))," / ",N_r*GEOM.N_angles;,FORMAT='(//A,5I4,A,I4,A,I4)'
	PRINT,"AngleRunsCum: ",angleRunsCum;,FORMAT='(A,5I4)'
	
	N_runs = N_r
	N_r = N_r_cum
	PRINT, "Final N_r_cum = ",N_r_cum

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

	PRINT,''
	PRINT,'IDL Elevation Moments:   ',elevMean,elevVariance,elevSkewness*(elevVariance^1.5D),FORMAT='(A,3E16.7)'
	PRINT,'    Elevation Moments:   ',elevMoments[0],elevMoments[1],elevMoments[2],FORMAT='(A,3E16.7)'
	PRINT,'    Elevation Cumulants: ',elevCumulants[0],elevCumulants[1],elevCumulants[2],FORMAT='(A,3E16.7)'
	PRINT,''
	PRINT,'IDL Slope Moments:   ',slopeMean,slopeVariance,slopeSkewness*(slopeVariance^1.5D),FORMAT='(A,3E16.7)'
	PRINT,'    Slope Moments:   ',slopeMoments[0],slopeMoments[1],slopeMoments[2],FORMAT='(A,3E16.7)'
	PRINT,'    Slope Cumulants: ',slopeCumulants[0],slopeCumulants[1],slopeCumulants[2],FORMAT='(A,3E16.7)'
	PRINT,''
	;PRINT,'IDL Glint Moments:   ',glintMean,glintVariance,glintSkewness*(glintVariance^1.5D),FORMAT='(A,3E16.7)'
	PRINT,'Glint First Moments:   ',glintFirstMoments,FORMAT='(A/,'+STRING(GEOM.N_angles)+'E16.7/)'
	PRINT,'Glint Cumulants:       ',TRANSPOSE(glintCumulants),FORMAT='(A/,3('+STRING(GEOM.N_angles)+'E16.7/))'

	PRINT,''
	PRINT,"Elevation third moment from bicovariance: ",DOUBLE(elevThirdMomentFunction[0L,0L]),FORMAT='(A,E16.7)'
	PRINT,"  Elevation third moment from bispectrum: ",TOTAL(DOUBLE(elevBispectrum)),FORMAT='(A,E16.7)'
	PRINT,''
	PRINT,"    Slope third moment from bicovariance: ",DOUBLE(slopeThirdMomentFunction[0L,0L]),FORMAT='(A,E16.7)'
	PRINT,"      Slope third moment from bispectrum: ",TOTAL(DOUBLE(slopeBispectrum)),FORMAT='(A,E16.7)'
	PRINT,''
	PRINT,"    glint third moment from bicovariance: ",DOUBLE(glintThirdMomentFunction[0L,0L,*]),FORMAT='(A,5E16.7)'
	PRINT,"      glint third moment from bispectrum: "
	FOR angle=0L, GEOM.N_angles-1L DO BEGIN
		PRINT, TOTAL(DOUBLE(glintBispectrum[*,*,angle])),FORMAT='(E16.7)'
	ENDFOR
	PRINT,''

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

	PRINT,"Open output filename: ",filename

	outFile = fileinfo(fileName)
	help, fileinfo(fileName), /structure
	PRINT, outFile

	IF (NOT outFile.EXIST) THEN BEGIN
		;;; Create and open file using SD interface
		fileID = HDF_SD_START(fileName, /CREATE)
		;fileID = HDF_OPEN(fileName, /CREATE,/WRITE)
		PRINT, 'Created new HDF file: ',fileName
	ENDIF ELSE BEGIN
		;;; Create and open file using SD interface
		PRINT, 'HDF file ',fileName,' exists, opening...'
		fileID = HDF_SD_START(fileName, /RdWr)
		;fileID = HDF_OPEN(fileName, /WRITE)
		PRINT, 'Opened HDF file ',fileName,' for reading and writing'
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

		PRINT, 'Writing geometry angles and slopes...'

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
		
		PRINT, 'Writing the elevation, slope and glint moments...'

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
		
		PRINT, 'Writing the elevation, slope and glint cumulants...'

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
		
		PRINT, 'Writing indices of nonlinear components...'

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

		PRINT, 'Writing average power spectra ...'

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
		
		PRINT, 'Writing average second moment functions ...'
		
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
		
		PRINT, 'Writing average second cumulant functions ...'
		
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
		
		PRINT, 'Writing average bispectra ...'
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
		
		PRINT, 'Writing average bicoherence ...'
		
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
		
		PRINT, 'Writing average third moment functions ...'
		
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
		
		PRINT, 'Writing average third Cumulant functions ...'
		
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

		PRINT, 'Updating the elevation, slope and glint moments...'

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
		
		PRINT, 'Updating the free and bound indices...'
		
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

		PRINT, 'Updating the elevation, slope and glint power spectra...'

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

		PRINT, 'Updating the elevation, slope and glint bispectra...'

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

		PRINT, 'Updating the elevation, slope and glint third moment functions...'

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
	PRINT, 'Write Operation Completed'
	PRINT, ''






	
END
