
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Routine to return the glint cumulants from the glint moment(s), for 
;;; a number of geometries
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO glintCumulantsFromMoments,glintMoments,glintCumulants

	PRINT,'N_ELEMENTS: ',N_ELEMENTS(glintCumulants[0,*])
	FOR geometry=0L,N_ELEMENTS(glintCumulants[0,*])-1L DO BEGIN
		glintCumulants[0,geometry] = glintMoments[geometry]
		glintCumulants[1,geometry] = glintMoments[geometry]-glintMoments[geometry]^2.D
		glintCumulants[2,geometry] = glintMoments[geometry] $
			-3.D*glintMoments[geometry]^2.D $
			+2.D*glintMoments[geometry]^3.D
	ENDFOR
;	glintCumulants[*,0] = glintMoments[*]
;	glintCumulants[*,1] = glintMoments[*]-glintMoments[*]^2.D
;	glintCumulants[*,2] = glintMoments[*]-3.D*glintMoments[*]*glintMoments[*]+2.D*glintMoments[*]^3.D
END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Routine to return the cumulants from the moments
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO cumulantsFromMoments,dataMoments,dataCumulants

	dataCumulants[0] = dataMoments[0]
	dataCumulants[1] = dataMoments[1]-dataMoments[0]^2.D
	dataCumulants[2] = dataMoments[2]-3.D*dataMoments[0]*dataMoments[1]+2.D*dataMoments[0]^3.D

END

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

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Routine to use the symmetry properties of the bispectrum of a real 
;;; sequence to populate an entire NxN bicoherence array from the primary 
;;; octant. Takes as input an NxN real array, and the array size N, and 
;;; returns the fully populated array in the input array
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO bicoherenceSymmetry,bicoherence,N

	FOR j=0L,N/4L DO BEGIN
		FOR i=j,(N/2L-j) DO BEGIN
			bicoherence[(j LT 0L) ? N+(j) : j , (i LT 0L) ? N+(i) : i] = bicoherence[i,j]
			bicoherence[(j LT 0L) ? N+(j) : j , (-i-j LT 0L) ? N+(-i-j) : -i-j] = bicoherence[i,j]
			bicoherence[(-i-j LT 0L) ? N+(-i-j) : -i-j , (j LT 0L) ? N+(j) : j] = bicoherence[i,j]
			bicoherence[(-i-j LT 0L) ? N+(-i-j) : -i-j , (i LT 0L) ? N+(i) : i] = bicoherence[i,j]
			bicoherence[(i LT 0L) ? N+(i) : i , (-i-j LT 0L) ? N+(-i-j) : -i-j] = bicoherence[i,j]

			bicoherence[(-i LT 0L) ? N+(-i) : -i , (-j LT 0L) ? N+(-j) : -j   ] = bicoherence[i,j]
			bicoherence[(-j LT 0L) ? N+(-j) : -j , (-i LT 0L) ? N+(-i) : -i   ] = bicoherence[i,j]
			bicoherence[(-j LT 0L) ? N+(-j) : -j , (i+j LT 0L) ? N+(i+j) : i+j] = bicoherence[i,j]
			bicoherence[(i+j LT 0L) ? N+(i+j) : i+j , (-j LT 0L) ? N+(-j) : -j] = bicoherence[i,j]
			bicoherence[(i+j LT 0L) ? N+(i+j) : i+j , (-i LT 0L) ? N+(-i) : -i] = bicoherence[i,j]
			bicoherence[(-i LT 0L) ? N+(-i) : -i , (i+j LT 0L) ? N+(i+j) : i+j] = bicoherence[i,j]
		ENDFOR
	ENDFOR
END
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
