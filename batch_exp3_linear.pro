;;;

;;; Arguments are...

;;; Data length
;;; Tricovariance array size
;;; spatial increment in meters
;;; Realisations in simulation
;;; Power Spectrum type
;;; Power spectrum exponent
;;; glint power spectrum window
;;; glint bispectrum window
;;; Nonlinear coupling off (0) or on (1) in simulated data

;-- Exponent Three - Linear

cumulantFunctionSimulate,4096,512,0.02,20000,'phillips',3,0,0,0
cumulantFunctionSimulate,4096,512,0.02,20000,'phillips',3,0,1,0
cumulantFunctionSimulate,4096,512,0.02,20000,'phillips',3,1,0,0
cumulantFunctionSimulate,4096,512,0.02,20000,'phillips',3,1,1,0
cumulantFunctionSimulate,4096,512,0.02,20000,'phillips',3,2,0,0
cumulantFunctionSimulate,4096,512,0.02,20000,'phillips',3,2,1,0
cumulantFunctionSimulate,4096,512,0.02,20000,'phillips',3,3,0,0
cumulantFunctionSimulate,4096,512,0.02,20000,'phillips',3,3,1,0

cumulantFunctionSimulate,4096,512,0.04,20000,'phillips',3,0,0,0
cumulantFunctionSimulate,4096,512,0.04,20000,'phillips',3,0,1,0
cumulantFunctionSimulate,4096,512,0.04,20000,'phillips',3,1,0,0
cumulantFunctionSimulate,4096,512,0.04,20000,'phillips',3,1,1,0
cumulantFunctionSimulate,4096,512,0.04,20000,'phillips',3,2,0,0
cumulantFunctionSimulate,4096,512,0.04,20000,'phillips',3,2,1,0
cumulantFunctionSimulate,4096,512,0.04,20000,'phillips',3,3,0,0
cumulantFunctionSimulate,4096,512,0.04,20000,'phillips',3,3,1,0
