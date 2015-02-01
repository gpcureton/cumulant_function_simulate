# Sunglint Cumulant Function Simulation

This program simulates multiple elevation and slope realisations with
given spectra, introducing quadratic phase coupling at selected wavenumbers.
Glint realisations are computed from the slope realisations. The average 
elevation, slope and glint moments and cumulants, and moment and cumulant 
functions are computed. The results are written to a HDF5 file.
    
## Dependencies

* pytables
* h5py
* numpy
* scipy

##Usage:

```[bash]
usage: cumulantFunctionSimulate.py [-h] [-n N] [-N NN] [-d DELTA_X]
                                   [-r NUMREALS] [-S SPECTRUMTYPE] [-l]
                                   [-o OUTFILE] [-v]

optional arguments:
  -h, --help            show this help message and exit
  -n N, --numdatapoints N
                        Number of points of the glint dataset. Must be an
                        integer power of 2. [default: 1024]
  -N NN, --num2Dpoints NN
                        Number of points of the bispectrum array side. Must be
                        an integer power of 2. [default: 64]
  -d DELTA_X, --deltax DELTA_X
                        The spatial increment in meters. [default: 0.02]
  -r NUMREALS, --num_realisations NUMREALS
                        Number of realisations. [default: 100]
  -S SPECTRUMTYPE, --spectrum_type SPECTRUMTYPE
                        Form of the elevation power spectrum. Possible values
                        are... 'phillips_3', 'phillips_4', 'gaussian'.
                        [default: 'phillips_3']
  -l, --nonlinear_modes
                        Switch on nonlinear modes. [default: False]
  -o OUTFILE, --output_file OUTFILE
                        The full path of the output HDF5 file. [default:
                        'outSimulatedGlint.h5']
  -v, --verbose         each occurrence increases verbosity 1 level from
                        ERROR: -v=WARNING -vv=INFO -vvv=DEBUG
```