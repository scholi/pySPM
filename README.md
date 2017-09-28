# pySPM
pySPM is a python library (python3, but should be compatible with python2) in order to read, handle and plot Scanning Probe Microscopy (SPM) images as well as ToF-SIMS data.

For now it support the following formats:
* Nanoscan .xml file format
* Bruker
* Iontof ToF-SIMS fileformats:
	* ITA
	* ITM
	* ITS

## Important
This library is offered as it is and is still in development. Please note that reading the raw data was done by reverse engineering and guessing and not with a manual as the file format is proprietary. It seems to work well with the data used by the developper of this library, but there is **NO GUARANTY** that this library will work correctly with your own specific data.

If you find bugs and issues, please report them to the developpe: https://github.com/scholi/pySPM/issues

## Documentation
The documentation is still in its early stage
[read the documentation](../master/doc/pySPM%20Documentation.ipynb)

## Citing
If you use this library for your work, please think about citing it.
[![DOI](https://zenodo.org/badge/64209401.svg)](https://zenodo.org/badge/latestdoi/64209401)
