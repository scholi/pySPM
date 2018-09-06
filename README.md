# pySPM
pySPM is a python library (python3, but should be compatible with python2) in order to read, handle and plot Scanning Probe Microscopy (SPM) images as well as ToF-SIMS data.

For now it support the following formats:
* Nanoscan .xml file format
* Bruker
* Iontof ToF-SIMS fileformats:
	* ITA
	* ITM
	* ITS
* Nanonis SXM file

## Important
This library is offered as it is and is still in development. Please note that reading the raw data was done by reverse engineering and guessing and not with a manual as the file format is proprietary. It seems to work well with the data used by the developper of this library, but there is **NO GUARANTY** that this library will work correctly with your own specific data.

If you find bugs and issues, please report them to the developpe: https://github.com/scholi/pySPM/issues

## News
### Structure reformatting
A setup.py is prestent in order to install the package easily. => in order to use the library do ```pip install -e . ```

### Nice Spectra Plotting
```python
import pySPM

filename = "..."
TOF = pySPM.ITA(filename)
TOF.showSpectrumAround(pySPM.utils.get_mass('C2H3NO'), pretty=True, formula=True)
```

![Spectra](../master/doc/Spectra.png)

## Dependencies
This library requires the following packages
* mendatory
    * numpy
    * scipy
    * matplotlib
* for PCA
    * scikit-learn
    * pandas
* for GUI
    * pyQT5
* displaying progressbar (while passing the prog=True parameter to functions)
    * tqdm
    
## Installation
### for regular users
```bash
python setup.py install
```

### For developpers and hackers
If you wish to adjust the library to your need, the best is to install it in editable mode as follow from the root pySPM directory:
```bash
pip install -e .
```

## Documentation
The documentation is still in its early stage
[read the documentation](../master/doc/pySPM%20Documentation.ipynb)

## Citing
If you use this library for your work, please think about citing it.
[![DOI](https://zenodo.org/badge/64209401.svg)](https://zenodo.org/badge/latestdoi/64209401)
