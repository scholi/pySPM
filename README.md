[![Downloads](https://pepy.tech/badge/pyspm)](https://pepy.tech/project/pyspm)
[![Build](https://travis-ci.org/scholi/pySPM.svg?branch=master)](https://travis-ci.org/scholi/pySPM)
[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![Python 2.7](https://img.shields.io/badge/python-2.7-yellow.svg)](https://www.python.org/downloads/release/python-2715/)
[![Python 3.6](https://img.shields.io/badge/python-3.4+-orange.svg)](https://www.python.org/download/releases/3.4.0/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.998575.svg)](https://doi.org/10.5281/zenodo.998575)

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
### pySPM is now availabe on pypi
The installation for end-user is now very easy. If you have [pip](https://pypi.org/project/pip/), then you can install the library in a single line:

```bash
pip install pySPM
```

and to upgrade an old version:

```bash
pip install -U pySPM
```

### Toy dataset
As the data are big and not necessary for the library another package [pySPM_data](https://github.com/scholi/pySPM_data) was created with several AFM and ToF-SIMS data.

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

### Python 2.7 compatible
The library is now compatibe with Python 3 and Python 2.7

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
#### With pip (easiest)
Just open a terminal (on Windows hit key `[WINDOWS]+R`, then type cmd, then
`[ENTER]`)
```bash
pip install pySPM
```

#### By manual installing
Download the library (zip) or git file. Unzip it and run
```bash
python setup.py install
```

#### For developpers and hackers
If you wish to adjust the library to your need, the best is to install it in editable mode as follow from the root pySPM directory:
```bash
pip install -e .
```

## Documentation
The documentation is still in its early stage
[read the documentation](https://nbviewer.jupyter.org/github/scholi/pySPM/blob/master/doc/pySPM%20Documentation.ipynb)

## Citing
If you use this library for your work, please think about citing it.
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.998575.svg)](https://doi.org/10.5281/zenodo.998575)
