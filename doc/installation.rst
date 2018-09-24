An installation guide for pySPM
===============================

Installing python
-----------------

There is several way to install python. If you are an advanced user we recommend you to install python 3.7 or newer from https://www.python.org then to install the required packages.

If you are a newcomer a quick solution to install python with all the required libraries is to install anaconda https://www.anaconda.com

Required libraries
------------------
Before installing pySPM, few libraries should be installed first.
The easiest way of doing it is to use pip. Go in a terminal (on windows hit WIN+R, then type "cmd" without quotes).

Then run:
pip install numpy scipy scikit-image sklearn tqdm matplotlib pyQt5 pandas seaborn

If you are a windows user, please download the latest numpy package for your version here: https://www.lfd.uci.edu/~gohlke/pythonlibs/#numpy
Then in the terminal, go in your download folder (usually "cd Downloads"), then type:

pip install -U "numpy-1.15.1+mkl-cp37-cp37m-win_amd64.whl"

in order to install the latest numpy package (please adapt the name with the verison you downloaded).

Installing the pySPM library
----------------------------

You can find and download the library here: https://github.com/scholi/pySPM
https://github.com/scholi/pySPM/archive/master.zip

Then you can install it with pip:

pip install "pySPM-master.zip"

Update the pySPM library
------------------------

You can find and download the library here: https://github.com/scholi/pySPM
https://github.com/scholi/pySPM/archive/master.zip

Then you can upgrade it with pip:

pip install --upgrade "pySPM-master.zip"
