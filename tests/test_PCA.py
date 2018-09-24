import pySPM
import numpy as np
from pySPM_data import get_data

def test_ITAcollection():
    C = pySPM.ITA_collection(get_data("BigSmiley.ita"))
    C.runPCA()
    assert C.PCA.loadings().shape == (25, 25)
