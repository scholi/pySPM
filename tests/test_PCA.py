import pySPM
import numpy as np

import os
data = os.path.join(os.path.dirname(__file__), "AuTi_Img_Bi1_p_4_0.ita")

def test_ITAcollection():
    C = pySPM.ITA_collection(data)
    C.runPCA()
    assert C.PCA.loadings().shape == (19, 19)
