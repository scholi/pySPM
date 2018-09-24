import pySPM
import numpy as np

def test_version():
    assert pySPM.__version__=='0.2.6'


def test_SPM_image():
    I = pySPM.SPM_image(BIN=np.random.randint(0,255,(256,256)))
