import pySPM
import numpy as np
from pySPM_data import get_data
import matplotlib.pyplot as plt

def test_ITA_loading():
    TOF = pySPM.ITA(get_data("BigSmiley.ita"))
    assert len(TOF.get_masses()) == 27
    assert TOF.Nscan == 100

def test_SI_image():
    TOF = pySPM.ITA(get_data("BigSmiley.ita"))
    assert TOF.img.shape == (512, 512)
    assert TOF.img.shape[0] == TOF.size['pixels']['y']
    assert TOF.img.shape[1] == TOF.size['pixels']['x']

def test_getAddedChannel():
    TOF = pySPM.ITA(get_data("BigSmiley.ita"))
    img1 = TOF.getAddedImageByMass(197)
    assert type(img1) is pySPM.SPM_image
    img2, CH = TOF.getAddedImageByName('Au')
    assert type(img2) is pySPM.SPM_image
    assert len(CH) == 1
    assert CH[0]['assign'] == 'Au-'
    assert np.all(img1.pixels == img2.pixels)


