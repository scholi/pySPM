import os

import numpy as np

import pySPM

data = os.path.join(os.path.dirname(__file__), "AuTi_Img_Bi1_p_4_0.ita")

import unittest


class TestITA(unittest.TestCase):
    def test_ITA_loading(self):
        pySPM.ITA(data)

    def test_SI_image(self):
        TOF = pySPM.ITA(data)
        assert TOF.img.pixels.shape[0] == TOF.size["pixels"]["y"]
        assert TOF.img.pixels.shape[1] == TOF.size["pixels"]["x"]

    def test_getAddedChannel(self):
        TOF = pySPM.ITA(data)
        img1 = TOF.getAddedImageByMass(107)
        img2, CH = TOF.getAddedImageByName("Ag", strict=True)
        assert len(CH) == 1
        assert CH[0]["assign"] == "Ag+"
        assert np.all(img1.pixels == img2.pixels)


if __name__ == "__main__":
    unittest.main()
