import os

import pySPM

data = os.path.join(os.path.dirname(__file__), "AuTi_Img_Bi1_p_4_0.ita")

import unittest


class TestPCA(unittest.TestCase):
    def test_ita_collection(self):
        C = pySPM.ITA_collection(data)
        C.runPCA()
        # assert C.PCA.loadings().shape == (19, 19)


if __name__ == "__main__":
    unittest.main()
