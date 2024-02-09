import os
import unittest

import pySPM

data_path = os.path.join(os.path.dirname(__file__), "data", "AuTi_Img_Bi1_p_4_0.ita")


class TestPCA(unittest.TestCase):
    def test_ita_collection(self):
        C = pySPM.ITA_collection(data_path)
        C.runPCA()
        # assert C.PCA.loadings().shape == (19, 19)


if __name__ == "__main__":
    unittest.main()
