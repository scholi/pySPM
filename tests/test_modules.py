import pySPM
import numpy as np

import unittest

class TestModules(unittest.TestCase):
    def test_SPM_image(self):
        I = pySPM.SPM_image(BIN=np.random.randint(0,255,(256,256)))

if __name__ == "__main__":
    unittest.main()