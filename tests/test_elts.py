import unittest

import pySPM


class TestElts(unittest.TestCase):
    def test_basic(self):
        assert pySPM.utils.simplify_formula("H^29SiH") == "H2^29Si"
        assert pySPM.utils.get_mass("C") == 12
        assert pySPM.utils.is_main_isotope("C", 12)
        assert not pySPM.utils.is_main_isotope("C", 13)


if __name__ == "__main__":
    unittest.main()
