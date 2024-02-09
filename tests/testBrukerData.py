import unittest
import pySPM
import os

path = os.path.abspath(os.path.join(__file__, "../../data"))

def get_Bruker_size(file, channel="Height Sensor"):
    return pySPM.Bruker(os.path.join(path, file)).get_channel(channel, mock_data=True).size

class TestBrukerData(unittest.TestCase):


    def test_data1(self):
        s = get_Bruker_size("bruker_data_header.001", "Height")
        assert s['pixels']['x']==512
        assert s['pixels']['y']==32
        assert s['real']['x']==65.0
        assert s['real']['y']>=4 and s['real']['y']<4.1

    def test_data3(self):
        s = get_Bruker_size("bruker_data_header.003")
        assert s['pixels']['x']==256
        assert s['pixels']['y']==256
        assert s['real']['x']==2
        assert s['real']['y']==2

    def test_data2(self):
        s = get_Bruker_size("bruker_data_header.002")
        assert s['pixels']['x']==512
        assert s['pixels']['y']==358
        assert s['real']['x']==300.0
        assert s['real']['y']>209 and s['real']['y']<210

    def test_data0(self):
        s = get_Bruker_size("bruker_data_header.000")
        assert s['pixels']['x']==512
        assert s['pixels']['y']==512
        assert s['real']['x']==93.0
        assert s['real']['y']==93

if __name__=="__main__":
    unittest.main()