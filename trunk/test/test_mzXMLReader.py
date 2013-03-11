
import unittest
import os

class TestmzXMLReader(unittest.TestCase):

    def setUp(self):
        dirname = os.path.dirname(os.path.abspath(__file__))
        joined = os.path.join(dirname, "data")
        # the test file is AFA1451A160MA2 of the dataset PAe001446_mzXML_201102101633.tar.gz from Peptide Atlas
        self.filename = os.path.join(joined, "testfile.small.mzXML")

    def test_readfile(self):
        import msproteomicstoolslib.format.mzXMLreader as mzXMLReader
        reader = mzXMLReader.mzXMLReader(self.filename, True)

    def test_readscan(self):
        import msproteomicstoolslib.format.mzXMLreader as mzXMLReader
        reader = mzXMLReader.mzXMLReader(self.filename, True)
        scan = reader.read_scan(5)

    def test_readpeaks(self):
        import msproteomicstoolslib.format.mzXMLreader as mzXMLReader
        reader = mzXMLReader.mzXMLReader(self.filename, True)
        scan = reader.read_scan(5, True)
        self.assertEqual( len(scan.peaks), 615)
        self.assertAlmostEqual( scan.peaks[0].int, 271.413513184)
        self.assertAlmostEqual( scan.peaks[0].mz, 350.3074646)
        self.assertAlmostEqual( scan.max_peak().int, 10638.6044922)
        self.assertAlmostEqual( scan.max_peak().mz, 390.929901123)

if __name__ == '__main__':
    unittest.main()
