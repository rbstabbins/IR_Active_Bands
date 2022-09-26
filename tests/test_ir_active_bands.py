"""Unit tests for the IR band combination calculator tool.

Part of the Spectral Parameters Toolkit
Author: Roger Stabbins, NHM
Date: 05-09-2022
"""
import unittest
import os
from pathlib import Path
import pandas as pd
import ir_active_bands as irab

class TestIRActiveBands(unittest.TestCase):
    """Class to test ir_active_bands
    """

    def test_wvl2wvn(self):
        """Unit Test of the wavelength to wavenumber convertor
        """
        wavelength = 10 # microns
        result = irab.wvl2wvn(wavelength)
        expected = 1000 # wavenumbers
        self.assertAlmostEqual(result, expected, places=7)

    def test_wvn2wvl(self):
        """Unit Test of the wavenumber to wavelength convertor
        """
        wavenumber = 1000 # wavenumbers
        result = irab.wvn2wvl(wavenumber)
        expected = 10 # microns
        self.assertAlmostEqual(result, expected, places=7)

    def test_combine(self):
        """Unit Test of the band combination function
        """
        fundamentals = pd.DataFrame(data={'band-centre': [1, 10, 100]}, index=['v1','v2','v3'])
        result = irab.combine(fundamentals)
        expected = 0.9009009009009009 # in wavenumbers: 11100 = 10000 + 1000 + 100
        self.assertAlmostEqual(result, expected, places=7)

    def test_init_dict(self):
        """Unit Test of the IRActiveBands constructor
        """
        test_molecule = {'test_molecule': {'v1': 1.0}}
        result = irab.IRActiveBands(test_molecule)
        expected = pd.DataFrame(
            data={'band-centre': 1.0,
                    'type': 'fundamental',
                    'level': 1.0},
            index=['v1']
        )
        with self.subTest():
            pd.testing.assert_frame_equal(result.absorptions, expected)
        with self.subTest():
            pd.testing.assert_frame_equal(result.fundamentals, expected)

    def test_init_str(self):
        """Unit Test of the IRActiveBands constructor
        """
        test_molecule = 'OH'
        result = irab.IRActiveBands(test_molecule)
        expected = pd.DataFrame(
            data={'band-centre': 2.75,
                    'type': 'fundamental',
                    'level': 1.0},
            index=['v1']
        )
        with self.subTest():
            pd.testing.assert_frame_equal(result.absorptions, expected)
        with self.subTest():
            pd.testing.assert_frame_equal(result.fundamentals, expected)

    def test_compute_overtones(self):
        """Unit Test of the IRActiveBands.compute_overtones function
        """
        test_molecule = {'test_molecule': {'v1': 1.0}}
        test_absorptions = irab.IRActiveBands(test_molecule)
        result_overtones = test_absorptions.compute_overtones()
        result_absorptions = test_absorptions.absorptions

        overtone_1 = 1.0 / 2
        overtone_2 = 1.0 / 3
        expected_overtones = {'2v1': overtone_1, '3v1': overtone_2}
        expected_absorptions = pd.DataFrame(
            data={'band-centre': [1.0, overtone_1, overtone_2],
                    'type': ['fundamental', 'overtone', 'overtone'],
                    'level': [1.0, 0.9, 0.9]},
            index=['v1', '2v1', '3v1']
        )
        with self.subTest():
            self.assertDictEqual(result_overtones, expected_overtones)
        with self.subTest():
            pd.testing.assert_frame_equal(expected_absorptions, result_absorptions)

    def test_compute_combinations(self):
        """Unit Test of the IRActiveBands.compute_combinations function
        """
        fundamental = 1.0
        test_molecule = {'test_molecule': {'v1': fundamental}}
        test_absorptions = irab.IRActiveBands(test_molecule)
        test_absorptions.compute_overtones()
        test_absorptions.compute_combinations()
        result = test_absorptions.absorptions

        overtone_1 = fundamental / 2
        overtone_2 = fundamental / 3
        combi_1 = 1E4 / ((1E4/fundamental) + (1E4/overtone_1))
        combi_2 = 1E4 / ((1E4/fundamental) + (1E4/overtone_2))
        combi_3 = 1E4 / ((1E4/overtone_1) + (1E4/overtone_2))
        combi_4 = 1E4 / ((1E4/fundamental) + (1E4/overtone_1) + (1E4/overtone_2))
        expected = pd.DataFrame(
            data={'band-centre': [fundamental, overtone_1, overtone_2, combi_1, combi_2, combi_3, combi_4],
                    'type': [
                        'fundamental',
                        'overtone',
                        'overtone',
                        'fundamental+overtone',
                        'fundamental+overtone',
                        'overtone+overtone',
                        'fundamental+overtone+overtone'],
                    'level': [1.0, 0.9, 0.9, 0.7, 0.7, 0.6, 0.3]},
            index=['v1', '2v1', '3v1', 'v1+2v1', 'v1+3v1', '2v1+3v1', 'v1+2v1+3v1']
        )
        pd.testing.assert_frame_equal(expected, result)

    def test_export_absorptions(self):
        """Unit Test of the exporting of absorptions DataFrame to csv
        """
        fundamental = 1.0
        test_molecule = {'test_molecule': {'v1': fundamental}}
        test_absorptions = irab.IRActiveBands(test_molecule)
        test_absorptions.compute_overtones()
        test_absorptions.compute_combinations()
        test_path = os.path.dirname(os.path.realpath(__file__))
        test_absorptions.export_absorptions(test_path)
        expected_file = Path(test_path, 'test_molecule_absorptions').with_suffix('.csv')
        file_exists = os.path.exists(expected_file)

        self.assertEqual(file_exists, True)

        os.remove(expected_file)

if __name__ == '__main__':
    unittest.main()