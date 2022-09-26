"""A library for generating a table of overtones and combinations
of the fundamental vibrational modes of the IR active molecules most relevant
to Martian/planetary mineralogy: H2O, OH and CO3.

Author: Roger Stabbins, NHM
Date: 02-09-2022
"""
from typing import Union
from typing import Dict
from typing import List
import itertools
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import Normalize

UNITS = 'microns'
LEVELS = {
    'fundamental': 1.0,
    'overtone': 0.9,
    'fundamental+fundamental': 0.8,
    'fundamental+overtone': 0.7,
    'overtone+overtone': 0.6,
    'fundamental+fundamental+fundamental': 0.5,
    'fundamental+fundamental+overtone': 0.4,
    'fundamental+overtone+overtone': 0.3,
    'overtone+overtone+overtone': 0.2
}
LEVELS_SHRT = {
    'fundamental': 1.0,
    'overtone': 0.9,
    'f + f': 0.8,
    'f + o': 0.7,
    'o + o': 0.6,
    'f + f + f': 0.5,
    'f + f + o': 0.4,
    'f + o + o': 0.3,
    'o + o + o': 0.2
}
FUNDAMENTALS = {
    # Reference:
    #   Hunt, "Spectral signatures of particulate minerals in the visible and near infrared",
    #   Geophysics, 42, 3, pp. 501-511, April 1977,
    #   DOI:10.1190/1.1440721
    'H2O': {'v1': 3.106, 'v2': 6.08, 'v3': 2.903},
    'OH': {'v1': 2.75},
    'CO3': {'v1': 9.23, 'v2': 11.36, 'v3': 7.0, 'v4': 14.0}
}

class IRActiveBands:
    """A class for computing Infrared Active Band overtones and combinations from the fundamental
    absorptions of a given IR active molecule.
    """

    def __init__(self, molecule: Union[str, Dict]) -> None:
        """Initiate the class object by storing the given fundamentals in a DataFrame.

        :param fundamentals_dict: fundamental absorption band-centers of given material
        :type fundamentals_dict: Dict
        """
        if isinstance(molecule, str):
            assert molecule in FUNDAMENTALS.keys(), f'Molecule is not in dictionary. \
                Please add new molecule fundamentals to dictionary, \
                or use one of the following molecules: {FUNDAMENTALS.keys()}'
            fundamentals_dict = FUNDAMENTALS[molecule]
        else:
            assert isinstance(molecule, dict), 'Please provide either a molecule or \
                fundamentals dictionary'
            fundamentals_dict = list(molecule.values())[0]
            molecule = list(molecule.keys())[0]
        self.molecule = molecule
        self.fundamentals = pd.DataFrame(
            data=fundamentals_dict.values(),
            index=fundamentals_dict.keys(),
            columns=['band-centre'])
        self.fundamentals['type'] = 'fundamental'
        self.fundamentals['level'] = LEVELS['fundamental']
        self.absorptions = self.fundamentals
        self.overtones = None # populated in the compute_overtones() method

    def compute_combinations_and_show(self, range: List = None) -> pd.DataFrame:
        """Top level function for comuting the overtones and combinations, collecting these,
        and displaying the complete set of absorptions within the given range.

        :param range: restriction on spectral window in microns, optional
        :type range: List
        :return: table of fundamentals, overtones and combinations of absorption features.
        :rtype: pd.DataFrame
        """
        print("Computing overtones...")
        self.compute_overtones()
        print("Computing combinations...")
        self.compute_combinations()
        if range is not None:
            print("Dropping out-of-range absorptions...")
            self.filter_absorptions(range)
        print("Visualising absorptions...")
        self.visualise_absorptions()
        return self.absorptions

    def compute_overtones(self) -> Dict:
        """Compute overtone frequencies and add these to the absorptions DataFrame.
        Also returns the overtones in a dictionary.

        :param fundamentals: fundamental absorptions in wavelength to compute overtones for
        :type fundamentals: np.array
        :return: overtones in wavelength
        :rtype: np.array
        """
        fundamentals_dict = self.fundamentals['band-centre'].to_dict()
        overtones_dict = {}
        for identifier, wavelength in fundamentals_dict.items():
            for k in [2,3]:
                overtone_key = str(k) + identifier
                overtones_dict[overtone_key] = wavelength / k # divide by factor
        self.overtones = pd.DataFrame(
            data=overtones_dict.values(),
            index=overtones_dict.keys(),
            columns=['band-centre'])
        self.overtones['type'] = 'overtone'
        self.overtones['level'] = LEVELS['overtone']
        self.absorptions = pd.concat([self.absorptions, self.overtones])
        return overtones_dict

    def compute_combinations(self) -> Dict:
        """Compute the combinations afforded by the fundamentals and overtones, and add
        these to the DataFrame.
        Also returns the combinations in a dictionary

        :param absorptions: fundamental and overtone frequencies of molecule
        :type fundamentals: Dict
        :return: combinations absorption features
        :rtype: Dict
        """
        combinations = {}
        # check that the combinations haven't already been computed by checking that the only
        # 'types' are fundamentals and overtones:
        if not np.array_equal(['fundamental', 'overtone'], self.absorptions.type.unique()):
            raise ValueError('Band combinations have likely already been computed, \
                            as types other than "fundamental" and "overtone" are present.')
        # make combination list
        pairs = [[a,b] for a,b in itertools.combinations(self.absorptions.index, 2)]
        triplets = [[a,b,c] for a,b,c in itertools.combinations(self.absorptions.index, 3)]
        # iterate over combinations, and compute
        for pair in pairs:
            # get the absorptions given by the keys
            features = self.absorptions.loc[pair]
            band_centre = combine(features)
            label = '+'.join(pair)
            combination_type = '+'.join(features['type'])
            level = LEVELS[combination_type]
            combinations[id] = band_centre
            absorption = pd.DataFrame(data={
                    'band-centre': band_centre,
                    'type': combination_type,
                    'level': level},
                    index=[label])
            self.absorptions = pd.concat([self.absorptions, absorption])

        for triplet in triplets:
            # get the absorptions given by the keys
            features = self.absorptions.loc[triplet]
            band_centre = combine(features)
            label = '+'.join(triplet)
            combination_type = '+'.join(features['type'])
            level = LEVELS[combination_type]
            combinations[label] = band_centre
            absorption = pd.DataFrame(data={
                    'band-centre': band_centre,
                    'type': combination_type,
                    'level': level},
                    index=[label])
            self.absorptions = pd.concat([self.absorptions, absorption])
        # return dictionary
        return combinations

    def filter_absorptions(self, spectral_range: List) -> List:
        """Remove absorption features outside a range of interest

        :param absorptions: all absorptions including fundamentals, overtones and combinations
        :type absorptions: Dict
        :param spectral_range: wavelength range to consider
        :type spectral_range: List
        :return: absorption features that have been removed
        :rtype: List
        """
        dropped_features = self.absorptions[
            (self.absorptions['band-centre'] < spectral_range[0]) &
                    (self.absorptions['band-centre'] > spectral_range[1])
            ].index
        self.absorptions = self.absorptions[
            (self.absorptions['band-centre'] >= spectral_range[0]) &
                    (self.absorptions['band-centre'] <= spectral_range[1])]
        return dropped_features

    def export_absorptions(self, path: str) -> None:
        """Export the absorptions DataFrame to an Excel file at the path location

        :param file_name: name of the file
        :type file_name: str
        :param path: directory to write the file to
        :type path: str
        """
        file_name = self.molecule+'_absorptions'
        file_out = Path(path, file_name).with_suffix('.csv')
        self.absorptions.to_csv(file_out)

    def visualise_absorptions(self, **path: str) -> None:
        """Visualise the absorption features in a figure
        """
        band_centres = self.absorptions['band-centre']
        levels = self.absorptions['level']

        # converting levels to colours
        tab_cm = cm.get_cmap('tab10_r')
        norm = Normalize(vmin=0.1, vmax=1.0)
        colors = tab_cm(norm(levels))

        # make plot
        CM = 1/2.54
        fig, ax = plt.subplots(figsize=[15*CM, 10*CM])
        ax.bar(band_centres, height = levels, width = band_centres/150, color = colors)

        # add grid
        ax.grid(which='major', axis='y', lw=0.8)
        ax.grid(which='major', axis='x', lw=0.8)
        plt.minorticks_on()
        ax.grid(which='minor', axis='x', lw=0.5)

        # relabel y-axis with Levels
        ax.set_yticks(list(LEVELS_SHRT.values()))
        ax.set_yticklabels(list(LEVELS_SHRT.keys()))
        ax.set_axisbelow(True)

        # add title
        ax.set_xlabel(f'Wavelength ({UNITS})')
        ax.set_ylabel('Combination Type')
        ax.set_title(self.molecule+' Absorption Features')

        plt.yticks(rotation=45)
        plt.tight_layout()
        plt.show()
        print('done')

    def __repr__(self):
        """Print the absorptions dataframe
        """
        return self.absorptions.to_string()

    def __str__(self):
        """Print the absorptions dataframe
        """
        return self.absorptions.to_string()

def wvl2wvn(wavelength: Union[float, np.array]) -> Union[float, np.array]:
    """Convert from wavelength (microns) to wavenumber (cm^-1)

    :param wavelength: wavelength in microns
    :type wavelength: float
    :return: wavenumber in cm^-1
    :rtype: float
    """
    wavenumber = (1E6 * 1E-2) / wavelength
    return  wavenumber

def wvn2wvl(wavenumber: Union[float, np.array]) -> Union[float, np.array]:
    """Convert from wavenumber (cm^-1) to wavelength (microns)

    :param wavenumber: wavenumber in cm^-1
    :type wavenumber: float
    :return: wavelength in microns
    :rtype: float
    """
    wavelength = (1E6 * 1E-2) / wavenumber
    return wavelength

def combine(features: pd.DataFrame) -> float:
    """Combine the given absorptions frequencies

    :param absorptions: absorptions in wavelength to combine
    :type fundamentals: Tuple[]
    :return: combination frequency in wavelength
    :rtype: float
    """
    absorptions = wvl2wvn(features['band-centre'].to_numpy())
    combination = np.sum(absorptions)
    combination = wvn2wvl(combination)
    return combination

if __name__ == "__main__":

    h2o_absorptions = IRActiveBands('OH')
    print(h2o_absorptions)
    h2o_overtones = h2o_absorptions.compute_overtones()
    print(h2o_overtones)
    h2o_combinations = h2o_absorptions.compute_combinations()
    print(h2o_combinations)
    h2o_dropped_features = h2o_absorptions.filter_absorptions([0.9, 4.0])
    print(h2o_dropped_features)
    h2o_absorptions.visualise_absorptions()