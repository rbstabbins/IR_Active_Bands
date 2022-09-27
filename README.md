# IR Active Bands

A simple toolkit for computing and visualising the locations of harmonic overtone and combination
absorption features from the fundamental vibrational modes of infrared active molecules found in
minerals of interest to the exploration of planetary surfaces.

## Installation

You can install the package with:
```
pip install git+https://github.com/rbstabbins/IR_Active_Bands
```

## Usage

The tool is a Python package. Once installed, you can access it by typing:
```
from ir_active_bands import IRActiveBands
```
The package comes loaded with the fundamentals for H2O, OH and CO3. These can be accessed by typing, for example:
```
absorptions = IRActiveBands('H2O')
```
There is a notebook that provides a breakdown of the functions of the code, but in brief, you can simply compute and visualise the overtones and combinations of the given molecule with:
```
absorption_table = absorptions.compute_combinations_and_show(range = [0.9, 4.0])
```
The range here has been restricted to 0.9 - 4.0 microns. The function returns the computed absorptions in a DataFrame, and produces a plot illustrating the locations of the absorption features in a bar plot.

## Potential for Extension

This simple tool only takes the locations of fundamental vibrations, and computes up to the first 2 harmonic overtones, and pair and triplet fundamental and overtone combinations.

The code would benefit from:
* handling and computing absorption cross-sections, to indicate absorption intensities,
* handling and computing anharmonic overtones,
* including effects of bonds with other molecules.

## Further note

This is my first published Python package. It has very limited utility, but it is also very simple, so it's been useful to use as an introduction to package publication. I'm following the guidelines of the [Turing Institute Research Software Engineering with Python course](https://alan-turing-institute.github.io/rse-course/html/module06_software_projects/index.html).
