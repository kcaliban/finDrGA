# PepGA

PepGA is a pipeline applying the principle of [genetic algorithms](http://mgltools.scripps.edu/)
to try to find a peptide ligand given a receptor, multiple receptors,
multiple conformations of a receptor or multiple conformations of multiple receptors.
Starting out with a sample
picked from helices extracted from [RCSB](http://mgltools.scripps.edu/),
PepGA prepares and executes a molecular dynamics simulation using GROMACS.
The result is clustered, the biggest cluster extracted and subsequently
docked against the input receptors. After selection, recombination and
mutation using the docking results, the process is repeated with the new
resulting generation of peptides.

## Dependencies
* [AutoDock Vina](http://vina.scripps.edu/), for binding affinity calculations
* [GROMACS](http://www.gromacs.org/), for molecular dynamics
* [MGLTools](http://mgltools.scripps.edu/), for generation of files required by AutoDock Vina
* [PyMOL (Open Source)](http://mgltools.scripps.edu/)

## Installation

## Usage

Logging
Config
