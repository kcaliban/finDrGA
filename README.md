# PepGA

PepGA is a pipeline applying the principle of [genetic algorithms](https://en.wikipedia.org/wiki/Genetic_algorithm)
to try to find a peptide ligand given a receptor, multiple receptors,
multiple conformations of a receptor or multiple conformations of multiple receptors.
Starting out with a sample
picked from helices extracted from [RCSB](https://www.rcsb.org/),
PepGA prepares and executes a molecular dynamics simulation using GROMACS.
The result is clustered, the biggest cluster extracted and subsequently
docked against the input receptors. After selection, recombination and
mutation using the docking results, the process is repeated with the new
resulting generation of peptides.

## Dependencies

* [AutoDock Vina](http://vina.scripps.edu/), for binding affinity calculations
* [GROMACS](http://www.gromacs.org/), for molecular dynamics
* [MGLTools](http://mgltools.scripps.edu/), for generation of files required by AutoDock Vina
* [PyMOL (Open Source)](https://sourceforge.net/projects/pymol/), for generation of PDB files

## Installation

PepGA is written for Linux. Make sure you have all dependencies installed, then
move on to the following steps:

Go to any directory you like and clone this repository, change into its directory
and compile using make
```bash
cd ilikethisdirectory
git clone https://github.com/kcaliban/PepGA.git
cd PepGA
make
```

Before you can use PepGA you have to configure it. Take a look at `config.ini`
and change the settings accordingly, making sure all directories you
specify exist.

## Usage

After putting all your receptors into the directory you specified in
`config.ini`, simply change into the folder of PepGA, pray to the gods
of probability and enter:
```bash
./PepGA n p1 pc m
```
Where:
* `n` is the number of generations as an integer
* `p1` is the probability of a mutation at a random uniformly picked sequence for each individual
* `pc` is the percentage of best-performing individuals to be copied at each generation as a floating point number
* `m` is the number of individuals to pick randomly from the initialpdbs folder specified in `config.ini`, if applicable

### Initial PDBs: Helices

In the psmalpha folder you will find an example of settings, multiple conformations
of a receptor (Phenol-soluble modulin Alpha 3 in D-form) and a tarball containing
all naturally occuring alpha helices of size 12 extracted from RCSB.

## License

See LICENSE file

inih is written by Ben Hoyt, see src/inih/LICENSE
