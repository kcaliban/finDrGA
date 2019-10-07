# Dvelopr

Dvelopr is a distributed computing application applying the principle of
[genetic algorithms](https://en.wikipedia.org/wiki/Genetic_algorithm)
to find a peptide binder given a receptor.

## Description

Dvelopr performs molecular dynamics simulations using GROMACS
and calculates binding affinities using AutoDock Vina. The results serve
as the fitness of a possible binder and therefore as the basis for
application of the principles of Darwinian evolution.

Peptides with higher fitness have a higher probability of being recombined
into new peptides for the next generation. After as many individuals as
in the initial population have been gathered through recombination and a
customizable amount of duplication, each peptide is mutated at a
random point with a user-defined probability.

After recombination, duplication and mutation, the newly generated individuals
are simulated and binding affinities calculated again.
This process is repeated for a user-defined amount of generations.

## Dependencies

* An implementation of the MPI standard, for distribution on cluster nodes, but also required to run on a single computer. We used [MPICH](https://www.mpich.org/)
* [AutoDock Vina](http://vina.scripps.edu/), for binding affinity calculations
* [GROMACS](http://www.gromacs.org/), for molecular dynamics
* [MGLTools](http://mgltools.scripps.edu/), for generation of files required by AutoDock Vina
* [PyMOL (Open Source)](https://sourceforge.net/projects/pymol/), for generation of PDB files

## Installation

Dvelopr is written for Linux. Make sure you have all dependencies installed, then
move on to the following steps:

Go to any directory you like and clone this repository, change into its directory
and compile using make
```bash
cd ilikethisdirectory
git clone https://github.com/kcaliban/Dvelopr.git
cd Dvelopr
make
```

Before you can use Dvelopr you have to configure it. Take a look at `config.ini`
and change the settings accordingly, making sure all directories you
specify exist.

## Usage

### Arguments

Dvelopr takes four command-line arguments.
* -n : Number of generations
* -m : Size of population
* -p : Probability of random-point mutation for each individual
* -c : Percentage of previous generation to copy, as a floating-point number

### Single computer

Dvelopr is written for computer clusters, it can however be executed on a single
computer.

```bash
mpirun -np 1 Dvelopr -n 100 -m 50 -p 0.5 -c 0.2 : -np 1 PoolWorker
```

### Computer cluster

For computation on a computing cluster, you have to specify how many
individual computing nodes (not threads!) you can use:

```bash
mpirun -np 1 Dvelopr -n 100 -m 50 -p 0.5 -c 0.2 : -np NUMNODES PoolWorker
```

## License

See LICENSE file

Used libraries:
* inih is written by Ben Hoyt, see src/inih/LICENSE
* cxxopts is written by Jarryd Beck, see src/cxxopts/LICENSE
