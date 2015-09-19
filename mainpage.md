CapacMD {#mainpage}
---------

### Introduction
A Molecular Dynamics Simulation Engine Using Capacitance-Polarizability Force Field

Original paper:
[X. Li and H. &Aring;gren, J. Phys. Chem. C, 2015, 119, 19430-19437](http://pubs.acs.org/doi/abs/10.1021/acs.jpcc.5b04347)

![TOC](http://www.theochem.kth.se/~lixin/capacmd/toc-400.png)

### Download
Current version: [CapacMD-current.zip](http://www.theochem.kth.se/~lixin/capacmd/CapacMD-current.zip)

### Usage
`mpirun -np 4 ./mpi-main -gro init.gro -par param.txt -mds mdset.txt`

### Examples
1. An icosahedral gold nanoparticle in water.
2. Two gold nanoparticles in an external electric field.
