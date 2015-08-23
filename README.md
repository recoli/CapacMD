# CapacMD
A Molecular Dynamics Simulation Engine Using Capacitance-Polarizability Force Field
http://www.theochem.kth.se/~lixin/capacmd/index.html

Related paper:
X. Li and H. Agren, J. Phys. Chem. C, 2015, 119, 19430-19437.
http://pubs.acs.org/doi/abs/10.1021/acs.jpcc.5b04347

Usage:
mpirun -np 4 ./mpi-main -gro init.gro -par param.txt -mds mdset.txt
