# LSICE
Level Set Discrete Element Method for Sea Ice

LS-ICE Instructions

Compile Main_Periodic_Ocean program using:
make -f makefile_fsd;


Run LS-ICE by executing:
./Main_Periodic_Ocean #Cores #Nodes #Iteration/Seed;



Required libraries for running LS-ICE:

Eigen: https://gitlab.com/libeigen/eigen
LSMLIB: https://github.com/velexi-research/LSMLIB
OpenMP:https://www.openmp.org/specifications/
MPI: https://www.open-mpi.org/software/ompi/v5.0/
