#!/bin/bash -x
#  File: <name>
#  Author: J. G. Z. de Mattos <joao.gerd@inpe.br>
#  Created: qui 09 jun 2022 10:34:17
#  Last Update: qui 09 jun 2022 10:34:17
#  Notes:
rm -fr *.mod *.o

INCLUDE='-I/opt/libmisc/include -I/opt/w3lib-2.0.6/include -I/opt/hdf5-hdf5-1_14_1-2/include -I/opt/netcdf-c-4.9.2/include -I/opt/netcdf-fortran-4.6.0/include'
LDFLAGS='-L/opt/libmisc/lib -lmisc -L/opt/w3lib-2.0.6/lib -lw3 -L/opt/netcdf-fortran-4.6.0/lib -lnetcdff -L/opt/netcdf-c-4.9.2/lib -lnetcdf -lnetcdf -lm'
DFLAGS= #'-ggdb3 -O0 -fbacktrace -fcheck=bounds -fcheck=all -Wuninitialized -ffpe-trap=zero,invalid,overflow,underflow'
FLAGS="$DFLAGS "
gfortran -c ${FLAGS} coord_compute.f90 $INCLUDE
gfortran -c ${FLAGS} accessGrib.F90 $INCLUDE
gfortran -c ${FLAGS} accessNetcdf.f90 $INCLUDE
gfortran -c ${FLAGS} m_GrADSfiles.F90 $INCLUDE
gfortran -c ${FLAGS} fileAccess.f90 $INCLUDE

gfortran -o testeFA.x testeFA.f90 coord_compute.o accessGrib.o accessNetcdf.o m_GrADSfiles.o fileAccess.o ${LDFLAGS}
