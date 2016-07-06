# -----------------------------------------------------------------------------------------------------------
# Intel Compiler and MKL Libraries
# -----------------------------------------------------------------------------------------------------------
MKL_PATH               = $(HOME)/intel/mkl
IFORT_PATH             = $(HOME)/intel/composerxe

# -----------------------------------------------------------------------------------------------------------
# CEQ Library (by S. Pope, http://eccentric.mae.cornell.edu/~pope/CEQ/)
# Needed only for thermodynamic equilibrium calculations
# -----------------------------------------------------------------------------------------------------------
CEQ_LIBS = /software/chimica2/libraries/ceq/libceq.a

# -----------------------------------------------------------------------------------------------------------
# OpenBLAS (to be used if Intel MKL are not available)
# -----------------------------------------------------------------------------------------------------------
OPENBLAS_INCLUDE        = -I/software/chimica2/libraries/openblas/openblas-0.2.15
OPENBLAS_INCLUDE_LIB    = -L/software/chimica2/libraries/openblas/openblas-0.2.15
OPENBLAS_LIBS_INTEL     = -lopenblas
OPENBLAS_LIBS_GNU       = -lopenblas

# -----------------------------------------------------------------------------------------------------------
# LINPACK (needed by the following packages: DASPK, MEBDF, ODEPACK, RADAU)
# -----------------------------------------------------------------------------------------------------------
LINPACK_INCLUDE_LIB    = -L/software/chimica2/libraries/linpack/linpack-gcc-4.8.2/lib
LINPACK_LIBS_INTEL     = -lLINPACK_intel_64bit
LINPACK_LIBS_GNU       = -lLINPACK_gcc_64bit

# -----------------------------------------------------------------------------------------------------------
# SUNDIALS
# -----------------------------------------------------------------------------------------------------------
SUNDIALS_INCLUDE_GNU       = -I/software/chimica2/libraries/sundials/sundials-2.6.2-opensmoke++/serial-gcc-4.8.2/include
SUNDIALS_INCLUDE_LIB_GNU   = -L/software/chimica2/libraries/sundials/sundials-2.6.2-opensmoke++/serial-gcc-4.8.2/lib
SUNDIALS_INCLUDE_INTEL     = -I/software/chimica2/libraries/sundials/sundials-2.6.2-opensmoke++/serial-intel-gcc-4.8.2/include
SUNDIALS_INCLUDE_LIB_INTEL = -L/software/chimica2/libraries/sundials/sundials-2.6.2-opensmoke++/serial-intel-gcc-4.8.2/lib
SUNDIALS_LIBS_INTEL        = -lsundials_ida -lsundials_kinsol -lsundials_cvode -lsundials_nvecserial -lklu -lbtf -lamd -lsuitesparseconfig -lcolamd 
SUNDIALS_LIBS_GNU          = -lsundials_ida -lsundials_kinsol -lsundials_cvode -lsundials_nvecserial -lklu -lbtf -lamd -lsuitesparseconfig -lcolamd

# -----------------------------------------------------------------------------------------------------------
# DASPK
# -----------------------------------------------------------------------------------------------------------
DASPK_INCLUDE_LIB = -L/software/chimica2/libraries/daspk/daspk-gcc-4.8.2/lib
DASPK_LIBS_INTEL  = -lDASPK20_LAPACK_intel_64bit
DASPK_LIBS_GNU    = -lDASPK20_LAPACK_gcc_64bit

# -----------------------------------------------------------------------------------------------------------
# ODEPACK
# -----------------------------------------------------------------------------------------------------------
ODEPACK_INCLUDE_LIB = -L/software/chimica2/libraries/odepack/odepack-gcc-4.8.2/lib
ODEPACK_LIBS_INTEL  = -lODEPACK_LAPACK_intel_64bit
ODEPACK_LIBS_GNU    = -lODEPACK_LAPACK_gcc_64bit

# -----------------------------------------------------------------------------------------------------------
# RADAU
# -----------------------------------------------------------------------------------------------------------
RADAU_INCLUDE_LIB = -L/software/chimica2/libraries/radau/radau-gcc-4.8.2/lib
RADAU_LIBS_INTEL  = -lRADAU_LAPACK_intel_64bit
RADAU_LIBS_GNU    = -lRADAU_LAPACK_gcc_64bit

# -----------------------------------------------------------------------------------------------------------
# MEBDF
# -----------------------------------------------------------------------------------------------------------
MEBDF_INCLUDE_LIB = -L/software/chimica2/libraries/mebdf/medf-gcc-4.8.2/lib
MEBDF_LIBS_INTEL  = -lMEBDF_LAPACK_intel_64bit
MEBDF_LIBS_GNU    = -lMEBDF_LAPACK_gcc_64bit

# -----------------------------------------------------------------------------------------------------------
# UMFPACK
# -----------------------------------------------------------------------------------------------------------
UMFPACK_INCLUDE     = -I/software/chimica2/libraries/suitesparse/suitesparse-4.5.2-intel-4.8.2/include
UMFPACK_INCLUDE_LIB = -L/software/chimica2/libraries/suitesparse/suitesparse-4.5.2-intel-4.8.2/lib
UMFPACK_LIBS_INTEL  = -lklu -lumfpack -lamd -lcholmod -lccolamd -lcamd -lsuitesparseconfig
UMFPACK_LIBS_GNU    = -lklu -lumfpack -lamd -lcholmod -lccolamd -lcamd -lsuitesparseconfig

# -----------------------------------------------------------------------------------------------------------
# METIS
# -----------------------------------------------------------------------------------------------------------
METIS_INCLUDE     = -I/software/chimica2/libraries/suitesparse/suitesparse-4.5.2-intel-4.8.2/include
METIS_INCLUDE_LIB = -L/software/chimica2/libraries/suitesparse/suitesparse-4.5.2-intel-4.8.2/lib
METIS_LIBS_INTEL  = 
METIS_LIBS_GNU    = 

# -----------------------------------------------------------------------------------------------------------
# SUPERLU_SERIAL
# -----------------------------------------------------------------------------------------------------------
SUPERLU_SERIAL_INCLUDE     = -I/software/chimica2/libraries/superlu/superlu/superlu-5.2-gcc-4.8.2/include
SUPERLU_SERIAL_INCLUDE_LIB = -L/software/chimica2/libraries/superlu/superlu/superlu-5.2-gcc-4.8.2/lib
SUPERLU_SERIAL_LIBS_INTEL  = -lSuperLU52_intel_64bit
SUPERLU_SERIAL_LIBS_GNU    = -lSuperLU52_gcc_64bit

# -----------------------------------------------------------------------------------------------------------
# LIS (TODO)
# -----------------------------------------------------------------------------------------------------------
LIS_INCLUDE     = 
LIS_INCLUDE_LIB = 
LIS_LIBS_INTEL  = 
LIS_LIBS_GNU    = 

# -----------------------------------------------------------------------------------------------------------
# BzzMath Libraries
# -----------------------------------------------------------------------------------------------------------
BZZMATH_INCLUDE        = -I/home/alberto/Development/ExternalNumericalLibraries/BzzMath/BzzMath6-dev/hpp
BZZMATH_INCLUDE_LIB    = -L/home/alberto/Development/ExternalNumericalLibraries/BzzMath/BzzMath6-dev/lib/linux
BZZMATH_LIBS_INTEL     = -lBzzMath6_INTEL
BZZMATH_LIBS_GNU       = -lBzzMath6_GNU

