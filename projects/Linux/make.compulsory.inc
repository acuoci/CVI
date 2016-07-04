# --------------------------------------------------------------------------------------------------------------------
# COMPULSORY INFORMATION
# 1. COMPILER_NAME: gnu | intel
# --------------------------------------------------------------------------------------------------------------------
# Choose the compiler: gnu | intel
COMPILER_NAME  = intel
MULTI_THREADED = off
OPENMP         = off
STATIC_LINKING = off
PATH_EXE       = ../../exe/linux/$(COMPILER_NAME)64-serial-mkl

# --------------------------------------------------------------------------------------------------------------------
# COMPULSORY INCLUDE PATHS (please change according to the paths on your machine)
# --------------------------------------------------------------------------------------------------------------------
OPENSMOKE_PATH         =  $(HOME)/Development/MyCodes/OpenSMOKEpp/source
EIGEN_PATH 	       =  /software/chimica2/libraries/eigen/eigen-3.3-dev
RAPIDXML_PATH 	       =  /software/chimica2/libraries/rapidxml/rapidxml-1.13
# --------------------------------------------------------------------------------------------------------------------

# --------------------------------------------------------------------------------------------------------------------
# COMPULSORY BOOST LIBRARIES (please change according to the names of the libraries on your machine)
# --------------------------------------------------------------------------------------------------------------------
BOOST_INCLUDE       = -I/software/chimica2/libraries/boost/boost-1.60.0-gcc-4.8.2/include
BOOST_INCLUDE_LIB   = -L/software/chimica2/libraries/boost/boost-1.60.0-gcc-4.8.2/lib
BOOST_LIBS          = -lboost_date_time -lboost_filesystem -lboost_program_options -lboost_system \
                      -lboost_regex -lboost_timer -lboost_chrono
# --------------------------------------------------------------------------------------------------------------------

# --------------------------------------------------------------------------------------------------------------------
# OPTIONAL EXTERNAL LIBRARIES (if on the correct paths have to be provided in the Make.optional.inc file)
# 1. only one between OPENSMOKE_USE_MKL and OPENSMOKE_USE_OPENBLAS can be turned on
# 2. the MKL or the OPENBLAS support is strictly required by: SUNDIALS, ODEPACK, DASPK, RADAU, MEBDF, BZZMATH
# 2. the LINPACK Libraries must be provided (in the Make.optional.inc file) for: DASPK, MEBDF, ODEPACK, RADAU
# --------------------------------------------------------------------------------------------------------------------
OPENSMOKE_USE_MKL             =    on
OPENSMOKE_USE_OPENBLAS        =    off
OPENSMOKE_USE_SUNDIALS        =    off
OPENSMOKE_USE_UMFPACK         =    on
OPENSMOKE_USE_SUPERLU_SERIAL  =    on
OPENSMOKE_USE_LIS             =    off
OPENSMOKE_USE_ODEPACK         =    off
OPENSMOKE_USE_DASPK           =    off
OPENSMOKE_USE_RADAU           =    off
OPENSMOKE_USE_MEBDF           =    off
OPENSMOKE_USE_BZZMATH         =    on

