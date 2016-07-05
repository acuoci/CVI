############################################################################
#
#  Program:         OpenSMOKE++
#
#  Module:          make.inc
#
#  Purpose:         Top-level Definitions
#
#  Creation date:   July 2014
#
############################################################################
# ------------------------------------------------------------------------------------------------
# PLEASE DO NOT CHANGE THIS FILE
# ------------------------------------------------------------------------------------------------
include make.compulsory.inc
include make.optional.inc

# ------------------------------------------------------------------------------------------------
# COMPILER OPTIONS
# ------------------------------------------------------------------------------------------------
ifeq ($(STATIC_LINKING),on)
OPTIONS_INTEL    =  -O3 -m64 -std=c++11
OPTIONS_GNU      =  -O3 -m64 -std=c++0x
else
OPTIONS_INTEL    =  -O3 -m64 -fPIC -std=c++11
OPTIONS_GNU      =  -O3 -m64 -fPIC -std=c++0x
endif

ifeq ($(OPENMP),on)
OPTIONS_INTEL    =  -O3 -m64 -std=c++0x -qopenmp  -I$(IFORT_PATH)/include 
OPTIONS_GNU      =  -O3 -m64 -std=c++0x -fopenmp  
endif

ifeq ($(MULTI_THREADED),on)
OPTIONS_INTEL    =  -O3 -m64 -std=c++0x -qopenmp
OPTIONS_GNU      =  -O3 -m64 -std=c++0x -fopenmp
endif

ifeq ($(PLASMA_SUPPORT),on)
OPTIONS_INTEL    =  -O2 -m64 -std=c++0x -diag-disable vec -DADD_
OPTIONS_GNU      =  -O2 -m64 -std=c++0x
endif


# ------------------------------------------------------------------------------------------------
# LIBRARY INTEL MKL
# ------------------------------------------------------------------------------------------------
ifeq ($(OPENSMOKE_USE_MKL),on)
	
	#Pre-processor directive	
	OPENSMOKE_USE_MKL_PREPROCESSOR = -DOPENSMOKE_USE_MKL=1
	MKL_PATH_LIB64=$(MKL_PATH)/lib/intel64

	#Groups of libraries
	ifeq ($(STATIC_LINKING),on)
		
		MKL_LIBS_GROUP        = -Wl,--start-group  $(MKL_PATH_LIB64)/libmkl_intel_lp64.a $(MKL_PATH_LIB64)/libmkl_sequential.a   $(MKL_PATH_LIB64)/libmkl_core.a -Wl,--end-group
		MULTI_THREADED_LIBS_GROUP     = -Wl,--start-group  $(MKL_PATH_LIB64)/libmkl_intel_lp64.a $(MKL_PATH_LIB64)/libmkl_intel_thread.a $(MKL_PATH_LIB64)/libmkl_core.a -Wl,--end-group
	
		ifeq ($(COMPILER_NAME),intel)
			IFORT_PATH_LIB64=$(IFORT_PATH)/lib/intel64
			IFORT_LIBS_GROUP      = -Wl,--start-group  $(IFORT_PATH_LIB64)/libifcore.a $(IFORT_PATH_LIB64)/libirc.a $(IFORT_PATH_LIB64)/libimf.a -Wl,--end-group
		else
			IFORT_PATH_LIB64	= $(IFORT_PATH)/lib/intel64
			IFORT_LIBS_GROUP    = -Wl,--start-group  $(IFORT_PATH_LIB64)/libifcore.a $(IFORT_PATH_LIB64)/libirc.a $(IFORT_PATH_LIB64)/libimf.a -Wl,--end-group
		endif
	else

		MKL_LIBS_GROUP        		= -Wl,--no-as-needed -L$(MKL_PATH_LIB64) -lmkl_intel_lp64 -lmkl_core -lmkl_sequential 
		MULTI_THREADED_LIBS_GROUP     	= -Wl,--no-as-needed -L$(MKL_PATH_LIB64) -lmkl_intel_lp64 -lmkl_core -lmkl_gnu_thread
	
		ifeq ($(COMPILER_NAME),intel)
			IFORT_PATH_LIB64    = $(IFORT_PATH)/lib/intel64
			IFORT_LIBS_GROUP    = -Wl,--no-as-needed  -L$(IFORT_PATH_LIB64) -lifcore -lirc -limf -lsvml
		else
			IFORT_PATH_LIB64    = $(IFORT_PATH)/lib/intel64
			IFORT_LIBS_GROUP    = -Wl,--no-as-needed  -L$(IFORT_PATH_LIB64) -lifcore -lirc -limf -lsvml
		endif
	
	endif

	#MKL libraries
	MKL_LIBS_INTEL  = $(MKL_LIBS_GROUP) -lpthread -lm 
	MKL_LIBS_GNU    = $(MKL_LIBS_GROUP) -lpthread -lm 

	ifeq ($(MULTI_THREADED),on)
		MKL_LIBS_INTEL  = $(MULTI_THREADED_LIBS_GROUP) -lpthread -lm 
		MKL_LIBS_GNU    = $(MULTI_THREADED_LIBS_GROUP) -ldl -lpthread -lm 
	endif

	#IFORT libraries
	IFORT_LIBS_INTEL  = $(IFORT_LIBS_GROUP) -ldl
	IFORT_LIBS_GNU    = $(IFORT_LIBS_GROUP) -ldl

else

	#Pre-processor directive	
	OPENSMOKE_USE_MKL_PREPROCESSOR=

	#MKL libraries
	MKL_PATH        =
	MKL_LIBS_INTEL  =
	MKL_LIBS_GNU    =

	#IFORT libraries
	IFORT_PATH        =
	IFORT_LIBS_INTEL  =
	IFORT_LIBS_GNU    =
endif

# ------------------------------------------------------------------------------------------------
# LIBRARY OPENBLAS
# ------------------------------------------------------------------------------------------------
ifeq ($(OPENSMOKE_USE_OPENBLAS),on)
	OPENSMOKE_USE_OPENBLAS_PREPROCESSOR = -DOPENSMOKE_USE_OPENBLAS=1
else
	OPENSMOKE_USE_OPENBLAS_PREPROCESSOR =
	OPENBLAS_INCLUDE              =
        OPENBLAS_INCLUDE_LIB          =
	OPENBLAS_LIBS_INTEL           =
	OPENBLAS_LIBS_GNU             =
endif


# ------------------------------------------------------------------------------------------------
# LIBRARY DASPK
# ------------------------------------------------------------------------------------------------
ifeq ($(OPENSMOKE_USE_DASPK),on)
	OPENSMOKE_USE_DASPK_PREPROCESSOR = -DOPENSMOKE_USE_DASPK=1
else
	OPENSMOKE_USE_DASPK_PREPROCESSOR =
	DASPK_INCLUDE_LIB          =
	DASPK_LIBS_INTEL           =
	DASPK_LIBS_GNU             =
endif

# ------------------------------------------------------------------------------------------------
# LIBRARY MEBDF
# ------------------------------------------------------------------------------------------------
ifeq ($(OPENSMOKE_USE_MEBDF),on)
	OPENSMOKE_USE_MEBDF_PREPROCESSOR = -DOPENSMOKE_USE_MEBDF=1
else
	OPENSMOKE_USE_MEBDF_PREPROCESSOR =
        MEBDF_INCLUDE_LIB          =
	MEBDF_LIBS_INTEL     =
	MEBDF_LIBS_GNU       =
endif

# ------------------------------------------------------------------------------------------------
# LIBRARY ODEPACK
# ------------------------------------------------------------------------------------------------
ifeq ($(OPENSMOKE_USE_ODEPACK),on)
	OPENSMOKE_USE_ODEPACK_PREPROCESSOR = -DOPENSMOKE_USE_ODEPACK=1
else
	OPENSMOKE_USE_ODEPACK_PREPROCESSOR =
        ODEPACK_INCLUDE_LIB          =
	ODEPACK_LIBS_INTEL           =
	ODEPACK_LIBS_GNU             =
endif

# ------------------------------------------------------------------------------------------------
# LIBRARY RADAU
# ------------------------------------------------------------------------------------------------
ifeq ($(OPENSMOKE_USE_RADAU),on)
	OPENSMOKE_USE_RADAU_PREPROCESSOR = -DOPENSMOKE_USE_RADAU=1
else
	OPENSMOKE_USE_RADAU_PREPROCESSOR =
        RADAU_INCLUDE_LIB          =
	RADAU_LIBS_INTEL           =
	RADAU_LIBS_GNU             =
endif

# ------------------------------------------------------------------------------------------------
# LIBRARY SUNDIALS
# ------------------------------------------------------------------------------------------------
ifeq ($(OPENSMOKE_USE_SUNDIALS),on)
	OPENSMOKE_USE_SUNDIALS_PREPROCESSOR = -DOPENSMOKE_USE_SUNDIALS=1
	ifeq ($(COMPILER_NAME),intel)
		SUNDIALS_INCLUDE     = $(SUNDIALS_INCLUDE_INTEL)
        	SUNDIALS_INCLUDE_LIB = $(SUNDIALS_INCLUDE_LIB_INTEL)
	else
		SUNDIALS_INCLUDE     = $(SUNDIALS_INCLUDE_GNU)
        	SUNDIALS_INCLUDE_LIB = $(SUNDIALS_INCLUDE_LIB_GNU)
	endif
else
	OPENSMOKE_USE_SUNDIALS_PREPROCESSOR =
        SUNDIALS_INCLUDE     =
        SUNDIALS_INCLUDE_LIB =
	SUNDIALS_LIBS_INTEL  =
	SUNDIALS_LIBS_GNU    =
endif

# ------------------------------------------------------------------------------------------------
# LIBRARY UMFPACK
# ------------------------------------------------------------------------------------------------
ifeq ($(OPENSMOKE_USE_UMFPACK),on)
	OPENSMOKE_USE_UMFPACK_PREPROCESSOR = -DOPENSMOKE_USE_UMFPACK=1
else
	OPENSMOKE_USE_UMFPACK_PREPROCESSOR =
	UMFPACK_INCLUDE_LIB          =
	UMFPACK_LIBS_INTEL           =
	UMFPACK_LIBS_GNU             =
endif

# ------------------------------------------------------------------------------------------------
# LIBRARY SUPERLUSERIAL
# ------------------------------------------------------------------------------------------------
ifeq ($(OPENSMOKE_USE_SUPERLU_SERIAL),on)
	OPENSMOKE_USE_SUPERLU_SERIAL_PREPROCESSOR = -DOPENSMOKE_USE_SUPERLU_SERIAL=1
else
	OPENSMOKE_USE_SUPERLU_SERIAL_PREPROCESSOR =
	SUPERLU_SERIAL_INCLUDE_LIB          =
	SUPERLU_SERIAL_LIBS_INTEL           =
	SUPERLU_SERIAL_LIBS_GNU             =
endif

# ------------------------------------------------------------------------------------------------
# LIBRARY LIS
# ------------------------------------------------------------------------------------------------
ifeq ($(OPENSMOKE_LIS),on)
	OPENSMOKE_LIS_PREPROCESSOR = -DOPENSMOKE_LIS=1
else
	OPENSMOKE_USE_LIS_PREPROCESSOR =
	LIS_INCLUDE_LIB          =
	LIS_LIBS_INTEL           =
	LIS_LIBS_GNU             =
endif

# ------------------------------------------------------------------------------------------------
# LIBRARY BZZMATH
# ------------------------------------------------------------------------------------------------
ifeq ($(OPENSMOKE_USE_BZZMATH),on)
	OPENSMOKE_USE_BZZMATH_PREPROCESSOR1 = -DOPENSMOKE_USE_BZZMATH=1
	OPENSMOKE_USE_BZZMATH_PREPROCESSOR2 = -DBZZMATH_USE_MKL=0
else
	OPENSMOKE_USE_BZZMATH_PREPROCESSOR1 =
	OPENSMOKE_USE_BZZMATH_PREPROCESSOR2 =
        BZZMATH_INCLUDE               =
        BZZMATH_INCLUDE_LIB           =
	BZZMATH_LIBS_INTEL            =
	BZZMATH_LIBS_GNU              =
endif

# BzzMath with MKL features can be used only together with the Intel compiler (bug to be fixed)
ifeq ($(COMPILER_NAME),intel)
	ifeq ($(OPENSMOKE_USE_BZZMATH),on)
		ifeq ($(OPENSMOKE_USE_MKL),on)
			OPENSMOKE_USE_BZZMATH_PREPROCESSOR2 = -DBZZMATH_USE_MKL=0
		endif
		ifeq ($(OPENSMOKE_USE_OPENBLAS),on)
			OPENSMOKE_USE_BZZMATH_PREPROCESSOR2 = -DBZZMATH_USE_MKL=0
		endif
	endif
endif

OPENSMOKE_USE_BZZMATH_PREPROCESSOR = $(OPENSMOKE_USE_BZZMATH_PREPROCESSOR1) $(OPENSMOKE_USE_BZZMATH_PREPROCESSOR2)

# ------------------------------------------------------------------------------------------------
# LIBRARY LINPACK
# ------------------------------------------------------------------------------------------------
ifeq ($(OPENSMOKE_USE_DASPK),on)
	LINPACK_SUPPORT=on
endif
ifeq ($(OPENSMOKE_USE_MEBDF),on)
	LINPACK_SUPPORT=on
endif
ifeq ($(OPENSMOKE_USE_ODEPACK),on)
	LINPACK_SUPPORT=on
endif
ifeq ($(OPENSMOKE_USE_RADAU),on)
	LINPACK_SUPPORT=on
endif

ifneq ($(LINPACK_SUPPORT),on) 
        LINPACK_INCLUDE_LIB  =
	LINPACK_LIBS_INTEL   =
	LINPACK_LIBS_GNU     =
endif



# ------------------------------------------------------------------------------------------------
# LIBRARY PLASMA
# ------------------------------------------------------------------------------------------------
ifeq ($(PLASMA_SUPPORT),on)
PLASMA_LIBS_INTEL  = -Wl,--start-group  -lplasma -lcoreblas -lquark -llapacke -Wl,--end-group -ldl
PLASMA_LIBS_GNU    = -Wl,--start-group  -lplasma -lcoreblas -lquark -llapacke -Wl,--end-group -ldl
PLASMA_INCLUDE     =
PLASMA_INCLUDE_LIB =
endif
ifeq ($(PLASMA_SUPPORT),off)
PLASMA_INCLUDE     =
PLASMA_INCLUDE_LIB =
PLASMA_LIBS_INTEL  =
PLASMA_LIBS_GNU    =
endif

# ------------------------------------------------------------------------------------------------
# CPREPROCESSOR DIRECTIVES
# ------------------------------------------------------------------------------------------------
PREPROCESSOR_DIRECTIVES = $(OPENSMOKE_USE_MKL_PREPROCESSOR) $(OPENSMOKE_USE_OPENBLAS_PREPROCESSOR) $(OPENSMOKE_USE_SUNDIALS_PREPROCESSOR) $(OPENSMOKE_USE_ODEPACK_PREPROCESSOR) $(OPENSMOKE_USE_DASPK_PREPROCESSOR) $(OPENSMOKE_USE_RADAU_PREPROCESSOR) $(OPENSMOKE_USE_MEBDF_PREPROCESSOR) $(OPENSMOKE_USE_BZZMATH_PREPROCESSOR) $(OPENSMOKE_USE_UMFPACK_PREPROCESSOR) $(OPENSMOKE_USE_SUPERLU_SERIAL_PREPROCESSOR) $(OPENSMOKE_USE_LIS_PREPROCESSOR)

# ------------------------------------------------------------------------------------------------
# COMPILERS
# ------------------------------------------------------------------------------------------------
ifeq ($(COMPILER_NAME),intel)

	COMPILER_TAG   = intel
	STRING_NAME    =
	ifeq ($(MULTI_THREADED),on)
	STRING_NAME    =
	endif
	CXX            = icpc
	CXX_FLAGS      = $(OPTIONS_INTEL) -I/usr/include/x86_64-linux-gnu/c++/4.8/ $(PREPROCESSOR_DIRECTIVES)

	BZZMATH_LIBS   = $(BZZMATH_LIBS_INTEL)
	DASPK_LIBS     = $(DASPK_LIBS_INTEL)
	MEBDF_LIBS     = $(MEBDF_LIBS_INTEL)
	ODEPACK_LIBS   = $(ODEPACK_LIBS_INTEL)
	LINPACK_LIBS   = $(LINPACK_LIBS_INTEL)
	RADAU_LIBS     = $(RADAU_LIBS_INTEL)
	SUNDIALS_LIBS  = $(SUNDIALS_LIBS_INTEL)
	UMFPACK_LIBS   = $(UMFPACK_LIBS_INTEL)
	METIS_LIBS     = $(METIS_LIBS_INTEL)
	SUPERLU_SERIAL_LIBS   = $(SUPERLU_SERIAL_LIBS_INTEL)
	LIS_LIBS       = $(LIS_LIBS_INTEL)
	MKL_LIBS       = $(MKL_LIBS_INTEL)
	OPENBLAS_LIBS  = $(OPENBLAS_LIBS_INTEL)
	IFORT_LIBS     = $(IFORT_LIBS_INTEL)
	PLASMA_LIBS    = $(PLASMA_LIBS_INTEL)

endif	


ifeq ($(COMPILER_NAME),gnu)

	COMPILER_TAG   = gnu
	STRING_NAME    =
	ifeq ($(MULTI_THREADED),on)
	STRING_NAME    =
	endif
	CXX            = g++
	CXX_FLAGS      = $(OPTIONS_GNU) $(PREPROCESSOR_DIRECTIVES)

	BZZMATH_LIBS   = $(BZZMATH_LIBS_GNU)
	DASPK_LIBS     = $(DASPK_LIBS_GNU)
	MEBDF_LIBS     = $(MEBDF_LIBS_GNU)
	ODEPACK_LIBS   = $(ODEPACK_LIBS_GNU)
	LINPACK_LIBS   = $(LINPACK_LIBS_GNU)
	RADAU_LIBS     = $(RADAU_LIBS_GNU)
	SUNDIALS_LIBS  = $(SUNDIALS_LIBS_GNU)
	UMFPACK_LIBS   = $(UMFPACK_LIBS_GNU) 
        METIS_LIBS     = $(METIS_LIBS_GNU)
	SUPERLU_SERIAL_LIBS   = $(SUPERLU_SERIAL_LIBS_GNU)
	LIS_LIBS       = $(LIS_LIBS_GNU)
	MKL_LIBS       = $(MKL_LIBS_GNU)
	OPENBLAS_LIBS  = $(OPENBLAS_LIBS_GNU)
	IFORT_LIBS     = $(IFORT_LIBS_GNU)
	PLASMA_LIBS    = $(PLASMA_LIBS_GNU)

endif	

# ------------------------------------------------------------------------------------------------
# INCLUDE
# ------------------------------------------------------------------------------------------------
OPENSMOKESUITE_INCLUDE = -I../../src
OPENSMOKE_INCLUDE      = -I$(OPENSMOKE_PATH)
MKL_INCLUDE            = -I$(MKL_PATH)/include
RAPIDXML_INCLUDE       = -I$(RAPIDXML_PATH)
EIGEN_INCLUDE          = -I$(EIGEN_PATH)

# ------------------------------------------------------------------------------------------------
# GLOBAL
# ------------------------------------------------------------------------------------------------
GLOBAL_INCLUDE     = $(OPENSMOKE_INCLUDE) $(BZZMATH_INCLUDE) $(SUNDIALS_INCLUDE) $(UMFPACK_INCLUDE) $(METIS_INCLUDE) $(SUPERLU_SERIAL_INCLUDE) $(LIS_INCLUDE) $(BOOST_INCLUDE) $(MKL_INCLUDE) $(OPENBLAS_INCLUDE) $(RAPIDXML_INCLUDE) $(EIGEN_INCLUDE) $(PLASMA_INCLUDE) $(OPENSMOKESUITE_INCLUDE)

GLOBAL_LIB_INCLUDE = $(OPENBLAS_INCLUDE_LIB) $(BZZMATH_INCLUDE_LIB) $(BOOST_INCLUDE_LIB) $(DASPK_INCLUDE_LIB) $(MEBDF_INCLUDE_LIB) $(ODEPACK_INCLUDE_LIB) $(LINPACK_INCLUDE_LIB) $(RADAU_INCLUDE_LIB) $(SUNDIALS_INCLUDE_LIB) $(UMFPACK_INCLUDE_LIB) $(METIS_INCLUDE_LIB) $(SUPERLU_SERIAL_INCLUDE_LIB) $(LIS_INCLUDE_LIB) $(PLASMA_INCLUDE_LIB)

GLOBAL_LIBS        = $(BOOST_LIBS) $(BZZMATH_LIBS) -Wl,--start-group $(DASPK_LIBS) $(MEBDF_LIBS) $(ODEPACK_LIBS) $(LINPACK_LIBS) -Wl,--end-group $(RADAU_LIBS) $(SUNDIALS_LIBS) $(UMFPACK_LIBS) $(METIS_LIBS) $(SUPERLU_SERIAL_LIBS) $(LIS_LIBS) $(PLASMA_LIBS) $(MKL_LIBS) $(OPENBLAS_LIBS) $(IFORT_LIBS) -lgfortran

GLOBAL_LIBS_MINIMALIST = $(BOOST_LIBS) $(MKL_LIBS) $(OPENBLAS_LIBS) $(IFORT_LIBS) -lgfortran

GLOBAL_LIBS_BZZMATH_AND_SUNDIALS = $(BOOST_LIBS) $(BZZMATH_LIBS) $(SUNDIALS_LIBS) $(MKL_LIBS) $(OPENBLAS_LIBS) $(IFORT_LIBS) -lgfortran

# ------------------------------------------------------------------------------------------------
# COMPILER INSTRUCTIONS
# ------------------------------------------------------------------------------------------------
CCP_COMPILE = $(CXX) $(CXX_FLAGS) -c $(GLOBAL_INCLUDE)  

