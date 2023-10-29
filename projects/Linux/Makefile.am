AUTOMAKE_OPTIONS = foreign

if INTEL_COMPILER
AM_CXXFLAGS = -O3 -m64 -fPIC -std=c++17 @OPENMP_INTEL@
AM_FCFLAGS  = -O3 -m64
AM_F77FLAGS = -O3 -m64
else
AM_CXXFLAGS = -O3 -m64 -fPIC -std=c++17 @OPENMP_GNU@
AM_FCFLAGS  = -O3 -m64
AM_F77FLAGS = -O3 -m64
endif

AM_CPPFLAGS = @USE_MKL@ @USE_EXTERNAL_LIBRARIES@ \
	      @COMPULSORY_INCLUDE_DIRS@ @EXTERNAL_INCLUDE_DIRS@ @MKL_INCLUDE_DIR@

AM_LDFLAGS  = @BOOST_LIBS@ @MKL_LIBS@ -lgfortran

#List of libraries
noinst_LIBRARIES = libblocktdsge.a
libblocktdsge_a_SOURCES = /home/chimica2/cuoci/Development/MyCodes/OpenSMOKEppSolvers/src/numerics/blocktdsge/dgeblttrf.f \
			  /home/chimica2/cuoci/Development/MyCodes/OpenSMOKEppSolvers/src/numerics/blocktdsge/dgeblttrs.f \
			  /home/chimica2/cuoci/Development/MyCodes/OpenSMOKEppSolvers/src/numerics/blocktdsge/auxfunctions.f

#Main program
bin_PROGRAMS = OpenSMOKEpp_CVI.sh
OpenSMOKEpp_CVI_sh_SOURCES = ../../src/solver.cpp
OpenSMOKEpp_CVI_sh_LDADD = libblocktdsge.a @EXTERNAL_LIBS@



