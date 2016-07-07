
AC_INIT([OpenSMOKEpp_CVI], [0.1], [alberto.cuoci@polimi.it])

AM_INIT_AUTOMAKE

AC_PROG_CC
AC_PROG_CXX

AC_PROG_CC([gcc icc])
AC_PROG_CXX([g++ icpc])

CXXFLAGS=''

AM_CONDITIONAL([INTEL_COMPILER],[test "$CXX" = "icpc"])

# ---------------------------------------------------------------------------------------------------------------------------------------
# Recognize if OpenMP libraries are required by the user
# ---------------------------------------------------------------------------------------------------------------------------------------
OPENMP_INTEL=""
OPENMP_GNU=""

AC_ARG_ENABLE([openmp], AS_HELP_STRING([--enable-openmp], [Enable OpenMP support]))
AS_IF([test "x$enable_openmp" = "xyes"], [OPENMP_INTEL=" -qopenmp"; OPENMP_GNU=" -fopenmp "])

AC_SUBST(OPENMP_INTEL)
AC_SUBST(OPENMP_GNU)

# ---------------------------------------------------------------------------------------------------------------------------------------
# Add OpenSMOKE++ library (compulsory)
# ---------------------------------------------------------------------------------------------------------------------------------------
AC_ARG_WITH([opensmoke], [AS_HELP_STRING([--with-opensmoke=DIR], 
            [installation prefix for OpenSMOKE++ libraries])], 
            [with_opensmoke_include_path=$with_opensmoke/source],
	    AC_MSG_ERROR([Please use --with-opensmokepp to provide the path to OpenSMOKE++. Ask alberto.cuoci@polimi.it for the most recent version of OpenSMOKE++]))

AC_CHECK_FILE([${with_opensmoke_include_path}/OpenSMOKEpp], [], 
              AC_MSG_ERROR([Include file OpenSMOKEpp not found. Ask alberto.cuoci@polimi.it for the most recent version of OpenSMOKE++]))

OPENSMOKEPP_INCLUDE_DIR="-I$with_opensmoke_include_path -I../../src"

# ---------------------------------------------------------------------------------------------------------------------------------------
# Add Eigen library (compulsory)
# ---------------------------------------------------------------------------------------------------------------------------------------
AC_ARG_WITH([eigen], [AS_HELP_STRING([--with-eigen=DIR], 
            [installation prefix for Eigen libraries])], 
            [with_eigen_include_path=$with_eigen],
	    AC_MSG_ERROR([Please use --with-eigen to provide the path to Eigen. Try downloading the Eigen libraries from http://eigen.tuxfamily.org/index.php?title=Main_Page]))

AC_CHECK_FILE([${with_eigen_include_path}/Eigen/Eigen], [], 
              AC_MSG_ERROR([Include file Eigen not found. Try downloading the Eigen libraries from http://eigen.tuxfamily.org/index.php?title=Main_Page]))

EIGEN_INCLUDE_DIR="-I$with_eigen_include_path"

# ---------------------------------------------------------------------------------------------------------------------------------------
# Add RapidXML library (compulsory)
# ---------------------------------------------------------------------------------------------------------------------------------------
AC_ARG_WITH([rapidxml], [AS_HELP_STRING([--with-rapidxml=DIR], 
            [installation prefix for RapidXML libraries])], 
            [with_rapidxml_include_path=$with_rapidxml],
	    AC_MSG_ERROR([Please use --with-rapidxml to provide the path to RapidXML. Try downloading the RapidXML libraries from http://rapidxml.sourceforge.net]))

AC_CHECK_FILE([${with_rapidxml_include_path}/rapidxml.hpp], [], 
              AC_MSG_ERROR([Include file rapidxml.hpp not found. Try downloading the RapidXML libraries from http://rapidxml.sourceforge.net]))

RAPIDXML_INCLUDE_DIR="-I$with_rapidxml_include_path"

# ---------------------------------------------------------------------------------------------------------------------------------------
# Add Boost C++ library (compulsory, but system implementation can be used)
# ---------------------------------------------------------------------------------------------------------------------------------------

SYSTEM_BOOST="yes"

AC_ARG_WITH([boost-include], [AS_HELP_STRING([--with-boost-include=DIR], 
            [installation prefix (include) for Boost C++ libraries])], 
            [with_boost_include_path=$with_boost_include; SYSTEM_BOOST="no"])

AC_ARG_WITH([boost-lib], [AS_HELP_STRING([--with-boost-lib=DIR], 
            [installation prefix (lib) for Boost C++ libraries])], 
            [with_boost_lib_path=$with_boost_lib; SYSTEM_BOOST="no"])


AS_IF([test "$SYSTEM_BOOST" = "yes"], [ 

	AC_CHECK_HEADERS([regex.h], [], [Include file regex.h not found. Try downloading the Boost C++ libraries from http://www.boost.org], [])

	BOOST_INCLUDE_DIR=""
	BOOST_LIBS="-lboost_date_time -lboost_filesystem -lboost_program_options -lboost_system -lboost_regex -lboost_timer -lboost_chrono";
])


AS_IF([test "$SYSTEM_BOOST" = "no"], [

	AC_CHECK_FILE([${with_boost_include_path}/boost/regex.h], [], 
	              AC_MSG_ERROR([Include file regex.h not found. Try downloading the Boost C++ libraries from http://www.boost.org]))

	AC_CHECK_FILE([${with_boost_lib_path}/libboost_date_time.a], [], 
		      AC_CHECK_FILE([${with_boost_lib_path}/libboost_date_time.so], [], 
		      AC_MSG_ERROR([Library file libboost_date_time.a or libboost_date_time.so not found. Try downloading the Boost C++ libraries from http://www.boost.org])))

	AC_CHECK_FILE([${with_boost_lib_path}/libboost_filesystem.a], [], 
		      AC_CHECK_FILE([${with_boost_lib_path}/libboost_filesystem.so], [], 
		      AC_MSG_ERROR([Library file libboost_filesystem.a or libboost_filesystem.so not found. Try downloading the Boost C++ libraries from http://www.boost.org])))

	AC_CHECK_FILE([${with_boost_lib_path}/libboost_program_options.a], [], 
		      AC_CHECK_FILE([${with_boost_lib_path}/libboost_program_options.so], [], 
		      AC_MSG_ERROR([Library file libboost_program_options.a or libboost_program_options.so not found. Try downloading the Boost C++ libraries from http://www.boost.org])))

	AC_CHECK_FILE([${with_boost_lib_path}/libboost_system.a], [], 
		      AC_CHECK_FILE([${with_boost_lib_path}/libboost_system.so], [], 
		      AC_MSG_ERROR([Library file libboost_system.a or libboost_system.so not found. Try downloading the Boost C++ libraries from http://www.boost.org])))

	AC_CHECK_FILE([${with_boost_lib_path}/libboost_regex.a], [], 
		      AC_CHECK_FILE([${with_boost_lib_path}/libboost_regex.so], [], 
		      AC_MSG_ERROR([Library file libboost_regex.a or libboost_regex.so not found. Try downloading the Boost C++ libraries from http://www.boost.org])))

	AC_CHECK_FILE([${with_boost_lib_path}/libboost_timer.a], [], 
		      AC_CHECK_FILE([${with_boost_lib_path}/libboost_timer.so], [], 
		      AC_MSG_ERROR([Library file libboost_timer.a or libboost_timer.so not found. Try downloading the Boost C++ libraries from http://www.boost.org])))

	AC_CHECK_FILE([${with_boost_lib_path}/libboost_chrono.a], [], 
		      AC_CHECK_FILE([${with_boost_lib_path}/libboost_chrono.so], [], 
		      AC_MSG_ERROR([Library file libboost_chrono.a or libboost_chrono.so not found. Try downloading the Boost C++ libraries from http://www.boost.org])))

	BOOST_INCLUDE_DIR="-I$with_boost_include_path"
	BOOST_LIBS="-L$with_boost_lib_path -lboost_date_time -lboost_filesystem -lboost_program_options -lboost_system -lboost_regex -lboost_timer -lboost_chrono"
])

# ---------------------------------------------------------------------------------------------------------------------------------------
# Add Intel MKL library (optional)
# ---------------------------------------------------------------------------------------------------------------------------------------

USE_MKL=" -DOPENSMOKE_USE_MKL=0 "
MKL_INCLUDE_DIR=""
MKL_LIBS=""

AC_ARG_WITH([mkl], [AS_HELP_STRING([--with-mkl=DIR], 
            [installation prefix for Intel MKL libraries])], 
            [ with_mkl_include_path=$with_mkl/include; 
	      with_mkl_lib_path=$with_mkl/lib/intel64; 
	      USE_MKL=" -DOPENSMOKE_USE_MKL=1 "; 
	      MKL_INCLUDE_DIR=" -I$with_mkl_include_path "; 
	      MKL_LIBS=" -Wl,--no-as-needed -L$with_mkl_lib_path -lmkl_intel_lp64 -lmkl_core -lmkl_sequential ";
	      AC_CHECK_FILE( [${with_mkl_include_path}/mkl.h], [], 
                             AC_MSG_ERROR([Include file mkl.h not found. Try downloading the Intel MKL libraries from https://software.intel.com/en-us/intel-mkl])) ])


# ---------------------------------------------------------------------------------------------------------------------------------------
# Add OpenBLAS (optional, as alternative to MKL)
# ---------------------------------------------------------------------------------------------------------------------------------------

USE_OPENBLAS=" -DOPENSMOKE_USE_OPENBLAS=0 "
OPENBLAS_INCLUDE_DIR=""
OPENBLAS_LIBS=""

AC_ARG_WITH([openblas], [AS_HELP_STRING([--with-openblas=DIR], 
            [installation prefix for OpenBLAS libraries])], 
            [ with_openblas_include_path=$with_openblas/include; 
	      with_openblas_lib_path=$with_openblas/lib; 
	      USE_OPENBLAS=" -DOPENSMOKE_USE_OPENBLAS=1 "; 
	      OPENBLAS_INCLUDE_DIR=" -I$with_openblas_include_path "; 
	      OPENBLAS_LIBS=" -Wl,--no-as-needed -L$with_openblas_lib_path -lopenblas ";
	      AC_CHECK_FILE( [${with_openblas_include_path}/lapacke.h], [], 
                             AC_MSG_ERROR([Include file lapacke.h not found. Try downloading the OpenBLAS libraries from http://www.openblas.net/])) ])


# ---------------------------------------------------------------------------------------------------------------------------------------
# Add Intel compiler (optional)
# ---------------------------------------------------------------------------------------------------------------------------------------

INTEL_COMPILER_LIBS=""

AC_ARG_WITH([intel-compiler], [AS_HELP_STRING([--with-intel-compiler=DIR], 
            [installation prefix for Intel Compiler])], 
            [ with_intel_compiler_path=$with_intel_compiler/lib/intel64;
	      INTEL_COMPILER_LIBS=" -L$with_intel_compiler_path -Wl,--start-group -lifcore -lirc -limf -Wl,--end-group ";
	      AC_CHECK_FILE( [${with_intel_compiler_path}/libifcore.so], [], 
                             AC_MSG_ERROR([Library file ${with_intel_compiler_path}/libifcore.so not found. Try downloading the Intel Compiler from https://software.intel.com])) ])

# ---------------------------------------------------------------------------------------------------------------------------------------
# Add BzzMath library (optional)
# ---------------------------------------------------------------------------------------------------------------------------------------

USE_BZZMATH=" -DOPENSMOKE_USE_BZZMATH=0 "
BZZMATH_INCLUDE_DIR=""
BZZMATH_LIBS=""

AC_ARG_WITH([bzzmath-include], [AS_HELP_STRING([--with-bzzmath-include=DIR], 
            [installation prefix (include) for BzzMath libraries])], 
            [ with_bzzmath_include_path=$with_bzzmath_include; USE_BZZMATH=" -DOPENSMOKE_USE_BZZMATH=1 "; BZZMATH_INCLUDE_DIR=" -I$with_bzzmath_include_path ";

		AC_CHECK_FILE([${with_bzzmath_include_path}/BzzMath.hpp], [], 
              	AC_MSG_ERROR([Include file BzzMath.hpp not found. Try downloading the BzzMath libraries from http://super.chem.polimi.it/download/bzzmath-download])) ])

AC_ARG_WITH([bzzmath-lib], [AS_HELP_STRING([--with-bzzmath-lib=PATH], 
            [path to the BzzMath libraries])], 
            [ with_bzzmath_lib_path=$with_bzzmath_lib; USE_BZZMATH=" -DOPENSMOKE_USE_BZZMATH=1 "; BZZMATH_LIBS=" $with_bzzmath_lib_path ";

		AC_CHECK_FILE([${with_bzzmath_lib_path}], [], 
              	AC_MSG_ERROR([Library file ${with_bzzmath_lib_path} not found. Try downloading the BzzMath libraries from http://super.chem.polimi.it/download/bzzmath-download])) ])


# ---------------------------------------------------------------------------------------------------------------------------------------
# Add Sundials library (optional)
# ---------------------------------------------------------------------------------------------------------------------------------------

USE_SUNDIALS=" -DOPENSMOKE_USE_SUNDIALS=0 "
SUNDIALS_INCLUDE_DIR=""
SUNDIALS_LIBS=""

AC_ARG_WITH([sundials-include], [AS_HELP_STRING([--with-sundials-include=DIR], 
            [installation prefix (include) for Sundials libraries])], 
            [ with_sundials_include_path=$with_sundials_include; USE_SUNDIALS=" -DOPENSMOKE_USE_SUNDIALS=1 "; SUNDIALS_INCLUDE_DIR=" -I$with_sundials_include_path ";

		AC_CHECK_FILE([${with_sundials_include_path}/sundials/sundials_config.h], [], 
              	AC_MSG_ERROR([Include file sundials_config.h not found. Try downloading the Sundials libraries from http://computation.llnl.gov/projects/sundials-suite-nonlinear-differential-algebraic-equation-solvers])) ])

AC_ARG_WITH([sundials-lib], [AS_HELP_STRING([--with-sundials-lib=PATH], 
            [path to the Sundials libraries])], 
            [ with_sundials_lib_path=$with_sundials_lib; USE_SUNDIALS=" -DOPENSMOKE_USE_SUNDIALS=1 "; SUNDIALS_LIBS=" -L$with_sundials_lib_path -lsundials_ida -lsundials_kinsol -lsundials_cvode -lsundials_nvecserial -lklu -lbtf -lamd -lsuitesparseconfig -lcolamd ";

		AC_CHECK_FILE([${with_sundials_lib_path}/libsundials_nvecserial.so], [], 
              	AC_MSG_ERROR([Library file ${with_sundials_lib_path}/libsundials_nvecserial.so not found. Try downloading the Sundials libraries from http://computation.llnl.gov/projects/sundials-suite-nonlinear-differential-algebraic-equation-solvers])) ])


# ---------------------------------------------------------------------------------------------------------------------------------------
# Add UMFPACK library (optional)
# ---------------------------------------------------------------------------------------------------------------------------------------

USE_UMFPACK=" -DOPENSMOKE_USE_UMFPACK=0 "
UMFPACK_INCLUDE_DIR=""
UMFPACK_LIBS=""

AC_ARG_WITH([umfpack-include], [AS_HELP_STRING([--with-umfpack-include=DIR], 
            [installation prefix (include) for SuiteSparse (UMFPACK) libraries])], 
            [ with_umfpack_include_path=$with_umfpack_include; USE_UMFPACK=" -DOPENSMOKE_USE_UMFPACK=1 "; UMFPACK_INCLUDE_DIR=" -I$with_umfpack_include_path ";

		AC_CHECK_FILE([${with_umfpack_include_path}/umfpack.h], [], 
              	AC_MSG_ERROR([Include file umfpack.h not found. Try downloading the SuiteSparse (UMFPACK) libraries from http://faculty.cse.tamu.edu/davis/suitesparse.html])) ])

AC_ARG_WITH([umfpack-lib], [AS_HELP_STRING([--with-umfpack-lib=DIR], 
            [path to the SuiteSparse (UMFPACK) libraries])], 
            [ with_umfpack_lib_path=$with_umfpack_lib; USE_UMFPACK=" -DOPENSMOKE_USE_UMFPACK=1 "; UMFPACK_LIBS=" -L$with_umfpack_lib_path -lklu -lumfpack -lamd -lcholmod -lccolamd -lcamd -lsuitesparseconfig";

		AC_CHECK_FILE([${with_umfpack_lib_path}/libumfpack.so], [], 
              	AC_MSG_ERROR([Library file ${with_umfpack_lib_path}/libumfpack.so not found. Try downloading the SuiteSparse (UMFPACK) libraries from http://faculty.cse.tamu.edu/davis/suitesparse.html])) ])


# ---------------------------------------------------------------------------------------------------------------------------------------
# Add SuperLU(serial) library (optional)
# ---------------------------------------------------------------------------------------------------------------------------------------

USE_SUPERLU=" -DOPENSMOKE_USE_SUPERLU=0 "
SUPERLU_INCLUDE_DIR=""
SUPERLU_LIBS=""

AC_ARG_WITH([superlu-include], [AS_HELP_STRING([--with-superlu-include=DIR], 
            [installation prefix (include) for SuperLU (serial) libraries])], 
            [ with_superlu_include_path=$with_superlu_include; USE_SUPERLU=" -DOPENSMOKE_USE_SUPERLU=1 "; SUPERLU_INCLUDE_DIR=" -I$with_superlu_include_path ";

		AC_CHECK_FILE([${with_superlu_include_path}/slu_cdefs.h], [], 
              	AC_MSG_ERROR([Include file slu_cdefs.h not found. Try downloading the SuperLU libraries from http://crd-legacy.lbl.gov/~xiaoye/SuperLU/])) ])

AC_ARG_WITH([superlu-lib], [AS_HELP_STRING([--with-superlu-lib=PATH], 
            [path to the BzzMath libraries])], 
            [ with_superlu_lib_path=$with_superlu_lib; USE_SUPERLU=" -DOPENSMOKE_USE_SUPERLU=1 "; SUPERLU_LIBS=" $with_superlu_lib_path ";

		AC_CHECK_FILE([${with_superlu_lib_path}], [], 
              	AC_MSG_ERROR([Library file ${with_superlu_lib_path} not found. Try downloading the SuperLU libraries from http://crd-legacy.lbl.gov/~xiaoye/SuperLU/])) ])


# ---------------------------------------------------------------------------------------------------------------------------------------
# Add ODEPACK library (optional)
# ---------------------------------------------------------------------------------------------------------------------------------------

USE_ODEPACK=" -DOPENSMOKE_USE_ODEPACK=0 "
ODEPACK_LIBS=""

AC_ARG_WITH([odepack], [AS_HELP_STRING([--with-odepack=DIR], 
            [installation prefix (include) for ODEPACK library])], 
            [ with_odepack_path=$with_odepack; USE_ODEPACK=" -DOPENSMOKE_USE_ODEPACK=1 "; ODEPACK_LIBS=" $with_odepack_path ";

		AC_CHECK_FILE([${with_odepack_path}], [], 
              	AC_MSG_ERROR([Library file ${with_odepack_path} not found. Try downloading the ODEPACK library from https://computation.llnl.gov/casc/odepack/odepack_home.html])) ])


# ---------------------------------------------------------------------------------------------------------------------------------------
# Add DASPK library (optional)
# ---------------------------------------------------------------------------------------------------------------------------------------

USE_DASPK=" -DOPENSMOKE_USE_DASPK=0 "
DASPK_INCLUDE_DIR=""
DASPK_LIBS=""

AC_ARG_WITH([daspk], [AS_HELP_STRING([--with-daspk=DIR], 
            [installation prefix (include) for DASPK library])], 
            [ with_daspk_path=$with_daspk; USE_DASPK=" -DOPENSMOKE_USE_DASPK=1 "; DASPK_LIBS=" $with_daspk_path ";

		AC_CHECK_FILE([${with_daspk_path}], [], 
              	AC_MSG_ERROR([Library file ${with_daspk_path} not found. Try downloading the DASPK library from http://www.engineering.ucsb.edu/~cse/software.html])) ])


# ---------------------------------------------------------------------------------------------------------------------------------------
# Add RADAU library (optional)
# ---------------------------------------------------------------------------------------------------------------------------------------

USE_RADAU=" -DOPENSMOKE_USE_RADAU=0 "
RADAU_LIBS=""

AC_ARG_WITH([radau], [AS_HELP_STRING([--with-radau=DIR], 
            [installation prefix (include) for RADAU library])], 
            [ with_radau_path=$with_radau; USE_RADAU=" -DOPENSMOKE_USE_RADAU=1 "; RADAU_LIBS=" $with_radau_path ";

		AC_CHECK_FILE([${with_radau_path}], [], 
              	AC_MSG_ERROR([Library file ${with_radau_path} not found. Try downloading the RADAU library from https://www.unige.ch/~hairer/software.html])) ])


# ---------------------------------------------------------------------------------------------------------------------------------------
# Add MEBDF library (optional)
# ---------------------------------------------------------------------------------------------------------------------------------------

USE_MEBDF=" -DOPENSMOKE_USE_MEBDF=0 "
MEBDF_LIBS=""

AC_ARG_WITH([mebdf], [AS_HELP_STRING([--with-mebdf=DIR], 
            [installation prefix (include) for MEBDF library])], 
            [ with_mebdf_path=$with_mebdf; USE_MEBDF=" -DOPENSMOKE_USE_MEBDF=1 "; MEBDF_LIBS=" $with_mebdf_path ";

		AC_CHECK_FILE([${with_mebdf_path}], [], 
              	AC_MSG_ERROR([Library file ${with_mebdf_path} not found. Try downloading the MEBDF library from http://wwwf.imperial.ac.uk/~jcash/IVP_software/readme.html])) ])


# ---------------------------------------------------------------------------------------------------------------------------------------
# Add LINPACK library (optional)
# ---------------------------------------------------------------------------------------------------------------------------------------

LINPACK_LIBS=""

AC_ARG_WITH([linpack], [AS_HELP_STRING([--with-linpack=DIR], 
            [installation prefix (include) for LINPACK library])], 
            [ with_linpack_path=$with_linpack; LINPACK_LIBS=" $with_linpack_path ";

		AC_CHECK_FILE([${with_linpack_path}], [], 
              	AC_MSG_ERROR([Library file ${with_linpack_path} not found. Try downloading the LINPACK library from http://www.netlib.org/linpack/])) ])

# ---------------------------------------------------------------------------------------------------------------------------------------
# Add CEQ library (optional, only for thermodynamic equilibrium calculations)
# ---------------------------------------------------------------------------------------------------------------------------------------

CEQ_LIBS=""

AC_ARG_WITH([ceq], [AS_HELP_STRING([--with-ceq=DIR], 
            [installation prefix (include) for CEQ library])], 
            [ with_ceq_path=$with_ceq; CEQ_LIBS=" $with_ceq_path ";

		AC_CHECK_FILE([${with_ceq_path}], [], 
              	AC_MSG_ERROR([Library file ${with_ceq_path} not found. Try downloading the CEQ library from http://eccentric.mae.cornell.edu/~pope/CEQ/])) ])

# ---------------------------------------------------------------------------------------------------------------------------------------
# Add LIS library (optional)
# ---------------------------------------------------------------------------------------------------------------------------------------
LIS_LIBS=""

# ---------------------------------------------------------------------------------------------------------------------------------------
# Add METIS library (optional)
# ---------------------------------------------------------------------------------------------------------------------------------------
METIS_LIBS=""

# ---------------------------------------------------------------------------------------------------------------------------------------
# Add PLASMA library (optional)
# ---------------------------------------------------------------------------------------------------------------------------------------
PLASMA_LIBS=""


#Compulsory libraries
AC_SUBST(OPENSMOKEPP_INCLUDE_DIR)
AC_SUBST(EIGEN_INCLUDE_DIR)
AC_SUBST(RAPIDXML_INCLUDE_DIR)
AC_SUBST(BOOST_INCLUDE_DIR)
AC_SUBST(BOOST_LIBS)

#Intel MKL
AC_SUBST(MKL_INCLUDE_DIR)
AC_SUBST(MKL_LIBS)
AC_SUBST(USE_MKL)

#Intel compiler
AC_SUBST(INTEL_COMPILER_LIBS)

#OpenBLAS
AC_SUBST(OPENBLAS_INCLUDE_DIR)
AC_SUBST(OPENBLAS_LIBS)
AC_SUBST(USE_OPENBLAS)

#Linpack
AC_SUBST(LINPACK_LIBS)

#BzzMath libraries
AC_SUBST(BZZMATH_INCLUDE_DIR)
AC_SUBST(BZZMATH_LIBS)
AC_SUBST(USE_BZZMATH)

#UMFPACK libraries
AC_SUBST(UMFPACK_INCLUDE_DIR)
AC_SUBST(UMFPACK_LIBS)
AC_SUBST(USE_UMFPACK)

#SUPERLU libraries
AC_SUBST(SUPERLU_INCLUDE_DIR)
AC_SUBST(SUPERLU_LIBS)
AC_SUBST(USE_SUPERLU)

#SUNDIALS libraries
AC_SUBST(SUNDIALS_INCLUDE_DIR)
AC_SUBST(SUNDIALS_LIBS)
AC_SUBST(USE_SUNDIALS)

#DASPK
AC_SUBST(DASPK_LIBS)
AC_SUBST(USE_DASPK)

#RADAU
AC_SUBST(RADAU_LIBS)
AC_SUBST(USE_RADAU)

#MEBDF
AC_SUBST(MEBDF_LIBS)
AC_SUBST(USE_MEBDF)

#ODEPACK
AC_SUBST(ODEPACK_LIBS)
AC_SUBST(USE_ODEPACK)

#CEQ
AC_SUBST(CEQ_LIBS)

#OpenBLAS
AC_SUBST(OPENBLAS_LIBS)

#LIS
AC_SUBST(LIS_LIBS)

#Metis
AC_SUBST(METIS_LIBS)

#PLASMA
AC_SUBST(PLASMA_LIBS)


USE_EXTERNAL_LIBRARIES="$USE_BZZMATH $USE_DASPK $USE_RADAU $USE_MEBDF $USE_ODEPACK $USE_UMFPACK $USE_SUPERLU $USE_SUNDIALS $USE_OPENBLAS"
AC_SUBST(USE_EXTERNAL_LIBRARIES)

COMPULSORY_INCLUDE_DIRS="$OPENSMOKEPP_INCLUDE_DIR $EIGEN_INCLUDE_DIR $RAPIDXML_INCLUDE_DIR $BOOST_INCLUDE_DIR"
AC_SUBST(COMPULSORY_INCLUDE_DIRS)

EXTERNAL_INCLUDE_DIRS="$BZZMATH_INCLUDE_DIR $UMFPACK_INCLUDE_DIR $SUPERLU_INCLUDE_DIR $SUNDIALS_INCLUDE_DIR $OPENBLAS_INCLUDE_DIR"
AC_SUBST(EXTERNAL_INCLUDE_DIRS)

EXTERNAL_LIBS="$CEQ_LIBS $BZZMATH_LIBS $DASPK_LIBS $MEBDF_LIBS $ODEPACK_LIBS $LINPACK_LIBS $RADAU_LIBS $SUNDIALS_LIBS $UMFPACK_LIBS $METIS_LIBS $SUPERLU_LIBS $LIS_LIBS $PLASMA_LIBS $OPENBLAS_LIBS $INTEL_COMPILER_LIBS"
AC_SUBST(EXTERNAL_LIBS)


#
AM_INIT_AUTOMAKE([subdir-objects])
AC_CONFIG_FILES([Makefile])

AC_OUTPUT

