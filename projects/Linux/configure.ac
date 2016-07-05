
AC_INIT([CVI], [0.1], [alberto.cuoci@polimi.it])

AM_INIT_AUTOMAKE

AC_PROG_CC
AC_PROG_CXX

CXXFLAGS='-O3 -m64 -fPIC -std=c++0x'



# ---------------------------------------------------------------------------------------------------------------------------------------
# Add OpenSMOKE++ library (compulsory)
# ---------------------------------------------------------------------------------------------------------------------------------------
AC_ARG_WITH([opensmoke], [AS_HELP_STRING([--with-opensmoke=DIR], 
            [installation prefix for OpenSMOKE++ libraries])], 
            [with_opensmoke_include_path=$with_opensmoke/source],
	    AC_MSG_ERROR([Please use --with-opensmokepp to provide the path to OpenSMOKE++. Ask alberto.cuoci@polimi.it for the most recent version of OpenSMOKE++]))

AC_CHECK_FILE([${with_opensmoke_include_path}/OpenSMOKEpp], [], 
              AC_MSG_ERROR([Include file OpenSMOKEpp not found. Ask alberto.cuoci@polimi.it for the most recent version of OpenSMOKE++]))

# ---------------------------------------------------------------------------------------------------------------------------------------
# Add Eigen library (compulsory)
# ---------------------------------------------------------------------------------------------------------------------------------------
AC_ARG_WITH([eigen], [AS_HELP_STRING([--with-eigen=DIR], 
            [installation prefix for Eigen libraries])], 
            [with_eigen_include_path=$with_eigen],
	    AC_MSG_ERROR([Please use --with-eigen to provide the path to Eigen. Try downloading the Eigen libraries from http://eigen.tuxfamily.org/index.php?title=Main_Page]))

AC_CHECK_FILE([${with_eigen_include_path}/Eigen/Eigen], [], 
              AC_MSG_ERROR([Include file Eigen not found. Try downloading the Eigen libraries from http://eigen.tuxfamily.org/index.php?title=Main_Page]))

# ---------------------------------------------------------------------------------------------------------------------------------------
# Add RapidXML library (compulsory)
# ---------------------------------------------------------------------------------------------------------------------------------------
AC_ARG_WITH([rapidxml], [AS_HELP_STRING([--with-rapidxml=DIR], 
            [installation prefix for RapidXML libraries])], 
            [with_rapidxml_include_path=$with_rapidxml],
	    AC_MSG_ERROR([Please use --with-rapidxml to provide the path to RapidXML. Try downloading the RapidXML libraries from http://rapidxml.sourceforge.net]))

AC_CHECK_FILE([${with_rapidxml_include_path}/rapidxml.hpp], [], 
              AC_MSG_ERROR([Include file rapidxml.hpp not found. Try downloading the RapidXML libraries from http://rapidxml.sourceforge.net]))

# ---------------------------------------------------------------------------------------------------------------------------------------
# Add Boost C++ library (compulsory)
# ---------------------------------------------------------------------------------------------------------------------------------------
AC_ARG_WITH([boost-include], [AS_HELP_STRING([--with-boost-include=DIR], 
            [installation prefix (include) for Boost C++ libraries])], 
            [with_boost_include_path=$with_boost_include])

AC_ARG_WITH([boost-lib], [AS_HELP_STRING([--with-boost-lib=DIR], 
            [installation prefix (lib) for Boost C++ libraries])], 
            [with_boost_lib_path=$with_boost_lib])

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
# Add OpenBLAS library (optional)
# ---------------------------------------------------------------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------------------------------------------------------------
# Add Sundials library (optional)
# ---------------------------------------------------------------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------------------------------------------------------------
# Add UMFPACK library (optional)
# ---------------------------------------------------------------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------------------------------------------------------------
# Add SuperLU(serial) library (optional)
# ---------------------------------------------------------------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------------------------------------------------------------
# Add LIS library (optional)
# ---------------------------------------------------------------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------------------------------------------------------------
# Add ODEPACK library (optional)
# ---------------------------------------------------------------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------------------------------------------------------------
# Add DASPK library (optional)
# ---------------------------------------------------------------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------------------------------------------------------------
# Add RADAU library (optional)
# ---------------------------------------------------------------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------------------------------------------------------------
# Add MEBDF library (optional)
# ---------------------------------------------------------------------------------------------------------------------------------------

#Compulsory libraries
OPENSMOKEPP_INCLUDE_DIR="-I$with_opensmoke_include_path"
EIGEN_INCLUDE_DIR="-I$with_eigen_include_path"
RAPIDXML_INCLUDE_DIR="-I$with_rapidxml_include_path"
BOOST_INCLUDE_DIR="-I$with_boost_include_path"
BOOST_LIBS="-L$with_boost_lib_path -lboost_date_time -lboost_filesystem -lboost_program_options -lboost_system -lboost_regex -lboost_timer -lboost_chrono"


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

#BzzMath libraries
AC_SUBST(BZZMATH_INCLUDE_DIR)
AC_SUBST(BZZMATH_LIBS)
AC_SUBST(USE_BZZMATH)

#
AM_INIT_AUTOMAKE([subdir-objects])
AC_CONFIG_FILES([Makefile])

AC_OUTPUT

