# OpenSMOKE++ CVI
Chemical Vapour Infiltration with OpenSMOKE++

This is the initial release of OpenSMOKE_CVI solver for numerical modelling of **Chemical Vapour Infiltration** in 1D and 2D geometries. 
The code solves the governing equations of transport of species and porosity through the finite difference method. The equations (corresponding to a structured DAE system) are solved using a fully-coupled method. 

Two types of simulations can be currently carried out:
- diffusion and reaction in a cylindrical capillary
- diffusion and reaction in a rectangular domain (planar and cylindrical symmetry) surrounded by a gaseous phase (1 or 3 sides)

# Main features
- arbitrarily complex kinetic mechanisms in the gaseous phase (based on the standard CHEMKIN format)
- hard-coded heterogeneous mechanisms are implemented (only global reactions)
- coupling with several DAE solvers (OpenSMOKE++, BzzMath, DASPK, IDA from Sundials)
- band, tridiagonal-block and sparse linear algebra (Intel MKL Pardiso, UMFPACK from SuiteSparse library, SuperLU serial)

# Dependencies (compulsory)
- OpenSMOKE++ (https://github.com/acuoci/OpenSMOKE.git)
- Eigen C++ (http://eigen.tuxfamily.org/index.php?title=Main_Page)
- RapidXML (http://rapidxml.sourceforge.net/)
- Boost C++ (http://www.boost.org/)
- Blas/Lapack (Intel MKL implementation is recommended for best performances)

# Dependencies (optional)
- Sundials (http://computation.llnl.gov/projects/sundials-suite-nonlinear-differential-algebraic-equation-solvers)
- BzzMath (http://super.chem.polimi.it/download/bzzmath-download/)
- SuperLU (http://crd-legacy.lbl.gov/~xiaoye/SuperLU/)
- SuiteSparse (http://faculty.cse.tamu.edu/davis/suitesparse.html)

# Future developments and extensions
Future developments will include:
- addition of arbitrarily complex heterogeneous kinetic mechanisms in CHEMKIN format
- addition of axial dispersion to the plug flow reactor (i.e. transformation into a tubular reactor)
- improvements of DAE solvers though the adoption of shared (openmp) libraries for solving linear systems

