/*----------------------------------------------------------------------*\
|    ___                   ____  __  __  ___  _  _______                  |
|   / _ \ _ __   ___ _ __ / ___||  \/  |/ _ \| |/ / ____| _     _         |
|  | | | | '_ \ / _ \ '_ \\___ \| |\/| | | | | ' /|  _| _| |_ _| |_       |
|  | |_| | |_) |  __/ | | |___) | |  | | |_| | . \| |__|_   _|_   _|      |
|   \___/| .__/ \___|_| |_|____/|_|  |_|\___/|_|\_\_____||_|   |_|        |
|        |_|                                                              |
|                                                                         |
|   Author: Alberto Cuoci <alberto.cuoci@polimi.it>                       |
|   CRECK Modeling Group <http://creckmodeling.chem.polimi.it>            |
|   Department of Chemistry, Materials and Chemical Engineering           |
|   Politecnico di Milano                                                 |
|   P.zza Leonardo da Vinci 32, 20133 Milano                              |
|                                                                         |
|-------------------------------------------------------------------------|
|                                                                         |
|   This file is part of OpenSMOKE++ Suite.                               |
|                                                                         |
|   Copyright(C) 2016  Alberto Cuoci                                      |
|   Source-code or binary products cannot be resold or distributed        |
|   Non-commercial use only                                               |
|   Cannot modify source-code for any purpose (cannot create              |
|   derivative works)                                                     |
|                                                                         |
\*-----------------------------------------------------------------------*/

#ifndef OpenSMOKE_Reactor2D_H
#define OpenSMOKE_Reactor2D_H

// BzzMath
#if OPENSMOKE_USE_BZZMATH == 1
#include "BzzMath.hpp"
#endif

// OpenSMOKE++ Definitions
#include "OpenSMOKEpp"

// CHEMKIN maps
#include "maps/Maps_CHEMKIN"

// CHEMKIN maps
#include "PorousMedium.h"

// 1D grid
#include "utilities/grids/adaptive/Grid1D.h"

// Numerical parameters
#include "math/native-ode-solvers/parameters/OdeSolver_Parameters.h"
#include "math/native-dae-solvers/parameters/DaeSolver_Parameters.h"

// Utilities
#include "Utilities.h"

namespace CVI
{
	enum GaseousPhase { GASEOUS_PHASE_FROM_PLUG_FLOW, GASEOUS_PHASE_FROM_CFD };
	enum EquationsSet { EQUATIONS_SET_COMPLETE, EQUATIONS_SET_ONLYTEMPERATURE };
	enum WallType { IMPERMEABLE, NOT_IMPERMEABLE };

	//!  A class to solve the reaction-diffusion equations in 1D 
	/*!
	This class provides the tools to solve the reaction-diffusion equations in 1D 
	*/

	class Reactor2D
	{
	public:

		/**
		*@brief Default constructor
		*@param thermodynamicsMap	reference to the thermodynamic map
		*@param kineticsMap			reference to the kinetic map
		*@param transportMap		reference to the transport map
		*@param porousMedium		reference to the porous medium
		*@param porosityDefect		reference to the porosity defect class
		*@param heterogeneousMechanism		reference to the heterogeneous mechanism
		*@param grid_x				reference to 1D grid along the x direction
		*@param grid_y				reference to 1D grid along the y direction
		*@param plugFlowReactor		reference to the plug flow reactor
		*/
		Reactor2D(	OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap,
					OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap,
					OpenSMOKE::TransportPropertiesMap_CHEMKIN& transportMap,
					OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN& thermodynamicsSurfaceMap,
					OpenSMOKE::KineticsMap_Surface_CHEMKIN&	kineticsSurfaceMap,
					CVI::PorousMedium& porousMedium,
					CVI::PorosityDefect& porosityDefect,
					CVI::HeterogeneousMechanism& heterogeneousMechanism,
					CVI::HeterogeneousDetailedMechanism& heterogeneousDetailedMechanism,
					OpenSMOKE::Grid1D& grid_x, OpenSMOKE::Grid1D& grid_y,
					CVI::PlugFlowReactorCoupled& plugFlowReactor,
					const bool detailed_heterogeneous_kinetics,
					const std::vector<bool>& site_non_conservation,
					const std::string gas_dae_species,
					const std::string surface_dae_species,
					const boost::filesystem::path output_folder);

		/**
		*@brief Sets the planar symmetry
		*@param flag if true, planar symmetry is adopted, otherwise cylindrical symmetry
		*/
		void SetPlanarSymmetry(const bool flag);

		/**
		*@brief Sets the conditions along the gas side
		*@param T_gas		gas side temperature [K]
		*@param P_gas		gas side pressure [Pa]
		*@param omega_gas	gas side mass fractions
		*/
		void SetGasSide(const double T_gas, const double P_gas, const std::vector<Eigen::VectorXd>& omega_gas);

		/**
		*@brief Sets the conditions along the gas side
		*@param disk_from_cfd	gas side mass fractions
		*/
		void SetGasSide(const CVI::DiskFromCFD& disk_from_cfd);

		/**
		*@brief Sets the conditions along the gas side
		*@param profile_temperature		gas side temperature temporal profile [K]
		*@param disk_from_cfd			gas side mass fractions
		*/
		void SetGasSide(OpenSMOKE::FixedProfile* profile_temperature, const CVI::DiskFromCFD& disk_from_cfd);

		/**
		*@brief Sets the initial conditions in the porous medium
		*@param path_to_backup_file		path to backup file (it can be empty)
		*@param T_initial				initial temperature [K]
		*@param P_initial				initial pressure [Pa]
		*@param omega_initial			initial mass fractions
		*/
		void SetInitialConditions(const boost::filesystem::path path_to_backup_file, const double T_initial, const double P_initial, const Eigen::VectorXd& omega_initial, const Eigen::VectorXd& Gamma0, const Eigen::VectorXd& Z0);

		/**
		*@brief Sets the total time of integration
		*@param time_total total time of integration [s]
		*/
		void SetTimeTotal(const double time_total);

		/**
		*@brief Sets the interval time for solving successive DAE systems
		*@param time_interval interval time [s]
		*/
		void SetDaeTimeInterval(const double time_interval);

		/**
		*@brief Sets the total maximum time for integrating ODE systems for determining the initial conditions
		*@param time_interval interval time [s]
		*/
		void SetOdeEndTime(const double time_interval);

		/**
		*@brief Sets the interval time for writing Tecplot output
		*@param time_interval interval time [s]
		*/
		void SetTecplotTimeInterval(const double time_interval);

		/**
		*@brief Sets steps video
		*@param steps_video number of steps to update info on the screen
		*/
		void SetStepsVideo(const int steps_video);

		/**
		*@brief Sets steps file
		*@param steps_file number of steps to update info on file
		*/
		void SetStepsFile(const int steps_file);

		/**
		*@brief Sets steps update plug flow
		*@param steps_update_plug_flow number of steps to update the plug flow
		*/
		void SetStepsUpdatePlugFlow(const int steps_update_plug_flow);

		/**
		*@brief Sets a uniform velocity field
		*@param vx x-component of velocity [m/s]
		*@param vy y-component of velocity [m/s]
		*/
		void SetUniformVelocity(const double vx, const double vy);

		/**
		*@brief Sets the output file name where the integral formation rates are written
		*@param filename file name
		*/
		void SetOutputFile(const std::string filename);

		/**
		*@brief Specifies if local adjustments of surface compostion have to be applied in case of initialization from backup files
		*@param readjust_backup if true, local adjustments are carried out
		*/
		void SetReadjustBackup(const bool readjust_backup);

		/**
		*@brief List of impermeable walls
		*@param impermeable_walls list of impermeable walls
		*/
		void SetImpermeableWalls(const std::vector<std::string>& impermeable_walls);

		/**
		*@brief Returns the differential equations
		*@param t current time [s]
		*@param y current solution
		*@param dy current time derivatives
		*/
		void Equations(const double t, const double* y, double* dy);

		/**
		*@brief Returns the differential equations for the complete set
		*@param t current time [s]
		*@param y current solution
		*@param dy current time derivatives
		*/
		void EquationsComplete(const double t, const double* y, double* dy);

		/**
		*@brief Returns the differential equations for temperature equation only
		*@param t current time [s]
		*@param y current solution
		*@param dy current time derivatives
		*/
		void EquationsOnlyTemperature(const double t, const double* y, double* dy);

		/**
		*@brief Solves the reactor equations
		*@param dae_parameters parameters governing the solution of the DAE system
		*@return the returned value is >0 in case of success, otherwise is <0
		*/
		int SolveFromScratch(DaeSMOKE::DaeSolver_Parameters& dae_parameters, OdeSMOKE::OdeSolver_Parameters& ode_parameters);

		/**
		*@brief Prints info on the screen
		*@param t current time [s]
		*@param y current vector of unknowns
		*/
		void Print(const double t, const double* y);

		/**
		*@brief Returns the total number of equations of the DAE system
		*/
		unsigned int NumberOfEquations() const { return ne_; }

		/**
		*@brief Returns the block dimension
		*/
		unsigned int BlockDimensions() const { return block_; }

		/**
		*@brief Returns the upper band size
		*/
		unsigned int UpperBand() const { return band_size_; }

		/**
		*@brief Returns the lower band size
		*/
		unsigned int LowerBand() const { return band_size_; }

		/**
		*@brief Returns if an equation is differential or algebraic
		*@param v the vector containing if an equation is differential (1.) or algebraic (0.)
		*/
		void AlgebraicDifferentialVector(double* v);

		/**
		*@brief Returns the sparsity pattern of the Jacobian matrix
		*@param row  row indices of non-zero elements (0-based)
		*@param cols column indices of non-zero elements (0-based)
		*/
		void SparsityPattern(std::vector<unsigned int>& rows, std::vector<unsigned int>& cols);

		/**
		*@brief Transfer the current vector to the plug flow reactor unknowns
		*@param y vector of variables to be transferred
		*/
		void Recover_Unknowns(const double* y);

		/**
		*@brief Transfer the current derivatives to the provided vector
		*@param y vector where to transfer the derivatives
		*/
		void Recover_Residuals(double* dy);

		/**
		*@brief Current vector of unknowns
		*@param y vector the current unknowns are copied
		*/
		void UnknownsVector(double* v);

		/**
		*@brief Current vector of unknowns (corrected to ensure the physical constraints)
		*@param y vector the current unknowns are copied
		*/
		void CorrectedUnknownsVector(double* v);

		/**
		*@brief Returns the minimum constraints
		*@param v the minimum constraints for each variable
		*/
		void MinimumUnknownsVector(double* v);

		/**
		*@brief Returns the maximum constraints
		*@param v the maximum constraints for each variable
		*/
		void MaximumUnknownsVector(double* v);

		void CorrectDifferentialEquations(double* upv, double* resv);
		void CorrectAlgebraicEquations(double* yp);
		void DiagonalJacobian(const double t, double* y, double* J);
		void DiagonalJacobianForIDA(const double alfa, double* J);

		/**
		*@brief Sets the non conservation of sites equations
		*@param site_non_conservation boolean vector (true: non conservation is allowed)
		*/
		void SetSiteNonConservation(std::vector<bool>& site_non_conservation);

		void SetSurfaceOnTheFlyROPA(OpenSMOKE::SurfaceOnTheFlyROPA* ropa);

		void SetDerivativeMassFractions(const OpenSMOKE::derivative_type value);

		void SetDerivativeEffectiveDiffusivity(const OpenSMOKE::derivative_type value);

		void SetDerivativeBulkDensity(const OpenSMOKE::derivative_type value);

		int OdeEquations(const double t, const Eigen::VectorXd& y, Eigen::VectorXd& dy);

		int OdePrint(const double t, const Eigen::VectorXd& y);

	private:

		/**
		*@brief Allocates vectors and matrices
		*/
		void MemoryAllocation();

		/**
		*@brief Calculates the concentrations, density, molecular weight, transport properties and formation rates of gaseous species
		*/
		void Properties();

		/**
		*@brief Equations of conservation of mass fractions of species
		*@param t current time [s]
		*/
		void SubEquations_MassFractions(const double t);

		/**
		*@brief Equation describing the evolution of porosity
		*/
		void SubEquations_Porosity();

		/**
		*@brief Equations of conservation of surface fractions of species
		*/
		void SubEquations_SurfaceSpeciesFractions();

		/**
		*@brief Equation of conservation of temperature
		*@param t current time [s]
		*/
		void SubEquations_Temperature(const double t);

		/**
		*@brief Boundary conditions for mass fractions
		*/
		void SubEquations_MassFractions_BoundaryConditions(const double t);

		void SubEquations_MassFractions_BoundaryConditions_WestSide(const double t);
		void SubEquations_MassFractions_BoundaryConditions_EastSide(const double t);
		void SubEquations_MassFractions_BoundaryConditions_NorthSide(const double t);
		void SubEquations_MassFractions_BoundaryConditions_SouthSide(const double t);

		/**
		*@brief Boundary conditions for temperature
		*/
		void SubEquations_Temperature_BoundaryConditions(const double t);

		void SubEquations_Temperature_BoundaryConditions_WestSide(const double t);
		void SubEquations_Temperature_BoundaryConditions_EastSide(const double t);
		void SubEquations_Temperature_BoundaryConditions_NorthSide(const double t);
		void SubEquations_Temperature_BoundaryConditions_SouthSide(const double t);

		/**
		*@brief Prints the current solution on a file
		*/
		void PrintSolution(const double t, const std::string name_file);
		void PrintTecplot(const double t, const std::string name_file);
		void PrintTecplotGlobalKinetics(const double t, const std::string name_file);
		void PrintTecplotDetailedKinetics(const double t, const std::string name_file);

		/**
		*@brief Prints the diffusion coefficients of gaseous species on file
		*@param t current time [s]
		*@param name_file name of file where the coefficients will be written
		*/
		void PrintDiffusionCoefficients(const double t, const std::string name_file);

		/**
		*@brief Prints the formation rates of species due to the homogeneous reactions on file
		*@param t current time [s]
		*@param name_file name of file where the formation rates will be written
		*/
		void PrintHomogeneousRates(const double t, const std::string name_file);

		/**
		*@brief Prints the formation rates of species due to the heterogeneous reactions on file
		*@param t current time [s]
		*@param name_file name of file where the formation rates will be written
		*/
		void PrintHeterogeneousRates(const double t, const std::string name_file);

		/**
		*@brief Prints the formation rates of species due to the heterogeneous reactions on file
		*@param t current time [s]
		*@param name_file name of file where the formation rates will be written
		*/
		void PrintGlobalHeterogeneousRates(const double t, const std::string name_file);

		/**
		*@brief Prints the formation rates of species due to the heterogeneous reactions on file
		*@param t current time [s]
		*@param name_file name of file where the formation rates will be written
		*/
		void PrintDetailedHeterogeneousRates(const double t, const std::string name_file);

		/**
		*@brief Prints the source terms due to the densification to be transferred to CFD
		*@param t current time [s]
		*@param name_file name of file where the source terms will be written
		*/
		void PrintIntegralHomogeneousRates(const double t, const std::string name_file);

		/**
		*@brief Solves the system of DAE describing the 1D reactor
		*@param dae_parameters parameters governing the solution of the DAE system
		*@param t0 starting time [s]
		*@param tf final time [s]
		*@return the returned value is >0 in case of success, otherwise is <0
		*/
		int Solve(DaeSMOKE::DaeSolver_Parameters& dae_parameters, const double t0, const double tf);

		/**
		*@brief Sets the algebraic and differential equations
		*/
		void SetAlgebraicDifferentialEquations();

		void PrintROPA(const double t, const std::string name_file);

		double VolumeAveraged(const Eigen::VectorXd& v);

		double VolumeIntegral(const Eigen::VectorXd& v);

		double VolumeStandardDeviation(const double mean, const Eigen::VectorXd& v);

		void PrintLabelMonitoringFile();

		void PrintXMLFile(const std::string file_name, const double t);

		void SetInitialConditionsFromBackupFile(const boost::filesystem::path path_to_backup_file);

		void UpdateTemperatureBoundaryConditions(const double time);

		void UpdateTemperatureField(const double time);

	protected:

		// References
		OpenSMOKE::ThermodynamicsMap_CHEMKIN&			thermodynamicsMap_;			//!< reference to the thermodynamic map
		OpenSMOKE::KineticsMap_CHEMKIN&					kineticsMap_;				//!< reference to the kinetic map
		OpenSMOKE::TransportPropertiesMap_CHEMKIN&		transportMap_;				//!< reference to the transport properties map
		OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN&	thermodynamicsSurfaceMap_;
		OpenSMOKE::KineticsMap_Surface_CHEMKIN&			kineticsSurfaceMap_;
		CVI::PorousMedium&								porousMedium_;				//!< reference to the porous medium
		CVI::PorosityDefect&							porosityDefect_;			//!< reference to the porosity defect
		CVI::HeterogeneousMechanism&					heterogeneousMechanism_;	//!< reference to the heterogeneous mechanism
		CVI::HeterogeneousDetailedMechanism&			heterogeneousDetailedMechanism_;	//!< reference to the heterogeneous detailed mechanism
		OpenSMOKE::Grid1D&								grid_x_;					//!< reference to the 1D grid along the x direction
		OpenSMOKE::Grid1D&								grid_y_;					//!< reference to the 1D grid along the y direction
		CVI::PlugFlowReactorCoupled&					plugFlowReactor_;			//!< reference to the plug flow reactor

		// ROPA
		OpenSMOKE::SurfaceOnTheFlyROPA*					ropa_;
		bool ropa_analysis_;

		// DAE formulation
		bool dae_formulation_;
		bool detailed_heterogeneous_kinetics_;
		double rho_graphite_;
		unsigned int surface_dae_species_index_;
		unsigned int gas_dae_species_index_;

		int i_current;

		// Dimensions
		unsigned int nc_;						//!< total number of gaseous species
		unsigned int nr_;						//!< total number of homogeneous reactions
		unsigned int np_;						//!< total number of grid points
		unsigned int nx_;						//!< number of grid points along the x direction
		unsigned int ny_;						//!< number of grid points along the x direction
		unsigned int ne_;						//!< total number of equations
		unsigned int block_;					//!< block size
		unsigned int band_size_;				//!< band size (upper and lower are the same)

		unsigned int surf_np_;					//!< total number of site phases
		unsigned int surf_nc_;					//!< total number of surface species
		unsigned int surf_nr_;					//!< total number of heterogeneous reactions

		unsigned int bulk_np_;					//!< total number of bulk phases
		unsigned int bulk_nc_;					//!< total number of bulk species

		std::vector<bool> site_non_conservation_;			//!< site non conservation (true: non conservation is allowed)
		std::vector< Eigen::VectorXd >		Gamma_;			//!< [kmol/m2] (TOCHECK)
		std::vector< Eigen::VectorXd >		GammaFromEqn_;	//!< [kmol/m2] (TOCHECK)


		// Main variables
		Eigen::VectorXd					T_;		//!< current temperature [K]
		Eigen::VectorXd					P_;		//!< current pressure [Pa]
		std::vector<Eigen::VectorXd>	Y_;		//!< mass fractions
		std::vector<Eigen::VectorXd>	X_;		//!< mole fractions
		std::vector<Eigen::VectorXd>	Z_;		//!< surface fractions

		Eigen::VectorXd eigen_C_;				//!< concentrations of gaseous species [kmol/m3]
		Eigen::VectorXd eigen_R_;				//!< formation rates of gaseous species [kg/m3/s]
		Eigen::VectorXd eigen_Z_;				//!< surface fractions [-]
		Eigen::VectorXd eigen_a_;				//!< activities of bulk species [-]
		Eigen::VectorXd eigen_gamma_;			//!< site surface densities [kmol/m2]

		// Properties
		Eigen::VectorXd					rho_gas_;			//!< density of gaseous phase [kg/m3]
		Eigen::VectorXd					rho_bulk_;			//!< density of bulk phase [kg/m3]
		Eigen::VectorXd					rho_bulk_initial_;	//!< density of bulk phase (initial value) [kg/m3]
		Eigen::VectorXd					mw_;				//!< molecular weight of gaseous phase [kg/kmol]
		Eigen::VectorXd					epsilon_;			//!< porosity of porous medium [-]
		Eigen::VectorXd					permeability_;		//!< permeability of porous medium [m2]
		Eigen::VectorXd					eta_bulk_;			//!< tortuosity for ordinary diffusion [-]
		Eigen::VectorXd					eta_knudsen_;		//!< tortuosity for Knudsen diffusion [-]
		Eigen::VectorXd					eta_viscous_;		//!< tortuosity for viscous flow [-]
		Eigen::VectorXd					Sv_;				//!< available area per unit of volume [1/m]
		Eigen::VectorXd					rp_;				//!< radius of pores [m]

		// Reactions
		std::vector<Eigen::VectorXd>	omega_homogeneous_from_homogeneous_;		//!< formation rates of gaseous species [kg/m3/s] (only contribution from homogeneous reactions)
		std::vector<Eigen::VectorXd>	omega_homogeneous_from_heterogeneous_;		//!< formation rates of gaseous species [kg/m3/s] (only contribution from heterogeneous reactions)
		std::vector<Eigen::VectorXd>	omega_heterogeneous_from_heterogeneous_;	//!< formation rates of surface species [kg/m2/s] (only contribution from heterogeneous reactions)

		Eigen::VectorXd					omega_deposition_per_unit_volume_;			//!< deposition rate [kg/m3/s]
		Eigen::VectorXd					omega_deposition_per_unit_area_;			//!< deposition rate [kg/m2/s]
		Eigen::VectorXd					omega_loss_per_unit_volume_;				//!< loss for the homogeneous phase because of heterogeneous reactions [kg/m3/s]

		// Diffusion
		std::vector<Eigen::VectorXd>	gamma_star_;	//!< mass diffusion coefficients [m2/s]

		// Time derivatives
		std::vector<Eigen::VectorXd>	dY_over_dt_;			//!< time derivatives of mass fractions	[1/s]
		Eigen::VectorXd					dT_over_dt_;			//!< time derivative of temperature [K/s]
		Eigen::VectorXd					depsilon_over_dt_;		//!< time derivative of porosity [1/s]
		std::vector<Eigen::VectorXd>	dZ_over_dt_;			//!< time derivatives of fractions of surface species [1/s]
		std::vector<Eigen::VectorXd>	dGamma_over_dt_;		//!< time derivatives of surface densities [kmol/m2/s]

		// North gas side data
		std::vector<Eigen::VectorXd>	Y_gas_north_side_;		//!< mass fractions along the north gas side
		Eigen::VectorXd					T_gas_north_side_;		//!< temperature along the north gas side [K]
		Eigen::VectorXd					P_gas_north_side_;		//!< pressure along the north gas side [Pa]

		// South gas side data
		std::vector<Eigen::VectorXd>	Y_gas_south_side_;		//!< mass fractions along the south gas side
		Eigen::VectorXd					T_gas_south_side_;		//!< temperature along the south gas side [K]
		Eigen::VectorXd					P_gas_south_side_;		//!< pressure along the south gas side [Pa]

		// East gas side data
		std::vector<Eigen::VectorXd>	Y_gas_east_side_;		//!< mass fractions along the east gas side
		Eigen::VectorXd					T_gas_east_side_;		//!< temperature along the east gas side [K]
		Eigen::VectorXd					P_gas_east_side_;		//!< pressure along the east gas side [Pa]

		// West gas side data
		std::vector<Eigen::VectorXd>	Y_gas_west_side_;		//!< mass fractions along the west gas side
		Eigen::VectorXd					T_gas_west_side_;		//!< temperature along the west gas side [K]
		Eigen::VectorXd					P_gas_west_side_;		//!< pressure along the west gas side [Pa]

		// Algebraic/Differential equations
		std::vector<bool>	id_equations_;				//!< algebraic/differential equations
		Eigen::VectorXi		differential_equations_;	//!< list of differential equations
		Eigen::VectorXi		algebraic_equations_;		//!< list of algebraic equations

		// Set of equations to be solved
		EquationsSet	equations_set_;					//!< current set of equations to be solved
		
		// Output
		double t_old_;							//!< time at the end of the previous step [s]
		unsigned int n_steps_video_;					//!< number of steps for updating info on the screen
		unsigned int count_dae_video_;				//!< counter of steps for updating info on the screen
		unsigned int count_ode_video_;				//!< counter of steps for updating info on the screen
		unsigned int n_steps_file_;					//!< number of steps for updating info on the file
		unsigned int count_file_;					//!< counter of steps for updating info on the file
		unsigned int n_steps_update_plug_flow_;				//!< number of steps for updating plug flow
		unsigned int count_update_plug_flow_;				//!< counter of steps for updating plug flow
		std::ofstream fMonitoring_;					//!< name of file to monitor integral quantities over the time

		boost::filesystem::path output_folder_;						//!< name of output folder
		boost::filesystem::path output_tecplot_folder_;				//!< name of output folder for Tecplot files
		boost::filesystem::path output_plug_flow_folder_;			//!< name of output folder for Tecplot files
		boost::filesystem::path output_matlab_folder_;				//!< name of output folder for Matlab files
		boost::filesystem::path output_diffusion_folder_;			//!< name of output folder for diffusion coefficient files
		boost::filesystem::path output_heterogeneous_folder_;		//!< name of output folder for heterogeneous reaction files
		boost::filesystem::path output_homogeneous_folder_;			//!< name of output folder for homogeneous reaction files
		boost::filesystem::path output_ropa_folder_;				//!< name of output folder for ROPA (on the fly)
		boost::filesystem::path output_backup_folder_;				//!< name of output folder for backup files (only 2D with detailed kinetics)
		boost::filesystem::path output_disks_source_terms_folder_;	//!< name of output folder for source terms of disks (to be imported in CFD)
		std::string output_disk_file_name_;							//!< name of output file where integral formation rates are written

		// Post-processing	
		Eigen::VectorXd					delta_rhobulk_;											//!< total increment of bulk density [kg/m3]
		std::vector<Eigen::VectorXd>	delta_rhobulk_due_to_single_reaction_;					//!< increment of bulk density due to the single reactions [kg/m3]
		std::vector<Eigen::VectorXd>	delta_rhobulk_due_to_single_reaction_over_rhobulk_;		//!< increment of bulk density due to the single reactions (normalized)

		// Time
		double time_total_;
		double dae_time_interval_;
		double time_smoothing_;
		double ode_end_time_;
		double time_starting_point_;
		bool start_from_backup_;
		bool readjust_backup_;

		// Tecplot
		int count_tecplot_;
		double tecplot_time_interval_;

		Eigen::VectorXi list_points_south_;
		Eigen::VectorXi list_points_north_;
		Eigen::VectorXi list_points_east_;
		Eigen::VectorXi list_points_west_;

		GaseousPhase gaseous_phase_;

		// Additional options
		bool planar_symmetry_;	// planar vs cylindrical symmetry
		bool hole_;				// true if the hole is present in the geometry
		double vx_;				// x-component of velocity [m/s]
		double vy_;				// y-component of velocity [m/s]

		// Spatial derivatives
		OpenSMOKE::derivative_type derivative_type_mass_fractions_;
		OpenSMOKE::derivative_type derivative_type_effective_diffusivity_;
		OpenSMOKE::derivative_type derivative_type_bulk_density_;

		// Provisional
		void PrintOnTheScreen(const std::string, const int k, double* v);

		// To be used only if the BzzMath libraries are available
		#if OPENSMOKE_USE_BZZMATH == 1
		BzzDaeSparseObject dae_object_;
		#endif

		bool time_profiles_;
		OpenSMOKE::FixedProfile* profile_temperature_;
		//std::vector<OpenSMOKE::FixedProfile*> profiles_omega_;

		Eigen::VectorXd homogeneous_total_mass_source_;		// [kg]
		Eigen::VectorXd heterogeneous_total_mass_source_;	// [kg]
		Eigen::VectorXd total_mass_exchanged_;				// [kg]
		Eigen::VectorXd total_mass_produced_;				// [kg]
		Eigen::VectorXd total_mass_gas_old_;				// [kg]

		CVI::WallType north_wall_type_;
		CVI::WallType south_wall_type_;
		CVI::WallType east_wall_type_;
		CVI::WallType west_wall_type_;

	};
}

namespace OpenSMOKE
{
	class ODESystem_OpenSMOKE_Reactor2D
	{
	public:

		ODESystem_OpenSMOKE_Reactor2D() {};

		void SetReactor2D(CVI::Reactor2D* reactor2d)
		{
			reactor2d_ = reactor2d;
		}

	protected:

		unsigned int ne_;

		void MemoryAllocation()
		{
		}

		virtual void Equations(const Eigen::VectorXd &Y, const double t, Eigen::VectorXd &DY)
		{
			reactor2d_->OdeEquations(t, Y, DY);
		}

		virtual void Jacobian(const Eigen::VectorXd &Y, const double t, Eigen::MatrixXd &J) {};

		void Print(const double t, const Eigen::VectorXd &Y)
		{
			reactor2d_->OdePrint(t, Y);
		}

	private:

		CVI::Reactor2D* reactor2d_;
	};
}

#include "Reactor2D.hpp"

#endif /* OpenSMOKE_PorousMedium_H */
