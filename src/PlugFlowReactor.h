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

#ifndef OpenSMOKE_PlugFlowReactor_H
#define OpenSMOKE_PlugFlowReactor_H

// OpenSMOKE++ Definitions
#include "OpenSMOKEpp"

// CHEMKIN maps
#include "maps/Maps_CHEMKIN"

namespace CVI
{
	//!  A class to solve a homogeneous plug flow reactor
	/*!
	This class provides the tools to solve a homogeneous plug flow reactor
	*/

	class PlugFlowReactor
	{

	public:

		enum GeometricPattern {ONE_SIDE, THREE_SIDES};

	public:

		/**
		*@brief Default constructor
		*@param thermodynamicsMap	reference to the thermodynamic map
		*@param kineticsMap			reference to the kinetic map
		*@param v					the axial velocity [m/s]
		*@param Dh					the hydraulic diameter [m]
		*/
		PlugFlowReactor(OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>& thermodynamicsMap, OpenSMOKE::KineticsMap_CHEMKIN<double>& kineticsMap, const double v, const double Dh);

		/**
		*@brief Default constructor
		*@param thermodynamicsMap	reference to the thermodynamic map
		*@param kineticsMap			reference to the kinetic map
		*/
		PlugFlowReactor(OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>& thermodynamicsMap, OpenSMOKE::KineticsMap_CHEMKIN<double>& kineticsMap, OpenSMOKE::OpenSMOKE_Dictionary& dictionary);

		/**
		*@brief Sets the inlet/initial conditions
		*@param T		inlet/initial temperature [K]
		*@param P		inlet/initial pressure [Pa]
		*@param omega	inlet/initial composition in mass fractions
		*/
		void SetInitialConditions(const double T, const double P, const Eigen::VectorXd& omega);

		/**
		*@brief Sets the length of the inert zone
		*@param inert_length length of the inert zone [m]
		*/
		void SetInertLength(const double inert_length);

		/**
		*@brief Sets the asymptoti Nusselt number
		*@param NuInf the asymptoti Nusselt number
		*/
		void SetAsymptoticNusseltNumber(const double NuInf);

		/**
		*@brief Sets the internal boundary layer correction (on/off)
		*@param flag true if the internal boundary layer limitations have to be accounted for
		*/
		void SetInternalBoundaryLayerCorrection(const bool flag);

		/**
		*@brief Sets the geometric pattern: ONE_SIDE | THREE SIDES
		*@param pattern the geometric pattern
		*/
		void SetGeometricPattern(const GeometricPattern pattern);

		/**
		*@brief Returns the differential equations
		*@param t current time [s]
		*@param y current solution
		*@param dy current time derivatives
		*/
		void Equations(const double t, const double* y, double* dy);

		/**
		*@brief Solves the plug flow reactor equations
		*@param tau residence time [s]
		*@return true if solution was found successfully
		*/
		bool Solve(const double tau);

		/**
		*@brief Prints info on the screen
		*@param t current time [s]
		*@param y current vector of unknowns
		*/
		void Print(const double t, const double* y);

		/**
		*@brief Returns the total number of equations of the ODE system
		*/
		unsigned int NumberOfEquations() const { return ne_; }

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
		*@brief Direct access to the current unknowns
		*/
		const Eigen::VectorXd& Y() const { return Y_; }

		/**
		*@brief Direct access to the history of tau
		*/
		const std::vector<double>& history_tau() const { return history_tau_; }

		/**
		*@brief Direct access to the history of csi
		*/
		const std::vector<double>& history_csi() const { return history_csi_; }

		/**
		*@brief Direct access to the history of mass fractions
		*/
		const std::vector<Eigen::VectorXd> history_Y() const { return history_Y_; }

		/**
		*@brief Returns the axial velocity [m/s]
		*/
		double v() const { return v_; }

		/**
		*@brief Returns the hydraulic diameter [m]
		*/
		double Dh() const { return Dh_; }

		/**
		*@brief Returns the asymptotic Nusselt number [-]
		*/
		double NuInf() const { return NuInf_; }

		/**
		*@brief Returns the length of the inert zone [m]
		*/
		double inert_length() const { return inert_length_; }

		/**
		*@brief Returns true if the internal boundary layer limitations have to be applied
		*/
		bool internal_boundary_layer_correction() const { return internal_boundary_layer_correction_; }

		/**
		*@brief Returns the geometric pattern
		*/
		GeometricPattern geometric_pattern() const { return geometric_pattern_; }

		/**
		*@brief Returns the mass transfer coeffcient as a function of the axial coordinate
		*@param T temperature [K]
		*@param P pressure [Pa]
		*@param rho density [kg/m3]
		*@param x axial coordinate [m]
		*@return mass transfer coefficient [m/s]
		*/
		double mass_transfer_coefficient(const double T, const double P_Pa, const double rho, const double x);
		
	private:

		/**
		*@brief Allocates vectors and matrices
		*/
		void MemoryAllocation();

		/**
		*@brief Calculates the concentrations, density, molecular weight and formation rates of gaseous species
		*/
		void Properties();

		/**
		*@brief Equations of conservation of mass fractions of species
		*/
		void SubEquations_MassFractions();

		/**
		*@brief Equations of space
		*/
		void SubEquations_Space();

		/**
		*@brief Prints the current solution on a file
		*/
		void PrintSolution(const std::string name_file);

		/**
		*@brief Sets the default values for relevant variables
		*/
		void DefaultValues();

		/**
		*@brief Initializes the variables
		*/
		void Initialize();

	private:

		OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>&			thermodynamicsMap_;	//!< reference to the thermodynamic map
		OpenSMOKE::KineticsMap_CHEMKIN<double>&					kineticsMap_;		//!< reference to the kinetic map

		double T_;						//!< current temperature [K]
		double P_;						//!< current pressure [Pa]
		double csi_;					//!< axial coordinate [m]
		double v_;						//!< axial velocity [m/s]
		double Dh_;						//!< hydraulic diameter [m]
		double NuInf_;					//!< asymptotic Nusselt number [-]
		double inert_length_;			//!< length of inert zone [m]
		bool internal_boundary_layer_correction_;		//!< true if the boundary layer limitations have to be accounted for
		GeometricPattern geometric_pattern_;			//!< geometric pattern: ONE_SIDE | THREE_SIDES

		double rho_;					//!< current density [kg/m3]
		Eigen::VectorXd Y_;				//!< current mass fractions [-]
		Eigen::VectorXd dY_over_dt_;	//!< current time derivatives of mass fractions [1/s]
		double dcsi_over_dt_;			//!< current time derivative of space [m/s]

		unsigned int n_steps_video_;	//!< number of steps for updating info on the screen
		unsigned int count_video_;		//!< counter of steps for updating info on the screen

		unsigned int ns_;						//!< total number of gaseous species
		unsigned int ne_;						//!< total number of equations
		boost::filesystem::path output_folder_;	//!< name of output folder

		// History
		std::vector<double>				history_tau_;
		std::vector<double>				history_csi_;
		std::vector<Eigen::VectorXd>	history_Y_;

		// Auxiliary vectors
		OpenSMOKE::OpenSMOKEVectorDouble	aux_Y;	//!< vector containing the mass fractions
		OpenSMOKE::OpenSMOKEVectorDouble	aux_X;	//!< vector containing the mole fractions
		OpenSMOKE::OpenSMOKEVectorDouble	aux_C;	//!< vector containing the concentration of gaseous species [kmol/m3]
		OpenSMOKE::OpenSMOKEVectorDouble	aux_R;	//!< vector containing the formation rates of gaseous species [kg/m3/s]
	};
}

#include "PlugFlowReactor.hpp"

#endif /* OpenSMOKE_PlugFlowReactor_H */