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

		/**
		*@brief Default constructor
		*@param thermodynamicsMap	reference to the thermodynamic map
		*@param kineticsMap			reference to the kinetic map
		*/
		PlugFlowReactor(OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>& thermodynamicsMap, OpenSMOKE::KineticsMap_CHEMKIN<double>& kineticsMap);

		/**
		*@brief Sets the inlet/initial conditions
		*@param T		inlet/initial temperature [K]
		*@param P		inlet/initial pressure [Pa]
		*@param omega	inlet/initial composition in mass fractions
		*/
		void SetInitialConditions(const double T, const double P, const OpenSMOKE::OpenSMOKEVectorDouble& omega);

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
		*@brief Prints the current solution on a file
		*/
		void PrintSolution(const std::string name_file);

	private:

		OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>&			thermodynamicsMap_;	//!< reference to the thermodynamic map
		OpenSMOKE::KineticsMap_CHEMKIN<double>&					kineticsMap_;		//!< reference to the kinetic map

		double T_;						//!< current temperature [K]
		double P_;						//!< current pressure [Pa]
		double rho_;					//!< current density [kg/m3]
		Eigen::VectorXd Y_;				//!< current mass fractions [-]
		Eigen::VectorXd dY_over_dt_;	//!< current time derivatives of mass fractions [1/s]

		unsigned int n_steps_video_;	//!< number of steps for updating info on the screen
		unsigned int count_video_;		//!< counter of steps for updating info on the screen

		unsigned int ns_;						//!< total number of gaseous species
		unsigned int ne_;						//!< total number of equations
		boost::filesystem::path output_folder_;	//!< name of output folder

		// Auxiliary vectors
		OpenSMOKE::OpenSMOKEVectorDouble	aux_Y;	//!< vector containing the mass fractions
		OpenSMOKE::OpenSMOKEVectorDouble	aux_X;	//!< vector containing the mole fractions
		OpenSMOKE::OpenSMOKEVectorDouble	aux_C;	//!< vector containing the concentration of gaseous species [kmol/m3]
		OpenSMOKE::OpenSMOKEVectorDouble	aux_R;	//!< vector containing the formation rates of gaseous species [kg/m3/s]
	};
}

#include "PlugFlowReactor.hpp"

#endif /* OpenSMOKE_PlugFlowReactor_H */