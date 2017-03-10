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

#ifndef OpenSMOKE_PlugFlowReactorCoupled_H
#define OpenSMOKE_PlugFlowReactorCoupled_H

// OpenSMOKE++ Definitions
#include "OpenSMOKEpp"

// CHEMKIN maps
#include "maps/Maps_CHEMKIN"

namespace CVI
{
	class LinearProfile
	{
	public:

		void Set(OpenSMOKE::OpenSMOKEVectorDouble& x, OpenSMOKE::OpenSMOKEVectorDouble& y);

		double Get(const double x) const;


	private:

		unsigned int number_of_points_;
		double x0_;
		double xf_;
		double y0_;
		double yf_;

		OpenSMOKE::OpenSMOKEVectorDouble x_;
		OpenSMOKE::OpenSMOKEVectorDouble y_;
		OpenSMOKE::OpenSMOKEVectorDouble m_;
	};

	void LinearProfile::Set(OpenSMOKE::OpenSMOKEVectorDouble& x, OpenSMOKE::OpenSMOKEVectorDouble& y)
	{
		x_ = x;
		y_ = y;
		number_of_points_ = x_.Size();

		x0_ = x[1];
		y0_ = y[1];
		xf_ = x[number_of_points_];
		yf_ = y[number_of_points_];

		ChangeDimensions(number_of_points_-1, &m_, true);

		for(unsigned int i=2;i<=number_of_points_;i++)
			m_[i-1] = (y_[i]-y_[i-1])/(x_[i]-x_[i-1]);
	}

	double LinearProfile::Get(const double x) const
	{
		if (x0_ > 0 && ((x - x0_) / x0_ > -1.e6))
			OpenSMOKE::FatalErrorMessage("Profile class: the required point is outside the domain (too small)");

		if (number_of_points_ == 2)
			return y0_ + m_[1]*(x-x0_);

		for(unsigned int i=2;i<=number_of_points_;i++)
			if (x <= x_[i])
				return y_[i-1] + m_[i-1]*(x-x_[i-1]);

		// In case of small excess
		return y_[number_of_points_];
	}

	//!  A class to solve a homogeneous plug flow reactor
	/*!
	This class provides the tools to solve a homogeneous plug flow reactor
	*/

	class PlugFlowReactorCoupled
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
		PlugFlowReactorCoupled(OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap, OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap, const double v, const double Dh);

		/**
		*@brief Default constructor
		*@param thermodynamicsMap	reference to the thermodynamic map
		*@param kineticsMap			reference to the kinetic map
		*/
		PlugFlowReactorCoupled(OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap, OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap, OpenSMOKE::OpenSMOKE_Dictionary& dictionary);

		/**
		*@brief Sets the inlet/initial conditions
		*@param T		inlet/initial temperature [K]
		*@param P		inlet/initial pressure [Pa]
		*@param omega	inlet/initial composition in mass fractions
		*/
		void SetInitialConditions(const double T, const double P, const Eigen::VectorXd& omega);

		/**
		*@brief Sets coupling
		*@param coupling coupling on/off
		*/
		void SetCoupling(const bool coupling);

		/**
		*@brief Sets verbose output
		*@param flag verbose output on/off
		*/
		void SetVerboseOutput(const bool flag);

		/**
		*@brief Sets the inlet/initial conditions
		*@param csi_external space coordinates [m]
		*@param omega_external composition in mass fractions (points x species)
		*/
		void SetExternalMassFractionsProfile(const Eigen::VectorXd& csi_external, const Eigen::MatrixXd& omega_external);

		/**
		*@brief Sets the length of the inert zone
		*@param inert_length length of the inert zone [m]
		*/
		void SetInertLength(const double inert_length);

		/**
		*@brief Sets the width of the channel
		*@param width of the channel
		*/
		void SetChannelWidth(const double channel_width);

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
		*@brief Prints the plug flow profiles
		*@param file_path the file name
		*@param t the current reactor time
		*/
		void Print(const double t, const boost::filesystem::path file_path);

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
		*@brief Direct access to the current unknowns
		*/
		bool coupling() const { return coupling_; }

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
		*@brief Last residence time simulated [s]
		*/
		double last_residence_time() const { return last_residence_time_; }

		/**
		*@brief Returns the mass transfer coeffcient as a function of the axial coordinate
		*@param T temperature [K]
		*@param P pressure [Pa]
		*@param rho density [kg/m3]
		*@param x axial coordinate [m]
		*@return mass transfer coefficient [m/s]
		*/
		double mass_transfer_coefficient(const double T, const double P_Pa, const double rho, const double x);

		/**
		*@brief Inlet temperature [K]
		*/
		double inlet_temperature() const { return inlet_temperature_; }

		/**
		*@brief Inlet pressure [Pa]
		*/
		double inlet_pressure() const { return inlet_pressure_; }

		/**
		*@brief Inlet mass fractions [-]
		*/
		const Eigen::VectorXd& inlet_mass_fractions() const { return inlet_omega_; }

		/**
		*@brief Inlet mass fractions [-]
		*/
		bool verbose_output() const { return verbose_output_; }
		
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

		OpenSMOKE::ThermodynamicsMap_CHEMKIN&			thermodynamicsMap_;	//!< reference to the thermodynamic map
		OpenSMOKE::KineticsMap_CHEMKIN&					kineticsMap_;		//!< reference to the kinetic map

		bool coupling_;									//!< coupling with the carbon felt
		double T_;										//!< current temperature [K]
		double P_;										//!< current pressure [Pa]
		double csi_;									//!< axial coordinate [m]
		double v_;										//!< axial velocity [m/s]
		double Dh_;										//!< hydraulic diameter [m]
		double NuInf_;									//!< asymptotic Nusselt number [-]
		double inert_length_;							//!< length of inert zone [m]
		double channel_width_;							//!< width of the channel [m]
		bool internal_boundary_layer_correction_;		//!< true if the boundary layer limitations have to be accounted for
		GeometricPattern geometric_pattern_;			//!< geometric pattern: ONE_SIDE | THREE_SIDES
		double last_residence_time_;					//!< last residence time simulated [s]

		double inlet_temperature_;
		double inlet_pressure_;
		Eigen::VectorXd inlet_omega_;

		double rho_;					//!< current density [kg/m3]
		Eigen::VectorXd Y_;				//!< current mass fractions [-]
		Eigen::VectorXd dY_over_dt_;	//!< current time derivatives of mass fractions [1/s]
		double dcsi_over_dt_;			//!< current time derivative of space [m/s]

		std::vector<LinearProfile>	Y_external_;		//!< external mass fraction profiles 

		unsigned int n_steps_video_;	//!< number of steps for updating info on the screen
		unsigned int count_video_;		//!< counter of steps for updating info on the screen

		unsigned int ns_;						//!< total number of gaseous species
		unsigned int ne_;						//!< total number of equations
		boost::filesystem::path output_folder_;				//!< name of output folder
		bool verbose_output_;

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

#include "PlugFlowReactorCoupled.hpp"

#endif /* OpenSMOKE_PlugFlowReactorCoupled_H */
