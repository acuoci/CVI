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

#ifndef OpenSMOKE_PorousMedium_H
#define OpenSMOKE_PorousMedium_H

// OpenSMOKE++ Definitions
#include "OpenSMOKEpp"

// CHEMKIN maps
#include "maps/Maps_CHEMKIN"

namespace CVI
{
	enum PorousSubstrateType { POLYNOMIAL, RANDOM, RANDOM_HARDCORE, POLINOMIAL_ONEHALF, FROM_SPHERES_TO_CYLINDERS, DEUTSCHMANN_CORRELATION };

	//!  A class to manage properties of porous media
	/*!
	This class provides the tools to manage properties of porous media
	*/

	class PorousMedium
	{
	public:

		/**
		*@brief Default constructor
		*@param thermodynamicsMap		reference to the thermodynamic map
		*@param kineticsMap			reference to the kinetic map
		*@param transportMap			reference to the transport map
		*@param porous_substrate_type		porous substrate type
		*@param rf				initial radius of the fibers [m]
		*@param rho_fiber			fiber density [kg/m3]
		*@param epsilon0			initial porosity
		*/
		PorousMedium(	OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap,
						OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap,
						OpenSMOKE::TransportPropertiesMap_CHEMKIN& transportMap,
						PorousSubstrateType porous_substrate_type, const double rf, const double rho_fiber, const double epsilon0);

		PorousMedium(	OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap,
						OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap,
						OpenSMOKE::TransportPropertiesMap_CHEMKIN& transportMap,
						OpenSMOKE::OpenSMOKE_Dictionary& dictionary);

		/**
		*@brief Sets the temperature at which the properties have to be evaluated
		*@param T temperature [K]
		*/
		void SetTemperature(const double T);

		/**
		*@brief Sets the pressure at which the properties have to be evaluated
		*@param P pressure [Pa]
		*/
		void SetPressure(const double P_Pa);

		/**
		*@brief Sets the porosity at which the properties have to be evaluated
		*@param epsilon porosity
		*/
		void SetPorosity(const double epsilon);

		/**
		*@brief Calculates and returns the surface per unit of volume
		*@return the surface per unit of volume [1/m]
		*/
		double Sv();

		/**
		*@brief Calculates and returns the radius of pores
		*@return the radius of pores [m]
		*/
		double rp();

		/**
		*@brief Calculates and returns the tortuosity for ordinary diffusion
		*@return the tortuosity for ordinary diffusion [-]
		*/
		double eta_bulk();

		/**
		*@brief Calculates and returns the tortuosity for Knudsen diffusion
		*@return the tortuosity for Knudsen diffusion [-]
		*/
		double eta_knudsen();

		/**
		*@brief Calculates and returns the tortuosity for viscous flows
		*@return the tortuosity for viscous flows [-]
		*/
		double eta_viscous();

		/**
		*@brief Calculates and returns the permeability
		*@return the permeability [m2]
		*/
		double permeability();

		/**
		*@brief Calculates and returns the bulk density
		*@return the bulk density [kg/m3]
		*/
		double density_bulk(const double rho_graphite);

		double epsilon_threshold() { return epsilon_threshold_; }

		/**
		*@brief Calculates the effective mass diffusion coefficients resulting from the combination of ordinary and Knudsen diffusion
		*@param mole_fractions mole fractions of gaseous species
		*/
		void EffectiveDiffusionCoefficients(const Eigen::VectorXd& mole_fractions);

		/**
		*@brief Returns the effective mass diffusion coefficients accounting for ordinary and Knudsen diffusion
		*@return the effective mass diffusion coefficients accounting for ordinary and Knudsen diffusion [m2/s]
		*/
		const Eigen::VectorXd& gamma_effective() const { return gamma_effective_; }

		/**
		*@brief Returns the effective mass diffusion coefficients accounting only for ordinary diffusion
		*@return the effective mass diffusion coefficients accounting only for ordinary diffusion [m2/s]
		*/
		const Eigen::VectorXd& gamma_fick_effective() const { return gamma_fick_effective_; }

		/**
		*@brief Returns the effective mass diffusion coefficients accounting only for Knudsen diffusion
		*@return the effective mass diffusion coefficients accounting only for Knudseny diffusion [m2/s]
		*/
		const Eigen::VectorXd& gamma_knudsen_effective() const { return gamma_knudsen_effective_; }

		/**
		*@brief Returns the non-effective mass diffusion coefficients accounting only for ordinary diffusion
		*@return the non-effective mass diffusion coefficients accounting only for ordinary diffusion [m2/s]
		*/
		const Eigen::VectorXd& gamma_fick() const { return gamma_fick_; }

		/**
		*@brief Returns the non-effective mass diffusion coefficients accounting only for Knudsen diffusion
		*@return the non-effective mass diffusion coefficients accounting only for Knudseny diffusion [m2/s]
		*/
		const Eigen::VectorXd& gamma_knudsen() const { return gamma_knudsen_; }

		/**
		*@brief Returns the porosity
		*@return the porosity [-]
		*/
		double porosity() const { return epsilon_; }

	private:

		/**
		*@brief Calculates the Knudsen mass diffusion coefficients
		*/
		void KnudsenDiffusionCoefficients();

		/**
		*@brief Calculates the Knudsen effective mass diffusion coefficients
		*/
		void KnudsenEffectiveDiffusionCoefficients();

		/**
		*@brief Calculates the ordinary mass diffusion coefficients
		*/
		void FickDiffusionCoefficients(const Eigen::VectorXd& mole_fractions);

		/**
		*@brief Calculates the ordinary effective mass diffusion coefficients
		*/
		void FickEffectiveDiffusionCoefficients(const Eigen::VectorXd& mole_fractions);

		/**
		*@brief Prints useful data on the screen
		*/
		void Summary();

		void Initialize();

	private:

		OpenSMOKE::ThermodynamicsMap_CHEMKIN&			thermodynamicsMap_;			//!< reference to the thermodynamic map
		OpenSMOKE::KineticsMap_CHEMKIN&					kineticsMap_;				//!< reference to the kinetic map
		OpenSMOKE::TransportPropertiesMap_CHEMKIN&		transportMap_;				//!< reference to the trasport properties map

		unsigned int ns_;

		double T_;				//!< current temperature [K]
		double P_Pa_;			//!< current pressure [Pa]
		double epsilon_;		//!< current porosity
		double epsilon0_;		//!< initial porosity
		double rf_;				//!< current radius of fibers [m]
		double rho_fiber_;		//!< density of fibers [kg/m3]

		double mass_diffusion_multiplier_;

		PorousSubstrateType porous_substrate_type_;		//!< type of substrate

		Eigen::VectorXd gamma_knudsen_;					//!< non-effective mass diffusion coefficients (Knudsen) [m2/s]
		Eigen::VectorXd gamma_knudsen_effective_;		//!< effective mass diffusion coefficients (Knudsen) [m2/s]
		Eigen::VectorXd gamma_fick_;					//!< non-effective mass diffusion coefficients (ordinary or Fick) [m2/s]
		Eigen::VectorXd gamma_fick_effective_;			//!< effective mass diffusion coefficients (ordinary or Fick) [m2/s]
		Eigen::VectorXd gamma_effective_;				//!< effective mass diffusion coefficients (ordinary + Knudsen) [m2/s]

		double epsilon_threshold_;
	};
}

#include "PorousMedium.hpp"

#endif /* OpenSMOKE_PorousMedium_H */
