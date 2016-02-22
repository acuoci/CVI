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
	enum PorousSubstrateType { POLYNOMIAL, RANDOM, RANDOM_HARDCORE, POLINOMIAL_ONEHALF, FROM_SPHERES_TO_CYLINDERS };
	enum HydrogenInhibitionType { NONE, BECKER };
	enum HeterogeneousMechanism { IBRAHIM_PAOLUCCI };

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
		*@param kineticsMap				reference to the kinetic map
		*@param transportMap			reference to the transport map
		*@param porous_substrate_type	porous substrate type
		*@param rf						initial radius of the fibers [m]
		*@param epsilon0						initial porosity
		*@param heterogeneous_mechanism_type	heterogeneous mechanism type
		*@param hydrogen_inhibition_type		hydrogen mechanism type
		*/
		PorousMedium(	OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>& thermodynamicsMap,
						OpenSMOKE::KineticsMap_CHEMKIN<double>& kineticsMap,
						OpenSMOKE::TransportPropertiesMap_CHEMKIN<double>& transportMap,
						PorousSubstrateType porous_substrate_type, const double rf, const double epsilon0, 
						const HeterogeneousMechanism heterogeneous_mechanism_type,
						const HydrogenInhibitionType hydrogen_inhibition_type);

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
		*@brief Sets the density of carbon matrix
		*@param rho_carbon density [kg/m3]
		*/
		void SetCarbonDensity(const double rho_carbon);

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
		*@brief Returns the density of the carbon matrix
		*@return the density of the carbon matrix [kg/m3]
		*/
		double rho_carbon() const { return rho_carbon_; }

		/**
		*@brief Returns the molecular weight of the carbon matrix
		*@return the molecular weight of the carbon matrix [kg/kmol]
		*/
		double mw_carbon() const { return mw_carbon_; }

		/**
		*@brief Calculates the effective mass diffusion coefficients resulting from the combination of ordinary and Knudsen diffusion
		*@param mole_fractions mole fractions of gaseous species
		*/
		void EffectiveDiffusionCoefficients(const Eigen::VectorXd& mole_fractions);

		/**
		*@brief Calculates the formation rates of species due to the heterogeneous reactions
		*@param C concentrations of gaseous species [kmol/m3]
		*/
		void FormationRates(const Eigen::VectorXd& C);

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

		/**
		*@brief Returns the formation rates of gaseous species due to the heterogeneous reactions
		*@return the formation rates [kmol/m3/s]
		*/
		const Eigen::VectorXd& Rgas() const { return Rgas_; }

		/**
		*@brief Returns the reaction rates of heterogeneous reactions
		*@return the reaction rates [kmol/m3/s]
		*/
		const Eigen::VectorXd& r() const { return r_; }

		/**
		*@brief Returns the deposition rate per unit of area
		*@return the deposition rate [kmol/m2/s]
		*/
		double r_deposition_per_unit_area() const { return r_deposition_per_unit_area_; }

		/**
		*@brief Returns the deposition rate per unit of volume
		*@return the deposition rate [kmol/m3/s]
		*/
		double r_deposition_per_unit_volume() const { return r_deposition_per_unit_volume_; }

		/**
		*@brief Returns the inhibition coefficient for methane
		*@return the inhibition coefficient for methane
		*/
		double I_CH4()  const { return I_CH4_;  }

		/**
		*@brief Returns the inhibition coefficient for ethylene
		*@return the inhibition coefficient for ethylene
		*/
		double I_C2H4() const { return I_C2H4_; }

		/**
		*@brief Returns the inhibition coefficient for acetylene
		*@return the inhibition coefficient for acetylene
		*/
		double I_C2H2() const { return I_C2H2_; }

		/**
		*@brief Returns the inhibition coefficient for benzene
		*@return the inhibition coefficient for benzene
		*/
		double I_C6H6() const { return I_C6H6_; }

		/**
		*@brief Returns the index for methane
		*@return the index (0-based) for methane
		*/
		unsigned int index_CH4() const { return index_CH4_; }

		/**
		*@brief Returns the index for ethylene
		*@return the index (0-based) for ethylene
		*/
		unsigned int index_C2H4() const { return index_C2H4_; }

		/**
		*@brief Returns the index for acetylene
		*@return the index (0-based) for acetylene
		*/
		unsigned int index_C2H2() const { return index_C2H2_; }

		/**
		*@brief Returns the index for benzene
		*@return the index (0-based) for benzene
		*/
		unsigned int index_C6H6() const { return index_C6H6_; }

		/**
		*@brief Returns the index for hydrogen
		*@return the index (0-based) for hydrogen
		*/
		unsigned int index_H2() const { return index_H2_; }

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
		*@brief Calculates the hydrogen inhibition coefficients
		*/
		void HydrogenInhibitionCoefficients(const Eigen::VectorXd& C);

		/**
		*@brief Prints useful data on the screen
		*/
		void Summary();

	private:

		OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>&			thermodynamicsMap_;	//!< reference to the thermodynamic map
		OpenSMOKE::KineticsMap_CHEMKIN<double>&					kineticsMap_;		//!< reference to the kinetic map
		OpenSMOKE::TransportPropertiesMap_CHEMKIN<double>&		transportMap_;		//!< reference to the trasport properties map

		unsigned int ns_;	//!< total number of gaseous species
		unsigned int nr_;	//!< total number of heterogeneous reactions

		int index_CH4_;		//!< index of CH4  (0-based)
		int index_C2H4_;	//!< index of C2H4 (0-based)
		int index_C2H2_;	//!< index of C2H2 (0-based)
		int index_C6H6_;	//!< index of C6H6 (0-based)
		int index_H2_;		//!< index of H2   (0-based)

		double T_;				//!< current temperature [K]
		double P_Pa_;			//!< current pressure [Pa]
		double epsilon_;		//!< current porosity
		double epsilon0_;		//!< initial porosity
		double rf_;				//!< current radius of fibers [m]
		double rho_carbon_;		//!< density of carbon matrix [kg/m3]
		double mw_carbon_;		//!< molecular weight of carbon matrix [kg/kmol]

		PorousSubstrateType porous_substrate_type_;		//!< type of substrate

		Eigen::VectorXd gamma_knudsen_;					//!< non-effective mass diffusion coefficients (Knudsen) [m2/s]
		Eigen::VectorXd gamma_knudsen_effective_;		//!< effective mass diffusion coefficients (Knudsen) [m2/s]
		Eigen::VectorXd gamma_fick_;					//!< non-effective mass diffusion coefficients (ordinary or Fick) [m2/s]
		Eigen::VectorXd gamma_fick_effective_;			//!< effective mass diffusion coefficients (ordinary or Fick) [m2/s]
		Eigen::VectorXd gamma_effective_;				//!< effective mass diffusion coefficients (ordinary + Knudsen) [m2/s]

		HydrogenInhibitionType hydrogen_inhibition_type_;	//!< hydrogene inhibition type
		double I_CH4_;										//!< inhibition coefficient for methane
		double I_C2H4_;										//!< inhibition coefficient for ethylene
		double I_C2H2_;										//!< inhibition coefficient for acetylene
		double I_C6H6_;										//!< inhibition coefficient for benzene

		HeterogeneousMechanism heterogeneous_mechanism_type_;		//!< type of heterogeneous mechanism
		Eigen::VectorXd Rgas_;										//!< formation rate of gaseous species due to heterogeneous reactions [kmol/m3/s]
		Eigen::VectorXd r_;											//!< reaction rates of heterogeneous reactions [kmol/m3/s]
		double r_deposition_per_unit_area_;							//!< deposition rate per unit of surface [kmol/m2/s]
		double r_deposition_per_unit_volume_;						//!< deposition rate per unit of surface [kmol/m3/s]
	};
}

#include "PorousMedium.hpp"

#endif /* OpenSMOKE_PorousMedium_H */