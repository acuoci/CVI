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

#ifndef OpenSMOKE_HeterogeneousMechanism_H
#define OpenSMOKE_HeterogeneousMechanism_H

// OpenSMOKE++ Definitions
#include "OpenSMOKEpp"

// CHEMKIN maps
#include "maps/Maps_CHEMKIN"

namespace CVI
{
	enum HydrogenInhibitionType { NONE, BECKER };
	enum HeterogeneousMechanismType { ZIEGLER, HUTTINGER, VIGNOLES, ZIEGLER_EXTENDED, HUTTINGER_EXTENDED, VIGNOLES_EXTENDED};

	//!  A class to manage properties of heterogeneous reactions
	/*!
	This class provides the tools to manage properties of heterogeneous reactions
	*/

	class HeterogeneousMechanism
	{
	public:

		/**
		*@brief Default constructor
		*@param thermodynamicsMap		reference to the thermodynamic map
		*@param kineticsMap			reference to the kinetic map
		*@param transportMap			reference to the transport map
		*@param homogeneous_reactions		homogeneous reactions
		*@param heterogeneous_reactions		heterogeneous reactions
		*@param heterogeneous_mechanism_type	heterogeneous mechanism type
		*@param hydrogen_inhibition_type	hydrogen mechanism type
		*/
		HeterogeneousMechanism(	OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap,
						OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap,
						OpenSMOKE::TransportPropertiesMap_CHEMKIN& transportMap,
						const bool homogeneous_reactions, const bool heterogeneous_reactions, 
						const HeterogeneousMechanismType heterogeneous_mechanism_type,
						const HydrogenInhibitionType hydrogen_inhibition_type);

		HeterogeneousMechanism(	OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap,
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
		*@brief Sets the density of graphite
		*@param rho_graphite density [kg/m3]
		*/
		void SetGraphiteDensity(const double rho_graphite);

		/**
		*@brief Returns the density of graphite
		*@return the density of graphite [kg/m3]
		*/
		double rho_graphite() const { return rho_graphite_; }

		/**
		*@brief Returns the molecular weight of the carbon matrix
		*@return the molecular weight of the carbon matrix [kg/kmol]
		*/
		double mw_carbon() const { return mw_carbon_; }

		/**
		*@brief Calculates the formation rates of species due to the heterogeneous reactions
		*@param Sv ratio between avialble area and gaseous volume [1/m]
		*@param C concentrations of gaseous species [kmol/m3]
		*/
		void FormationRates(const double Sv, const Eigen::VectorXd& C);

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
		*@brief Returns the deposition rate per unit of surface for single reactions [kmol/m2/s]
		*@return the deposition rate per unit of surface for single reactions [kmol/m2/s]
		*/
		const Eigen::VectorXd& r_deposition_per_unit_area_per_single_reaction() const { return r_deposition_per_unit_area_per_single_reaction_; }

		/**
		*@brief Returns the deposition rate per unit of volume for single reactions [kmol/m3/s]
		*@return the deposition rate per unit of volume for single reactions [kmol/m3/s]
		*/
		const Eigen::VectorXd& r_deposition_per_unit_volume_per_single_reaction() const { return r_deposition_per_unit_volume_per_single_reaction_; }

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
		*@brief Returns the index for antracene
		*@return the index (0-based) for antracene
		*/
		unsigned int index_C14H10() const { return index_C14H10_; }

		/**
		*@brief Returns the index for antracene
		*@return the index (0-based) for antracene
		*/
		unsigned int index_C10H8() const { return index_C10H8_; }

		/**
		*@brief Returns the index for hydrogen
		*@return the index (0-based) for hydrogen
		*/
		unsigned int index_H2() const { return index_H2_; }

		/**
		*@brief Returns the flag for homogenous reactions in the gaseous phase
		*@return true if the homogeneous reactions are accounted for
		*/
		bool homogeneous_reactions() const { return homogeneous_reactions_; }

		/**
		*@brief Returns the flag for heterogeneous reactions in the gaseous phase
		*@return true if the heterogeneous reactions are accounted for
		*/
		bool heterogeneous_reactions() const { return heterogeneous_reactions_; }

		/**
		*@brief Returns the tags for the single reactions
		*@return the tags for the single reactions
		*/
		const std::vector<std::string>& tags() const { return tags_; }

	private:

		/**
		*@brief Calculates the hydrogen inhibition coefficients
		*/
		void HydrogenInhibitionCoefficients(const Eigen::VectorXd& C);

		/**
		*@brief Prints useful data on the screen
		*/
		void Summary();

		void Initialize();

	private:

		OpenSMOKE::ThermodynamicsMap_CHEMKIN&			thermodynamicsMap_;	//!< reference to the thermodynamic map
		OpenSMOKE::KineticsMap_CHEMKIN&					kineticsMap_;		//!< reference to the kinetic map
		OpenSMOKE::TransportPropertiesMap_CHEMKIN&		transportMap_;		//!< reference to the trasport properties map

		unsigned int ns_;	//!< total number of gaseous species
		unsigned int nr_;	//!< total number of heterogeneous reactions

		int index_CH4_;		//!< index of CH4  (0-based)
		int index_C2H4_;	//!< index of C2H4 (0-based)
		int index_C2H2_;	//!< index of C2H2 (0-based)
		int index_C6H6_;	//!< index of C6H6 (0-based)
		int index_C14H10_;	//!< index of C14H10 (0-based)
		int index_C10H8_;	//!< index of C10H8 (0-based)
		int index_H2_;		//!< index of H2   (0-based)

		double T_;				//!< current temperature [K]
		double P_Pa_;			//!< current pressure [Pa]
		double rho_graphite_;	//!< density of graphite [kg/m3]
		double mw_carbon_;		//!< molecular weight of carbon matrix [kg/kmol]

		HydrogenInhibitionType hydrogen_inhibition_type_;	//!< hydrogene inhibition type
		double I_CH4_;										//!< inhibition coefficient for methane
		double I_C2H4_;										//!< inhibition coefficient for ethylene
		double I_C2H2_;										//!< inhibition coefficient for acetylene
		double I_C6H6_;										//!< inhibition coefficient for benzene
		double I_C14H10_;									//!< inhibition coefficient for antracene
		double I_C10H8_;									//!< inhibition coefficient for naphtalene

		std::vector<std::string> tags_;						//!< tags for single reactions

		double heterogeneous_reaction_rates_multiplier_;			//!< heterogeneous reaction rates multiplier

		HeterogeneousMechanismType heterogeneous_mechanism_type_;		//!< type of heterogeneous mechanism
		bool homogeneous_reactions_;						//!< homogeneous reactions on/off
		bool heterogeneous_reactions_;						//!< heterogeneous reactions on/off
		Eigen::VectorXd Rgas_;							//!< formation rate of gaseous species due to heterogeneous reactions [kmol/m3/s]
		Eigen::VectorXd r_;							//!< reaction rates of heterogeneous reactions [kmol/m3/s]
		double r_deposition_per_unit_area_;					//!< deposition rate per unit of surface [kmol/m2/s]
		double r_deposition_per_unit_volume_;					//!< deposition rate per unit of volume [kmol/m3/s]
		Eigen::VectorXd r_deposition_per_unit_area_per_single_reaction_;	//!< deposition rate per unit of surface for single reactions [kmol/m2/s]
		Eigen::VectorXd r_deposition_per_unit_volume_per_single_reaction_;	//!< deposition rate per unit of volume for single reactions [kmol/m3/s]
	};
}

#include "HeterogeneousMechanism.hpp"

#endif /* OpenSMOKE_HeterogeneousMechanism_H */
