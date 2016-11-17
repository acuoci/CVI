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

#ifndef OpenSMOKE_HeterogeneousDetailedMechanism_H
#define OpenSMOKE_HeterogeneousDetailedMechanism_H

// OpenSMOKE++ Definitions
#include "OpenSMOKEpp"

// CHEMKIN maps
#include "maps/Maps_CHEMKIN"
#include "maps/ThermodynamicsMap_Surface_CHEMKIN.h"
#include "maps/KineticsMap_Surface_CHEMKIN.h"

namespace CVI
{
	//!  A class to manage properties of heterogeneous reactions
	/*!
	This class provides the tools to manage properties of heterogeneous reactions
	*/

	class HeterogeneousDetailedMechanism
	{
	public:

		/**
		*@brief Default constructor
		*@param thermodynamicsMap			reference to the thermodynamic map
		*@param kineticsMap					reference to the kinetic map
		*@param transportMap				reference to the transport map
		*@param thermodynamicsSurfaceMapXML reference to the surface thermodynamic map
		*@param kineticsSurfaceMapXML		reference to the surface kinetics map
		*@param homogeneous_reactions		homogeneous reactions on/off
		*@param heterogeneous_reactions		heterogeneous reactions on/off
		*/
		HeterogeneousDetailedMechanism(	OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>& thermodynamicsMap,
										OpenSMOKE::KineticsMap_CHEMKIN<double>& kineticsMap,
										OpenSMOKE::TransportPropertiesMap_CHEMKIN<double>& transportMap,
										OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN<double>& thermodynamicsSurfaceMapXML,
										OpenSMOKE::KineticsMap_Surface_CHEMKIN<double>&	kineticsSurfaceMapXML,
										const bool homogeneous_reactions, const bool heterogeneous_reactions);

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
		*@param Sv ratio between available area and gaseous volume [1/m]
		*@param C concentrations of gaseous species [kmol/m3]
		*@param Z fractions of surface species [-]
		*@param a bulk activities [-] (they are always assumed equal to 1)
		*@param Gamma current surface densities [kmol/m2]
		*/
		void FormationRates(const double Sv, const Eigen::VectorXd& C, const Eigen::VectorXd& Z, const Eigen::VectorXd& a, const Eigen::VectorXd& Gamma);
		
		/**
		*@brief Returns the formation rates of gaseous species due to the heterogeneous reactions
		*@return the formation rates [kmol/m3/s]
		*/
		const Eigen::VectorXd& Rgas() const { return Rgas_; }

		/**
		*@brief Returns the formation rates of surface species due to the heterogeneous reactions
		*@return the formation rates [kmol/m2/s]
		*/
		const Eigen::VectorXd& Rsurface() const { return Rsurface_; }

		/**
		*@brief Returns the reaction rates of heterogeneous reactions
		*@return the reaction rates [kmol/m3/s]
		*/
		const Eigen::VectorXd& r() const { return r_heterogeneous_; }

		/**
		*@brief Returns the deposition rate per unit of area
		*@return the deposition rate [kmol/m2/s]
		*/
		double r_deposition_per_unit_area() const { return r_heterogeneous_deposition_per_unit_area_; }

		/**
		*@brief Returns the deposition rate per unit of volume
		*@return the deposition rate [kmol/m3/s]
		*/
		double r_deposition_per_unit_volume() const { return r_heterogeneous_deposition_per_unit_volume_; }

		/**
		*@brief Returns the deposition rate per unit of surface for single reactions [kmol/m2/s]
		*@return the deposition rate per unit of surface for single reactions [kmol/m2/s]
		*/
		const Eigen::VectorXd& r_deposition_per_unit_area_per_single_reaction() const { return r_heterogeneous_deposition_per_unit_area_per_single_reaction_; }

		/**
		*@brief Returns the deposition rate per unit of volume for single reactions [kmol/m3/s]
		*@return the deposition rate per unit of volume for single reactions [kmol/m3/s]
		*/
		const Eigen::VectorXd& r_deposition_per_unit_volume_per_single_reaction() const { return r_heterogeneous_deposition_per_unit_volume_per_single_reaction_; }

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
		*@brief Prints useful data on the screen
		*/
		void Summary();

		void Initialize();

	private:

		OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>&			thermodynamicsMap_;				//!< reference to the thermodynamic map
		OpenSMOKE::KineticsMap_CHEMKIN<double>&					kineticsMap_;					//!< reference to the kinetic map
		OpenSMOKE::TransportPropertiesMap_CHEMKIN<double>&		transportMap_;					//!< reference to the trasport properties map
		OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN<double>&	thermodynamicsSurfaceMap_;		//!< reference to the surface thermodynamic map
		OpenSMOKE::KineticsMap_Surface_CHEMKIN<double>&			kineticsSurfaceMap_;			//!< reference to the surface kinetic map

		unsigned int nc_;	//!< total number of gaseous species
		unsigned int nr_;	//!< total number of homogeneous reactions

		unsigned int surf_np_;	//!< total number of site phases
		unsigned int surf_nc_;	//!< total number of surface species
		unsigned int surf_nr_;	//!< total number of heterogeneous reactions

		unsigned int bulk_np_;	//!< total number of bulk phases
		unsigned int bulk_nc_;	//!< total number of bulk species

		double T_;				//!< current temperature [K]
		double P_Pa_;			//!< current pressure [Pa]
		double rho_graphite_;	//!< density of graphite [kg/m3]
		double mw_carbon_;		//!< molecular weight of carbon matrix [kg/kmol]

		std::vector<std::string> tags_;						//!< tags for single reactions

		double heterogeneous_reaction_rates_multiplier_;					//!< heterogeneous reaction rates multiplier

		bool homogeneous_reactions_;										//!< homogeneous reactions on/off
		bool heterogeneous_reactions_;										//!< heterogeneous reactions on/off
		
		OpenSMOKE::OpenSMOKEVectorDouble Rgas_from_surface_;		//!< current formation rates (from heterogeneous reations) [kmol/m2/s]
		OpenSMOKE::OpenSMOKEVectorDouble Rsurface_from_surface_;	//!< current surface species formation rates [kmol/m2/s]
		OpenSMOKE::OpenSMOKEVectorDouble Rbulk_from_surface_;		//!< current bulk species formation rates [kmol/m2/s]
		OpenSMOKE::OpenSMOKEVectorDouble Rphases_from_surface_;		//!< current surface site phases formation rates [kmol/m2/s]
				
		Eigen::VectorXd Rgas_;		//!< formation rate of gaseous species due to heterogeneous reactions [kmol/m3/s]
		Eigen::VectorXd Rsurface_;	//!< formation rate of surface species due to homogeneous reactions [kmol/m2/s]

		Eigen::VectorXd r_heterogeneous_;													//!< reaction rates of heterogeneous reactions [kmol/m3/s]
		Eigen::VectorXd r_heterogeneous_deposition_per_unit_area_per_single_reaction_;		//!< deposition rate per unit of surface for single reactions [kmol/m2/s]
		Eigen::VectorXd r_heterogeneous_deposition_per_unit_volume_per_single_reaction_;	//!< deposition rate per unit of volume for single reactions [kmol/m3/s]

		double r_heterogeneous_deposition_per_unit_area_;									//!< total deposition rate per unit of surface [kmol/m2/s]
		double r_heterogeneous_deposition_per_unit_volume_;									//!< total deposition rate per unit of volume [kmol/m3/s]

	};
}

#include "HeterogeneousDetailedMechanism.hpp"

#endif /* OpenSMOKE_HeterogeneousDetailedMechanism	_H */
