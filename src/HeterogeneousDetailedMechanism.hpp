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

namespace CVI
{
	HeterogeneousDetailedMechanism::HeterogeneousDetailedMechanism(	
								OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>& thermodynamicsMap,
								OpenSMOKE::KineticsMap_CHEMKIN<double>& kineticsMap,
								OpenSMOKE::TransportPropertiesMap_CHEMKIN<double>& transportMap,
								OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN<double>& thermodynamicsSurfaceMap,
								OpenSMOKE::KineticsMap_Surface_CHEMKIN<double>&	kineticsSurfaceMap,
								const bool homogeneous_reactions, const bool heterogeneous_reactions) :

	thermodynamicsMap_(thermodynamicsMap),
	kineticsMap_(kineticsMap),
	transportMap_(transportMap),
	thermodynamicsSurfaceMap_(thermodynamicsSurfaceMap),
	kineticsSurfaceMap_(kineticsSurfaceMap),
	homogeneous_reactions_(homogeneous_reactions),
	heterogeneous_reactions_(heterogeneous_reactions)
	{
		Initialize();
	}

	void HeterogeneousDetailedMechanism::Initialize()
	{
		// Constants
		mw_carbon_ = 12.010999679565430;												// [kg/kmol]
		rho_graphite_ = thermodynamicsSurfaceMap_.vector_densities_bulk_species()[0];	// [kg/m3]

		// Number of species
		nc_ = thermodynamicsMap_.NumberOfSpecies();
		nr_ = kineticsMap_.NumberOfReactions();

		// Surface phases
		surf_np_ = thermodynamicsSurfaceMap_.number_of_site_phases(0);
		surf_nc_ = thermodynamicsSurfaceMap_.number_of_site_species();
		surf_nr_ = kineticsSurfaceMap_.NumberOfReactions();

		// Bulk phases
		bulk_np_ = thermodynamicsSurfaceMap_.number_of_bulk_phases(0);
		bulk_nc_ = thermodynamicsSurfaceMap_.number_of_bulk_species();

		// Memory allocation
		Rgas_.resize(nc_);
		Rsurface_.resize(surf_nc_);
		r_heterogeneous_.resize(surf_nr_);
		r_heterogeneous_deposition_per_unit_area_per_single_reaction_.resize(surf_nr_);
		r_heterogeneous_deposition_per_unit_volume_per_single_reaction_.resize(surf_nr_);

		OpenSMOKE::ChangeDimensions(nc_, &Rgas_from_surface_, true);		
		OpenSMOKE::ChangeDimensions(surf_nc_, &Rsurface_from_surface_, true);	
		OpenSMOKE::ChangeDimensions(bulk_nc_, &Rbulk_from_surface_, true);		
		OpenSMOKE::ChangeDimensions(surf_np_, &Rphases_from_surface_, true);
		
		// Summary on the screen
		Summary();
	}

	void HeterogeneousDetailedMechanism::Summary()
	{
		std::cout << std::endl;
		std::cout << "-------------------------------------------------------" << std::endl;
		std::cout << "                      Porous Medium                    " << std::endl;
		std::cout << "-------------------------------------------------------" << std::endl;
		std::cout << " * Graphite density [kg/m3]:         " << rho_graphite_  << std::endl;
		std::cout << "-------------------------------------------------------" << std::endl;
		std::cout << std::endl;
	}

	void HeterogeneousDetailedMechanism::SetTemperature(const double T)
	{
		T_ = T;
	}

	void HeterogeneousDetailedMechanism::SetPressure(const double P_Pa)
	{
		P_Pa_ = P_Pa;
	}

	void HeterogeneousDetailedMechanism::FormationRates(const double Sv, const Eigen::VectorXd& C, const Eigen::VectorXd& Z, const Eigen::VectorXd& a, const Eigen::VectorXd& Gamma)
	{
		// To improve
		OpenSMOKE::OpenSMOKEVectorDouble os_C(C.size());
		OpenSMOKE::OpenSMOKEVectorDouble os_Z(Z.size());
		OpenSMOKE::OpenSMOKEVectorDouble os_a(a.size());
		OpenSMOKE::OpenSMOKEVectorDouble os_Gamma(Gamma.size());
		for (unsigned int i = 0; i < C.size(); i++)		os_C[i + 1] = C(i);
		for (unsigned int i = 0; i < Z.size(); i++)		os_Z[i + 1] = Z(i);
		for (unsigned int i = 0; i < a.size(); i++)		os_a[i + 1] = a(i);
		for (unsigned int i = 0; i < Gamma.size(); i++)	os_Gamma[i + 1] = Gamma(i);

		// Calculates heterogeneous kinetics
		{
			thermodynamicsSurfaceMap_.SetPressure(P_Pa_);
			thermodynamicsSurfaceMap_.SetTemperature(T_);
			kineticsSurfaceMap_.SetPressure(P_Pa_);
			kineticsSurfaceMap_.SetTemperature(T_);
			kineticsSurfaceMap_.ReactionEnthalpiesAndEntropies();
			kineticsSurfaceMap_.ArrheniusKineticConstants();
			kineticsSurfaceMap_.ReactionRates(os_C, os_Z, os_a, os_Gamma);
			kineticsSurfaceMap_.FormationRates(&Rgas_from_surface_, &Rsurface_from_surface_, &Rbulk_from_surface_, &Rphases_from_surface_);
		}

		// Formation rate of homogeneous species [kmol/m3/s]
		for (unsigned int i = 0; i < nc_; i++)
			Rgas_(i) = Sv*Rgas_from_surface_[i+1];

		// Formation rate of homogeneous species [kmol/m2/s]
		for (unsigned int i = 0; i < surf_nc_; i++)
			Rsurface_(i) = Rsurface_from_surface_[i + 1];

		// Total heterogeneous deposition rate [kmol/m2/s]
		r_heterogeneous_deposition_per_unit_area_ = Rbulk_from_surface_[1];

		// Heterogeneous deposition rate [kmol/m3/s]
		r_heterogeneous_deposition_per_unit_volume_ = Sv*r_heterogeneous_deposition_per_unit_area_;

		// Homogeneous rate of disappearance
		r_heterogeneous_deposition_per_unit_volume_ = Sv*r_heterogeneous_deposition_per_unit_area_;

		// Single constributions (TODO)
	}
}
