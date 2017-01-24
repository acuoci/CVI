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

// OpenSMOKE++ Solvers
#include "Interface_OpenSMOKEppDae.h"

namespace CVI
{
	Reactor1D::Reactor1D(	OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>& thermodynamicsMap,
							OpenSMOKE::KineticsMap_CHEMKIN<double>& kineticsMap,
							OpenSMOKE::TransportPropertiesMap_CHEMKIN<double>& transportMap,
							OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN<double>& thermodynamicsSurfaceMap,
							OpenSMOKE::KineticsMap_Surface_CHEMKIN<double>&	kineticsSurfaceMap,
							CVI::PorousMedium& porousMedium,
							CVI::HeterogeneousMechanism& heterogeneousMechanism,
							CVI::HeterogeneousDetailedMechanism& heterogeneousDetailedMechanism,
							OpenSMOKE::Grid1D& grid,
							const bool detailed_heterogeneous_kinetics,
							const std::vector<bool>& site_non_conservation,
							const std::string dae_species,
							const boost::filesystem::path output_folder) :

	thermodynamicsMap_(thermodynamicsMap),
	kineticsMap_(kineticsMap),
	transportMap_(transportMap),
	porousMedium_(porousMedium),
	thermodynamicsSurfaceMap_(thermodynamicsSurfaceMap),
	kineticsSurfaceMap_(kineticsSurfaceMap),
	heterogeneousMechanism_(heterogeneousMechanism),
	heterogeneousDetailedMechanism_(heterogeneousDetailedMechanism),
	grid_(grid),
	detailed_heterogeneous_kinetics_(detailed_heterogeneous_kinetics),
	site_non_conservation_(site_non_conservation),
	output_folder_(output_folder)
	{
		n_steps_video_ = 10;
		count_dae_video_ = 1;
		count_ode_video_ = 1;
		n_steps_file_ = 3;
		count_file_ = n_steps_file_;
		ropa_analysis_ = false;

		time_total_ = 48.*3600.;
		dae_time_interval_ = 3600.;
		ode_end_time_ = 1.;

		derivative_type_mass_fractions_ = OpenSMOKE::DERIVATIVE_1ST_CENTERED;
		derivative_type_effective_diffusivity_ = OpenSMOKE::DERIVATIVE_1ST_CENTERED;
		derivative_type_bulk_density_ = OpenSMOKE::DERIVATIVE_1ST_CENTERED;

		if (detailed_heterogeneous_kinetics_ == true)
		{
			// Homogeneous species
			nc_ = thermodynamicsMap_.NumberOfSpecies();
			nr_ = kineticsMap_.NumberOfReactions();

			// Surface phases
			surf_np_ = thermodynamicsSurfaceMap_.number_of_site_phases(0);
			surf_nc_ = thermodynamicsSurfaceMap_.number_of_site_species();
			surf_nr_ = kineticsSurfaceMap_.NumberOfReactions();

			// Bulk phases
			bulk_np_ = thermodynamicsSurfaceMap_.number_of_bulk_phases(0);
			bulk_nc_ = thermodynamicsSurfaceMap_.number_of_bulk_species();

			// Block size
			block_ = nc_ + 1 + surf_np_ + surf_nc_;

			// Graphite density [kg/m3]
			rho_graphite_ = heterogeneousDetailedMechanism_.rho_graphite();

			// Dae species
			dae_formulation_ = true;
			dae_species_index_ = thermodynamicsSurfaceMap_.IndexOfSpecies(dae_species) - (nc_ + 1);
		}
		else
		{
			// Homogeneous species
			nc_ = thermodynamicsMap_.NumberOfSpecies();
			nr_ = kineticsMap_.NumberOfReactions();

			// Surface phases
			surf_np_ = 0;
			surf_nc_ = 0;
			surf_nr_ = heterogeneousMechanism_.r().size();

			// Bulk phases
			bulk_np_ = 0;
			bulk_nc_ = 0;

			// Block size
			block_ = nc_ + 1;

			// Graphite density [kg/m3]
			rho_graphite_ = heterogeneousMechanism_.rho_graphite();

			// Dae Species
			dae_formulation_ = false;
		}
		
		np_ = grid_.Np();
		ne_ = block_*np_;
		band_size_ = 2 * block_ - 1;

		planar_symmetry_ = true;

		output_matlab_folder_ = output_folder_ / "Matlab";
		OpenSMOKE::CreateDirectory(output_matlab_folder_);

		output_diffusion_folder_ = output_folder_ / "DiffusionCoefficients";
		OpenSMOKE::CreateDirectory(output_diffusion_folder_);

		output_heterogeneous_folder_ = output_folder_ / "HeterogeneousReactions";
		OpenSMOKE::CreateDirectory(output_heterogeneous_folder_);

		output_homogeneous_folder_ = output_folder_ / "HomogeneousReactions";
		OpenSMOKE::CreateDirectory(output_homogeneous_folder_);

		output_ropa_folder_ = output_folder_ / "ROPA";
		OpenSMOKE::CreateDirectory(output_ropa_folder_);

		fMonitoring_.open((output_folder_ / "monitor.out").string().c_str(), std::ios::out);
		fMonitoring_.setf(std::ios::scientific);
		PrintLabelMonitoringFile();

		MemoryAllocation();
		SetAlgebraicDifferentialEquations();
	}

	void Reactor1D::MemoryAllocation()
	{
		OpenSMOKE::ChangeDimensions(nc_, &aux_X, true);
		OpenSMOKE::ChangeDimensions(nc_, &aux_Y, true);
		OpenSMOKE::ChangeDimensions(nc_, &aux_C, true);
		OpenSMOKE::ChangeDimensions(nc_, &aux_R, true);
		aux_eigen.resize(nc_);
		
		// Densities [kg/m3]
		rho_gas_.resize(np_);
		rho_bulk_.resize(np_);
		
		// Porosity [-]
		epsilon_.resize(np_);

		// Permeability [m2]
		permeability_.resize(np_);

		// Tortuosities [-]
		eta_bulk_.resize(np_);
		eta_knudsen_.resize(np_);
		eta_viscous_.resize(np_);

		// Surface per unit of volume [1/m]
		Sv_.resize(np_);

		// Porous radius [m]
		rp_.resize(np_);

		// Molecular weight [kg/kmol]
		mw_.resize(np_);

		// Temperature and pressure
		T_.resize(np_);
		P_.resize(np_);

		// Mole fractions [-]
		X_.resize(np_);
		for (unsigned int i = 0; i < np_; i++)
			X_[i].resize(nc_);

		// Mass fractions [-]
		Y_.resize(np_);
		for (unsigned int i = 0; i < np_; i++)
			Y_[i].resize(nc_);

		// Formation rates in gas pahse [kg/m3/s]
		omega_homogeneous_from_homogeneous_.resize(np_);
		for (unsigned int i = 0; i < np_; i++)
			omega_homogeneous_from_homogeneous_[i].resize(nc_);

		// Heterogeneous formation rates [kg/m3/s]
		omega_homogeneous_from_heterogeneous_.resize(np_);
		for (unsigned int i = 0; i < np_; i++)
			omega_homogeneous_from_heterogeneous_[i].resize(nc_);

		// Deposition rate [kg/m3/s]
		omega_deposition_per_unit_volume_.resize(np_);

		// Loss for the homogeneous phase because the heterogeneous reactions [kg/m3/s]
		omega_loss_per_unit_volume_.resize(np_);

		// Deposition rate [kg/m2/s]
		omega_deposition_per_unit_area_.resize(np_);

		// Effective diffusion coefficiennts [m2/s]
		gamma_star_.resize(np_);
		for (unsigned int i = 0; i < np_; i++)
			gamma_star_[i].resize(nc_);

		// Diffusion fluxes [???]
		j_star_.resize(np_ - 1);
		for (unsigned int i = 0; i < np_ - 1; i++)
			j_star_[i].resize(nc_);

		// Correction diffusion flux [???]
		jc_star_.resize(np_ - 1);

		// Spatial derivatives: mole fractions [1/m]
		dX_over_dx_.resize(np_);
		for (unsigned int i = 0; i < np_; i++)
			dX_over_dx_[i].resize(nc_);

		// Spatial derivatives: mass fractions [1/m]
		dY_over_dx_.resize(np_);
		for (unsigned int i = 0; i < np_; i++)
			dY_over_dx_[i].resize(nc_);

		// Second order spatial derivatives: mass fractions [1/m2]
		d2Y_over_dx2_.resize(np_);
		for (unsigned int i = 0; i < np_; i++)
			d2Y_over_dx2_[i].resize(nc_);

		// Time derivatives: mass fractions [1/s]
		dY_over_dt_.resize(np_);
		for (unsigned int i = 0; i < np_; i++)
			dY_over_dt_[i].resize(nc_);

		// Spatial derivatives: mass diffusion coefficients [m2/s/m]
		dgamma_star_over_dx_.resize(np_);
		for (unsigned int i = 0; i < np_; i++)
			dgamma_star_over_dx_[i].resize(nc_);

		// Spatial derivatives: density of gaseous phase
		drho_gas_over_dx_.resize(np_);

		// Time derivatives: porosity [1/s]
		depsilon_over_dt_.resize(np_);

		// Reset to zero
		for (unsigned int i = 0; i < np_-1; i++)
			j_star_[i].setZero();
		for (unsigned int i = 0; i < np_; i++)
			omega_homogeneous_from_homogeneous_[i].setZero();
		for (unsigned int i = 0; i < np_; i++)
			omega_homogeneous_from_heterogeneous_[i].setZero();
		omega_deposition_per_unit_volume_.setZero();
		omega_deposition_per_unit_area_.setZero();
		omega_loss_per_unit_volume_.setZero();

		// Detailed heterogeneous kinetics
		if (detailed_heterogeneous_kinetics_ == true)
		{
			// Surface fractions [-]
			Z_.resize(np_);
			for (unsigned int i = 0; i < np_; i++)
				Z_[i].resize(surf_nc_);

			// Time derivatives: surface species fractions [1/s]
			dZ_over_dt_.resize(np_);
			for (unsigned int i = 0; i < np_; i++)
				dZ_over_dt_[i].resize(surf_nc_);

			Gamma_.resize(np_);
			for (unsigned int i = 0; i < np_; i++)
				Gamma_[i].resize(surf_np_);

			GammaFromEqn_.resize(np_);
			for (unsigned int i = 0; i < np_; i++)
				GammaFromEqn_[i].resize(surf_np_);

			dGamma_over_dt_.resize(np_);
			for (unsigned int i = 0; i < np_; i++)
			{
				dGamma_over_dt_[i].resize(surf_np_);
				dGamma_over_dt_[i].setZero();
			}

			// Heterogeneous formation rates for surface species [kg/m2/s]
			omega_heterogeneous_from_heterogeneous_.resize(np_);
			for (unsigned int i = 0; i < np_; i++)
				omega_heterogeneous_from_heterogeneous_[i].resize(surf_nc_);

			// Auxiliary vectors
			eigen_C_.resize(nc_);
			eigen_Z_.resize(surf_nc_);
			eigen_a_.resize(bulk_nc_);
			eigen_gamma_.resize(surf_np_);

			// Reset to zero
			for (unsigned int i = 0; i < np_; i++)
				omega_heterogeneous_from_heterogeneous_[i].setZero();
		}
	}

	void Reactor1D::SetSurfaceOnTheFlyROPA(OpenSMOKE::SurfaceOnTheFlyROPA* ropa)
	{
		ropa_ = ropa;
		ropa_analysis_ = true;
	}

	void Reactor1D::SetTimeTotal(const double time_total)
	{
		time_total_ = time_total;
	}

	void Reactor1D::SetDaeTimeInterval(const double time_interval)
	{
		dae_time_interval_ = time_interval;
	}

	void Reactor1D::SetOdeEndTime(const double time_interval)
	{
		ode_end_time_ = time_interval;
	}

	void Reactor1D::SetPlanarSymmetry(const bool flag)
	{
		planar_symmetry_ = flag;
	}

	void Reactor1D::SetDerivativeMassFractions(const OpenSMOKE::derivative_type value)
	{
		derivative_type_mass_fractions_ = value;
	}

	void Reactor1D::SetDerivativeEffectiveDiffusivity(const OpenSMOKE::derivative_type value)
	{
		derivative_type_effective_diffusivity_ = value;
	}

	void Reactor1D::SetDerivativeBulkDensity(const OpenSMOKE::derivative_type value)
	{
		derivative_type_bulk_density_ = value;
	}

	void Reactor1D::SetSiteNonConservation(std::vector<bool>& site_non_conservation)
	{
		site_non_conservation_ = site_non_conservation;
	}

	void Reactor1D::SetAlgebraicDifferentialEquations()
	{
		id_equations_.resize(ne_);

		unsigned int count = 0;

		// Internal plane
		{
			for (unsigned int j = 0; j < nc_; j++)
				id_equations_[count++] = false;
			id_equations_[count++] = true;
			for (unsigned int j = 0; j < surf_np_; j++)
				id_equations_[count++] = true;

			if (dae_formulation_ == false)
			{
				for (unsigned int j = 0; j < surf_nc_; j++)
					id_equations_[count++] = true;
			}
			else
			{
				for (unsigned int j = 0; j < dae_species_index_; j++)
					id_equations_[count++] = true;
				id_equations_[count++] = false;
				for (unsigned int j = dae_species_index_ + 1; j < surf_nc_; j++)
					id_equations_[count++] = true;
			}
		}

		// Internal points
		for (int i = 1; i < grid_.Ni(); i++)
		{
			for (unsigned int j = 0; j < nc_; j++)
				id_equations_[count++] = true;
			id_equations_[count++] = true;
			for (unsigned int j = 0; j < surf_np_; j++)
				id_equations_[count++] = true;

			if (dae_formulation_ == false)
			{
				for (unsigned int j = 0; j < surf_nc_; j++)
					id_equations_[count++] = true;
			}
			else
			{
				for (unsigned int j = 0; j < dae_species_index_; j++)
					id_equations_[count++] = true;
				id_equations_[count++] = false;
				for (unsigned int j = dae_species_index_ + 1; j < surf_nc_; j++)
					id_equations_[count++] = true;
			}
		}

		// Gas side
		{
			for (unsigned int j = 0; j < nc_; j++)
				id_equations_[count++] = false;
			id_equations_[count++] = true;
			for (unsigned int j = 0; j < surf_np_; j++)
				id_equations_[count++] = true;

			if (dae_formulation_ == false)
			{
				for (unsigned int j = 0; j < surf_nc_; j++)
					id_equations_[count++] = true;
			}
			else
			{
				for (unsigned int j = 0; j < dae_species_index_; j++)
					id_equations_[count++] = true;
				id_equations_[count++] = false;
				for (unsigned int j = dae_species_index_ + 1; j < surf_nc_; j++)
					id_equations_[count++] = true;
			}
		}
	}

	void Reactor1D::SetGasSide(const double T_gas, const double P_gas, const Eigen::VectorXd& omega_gas)
	{
		Y_gas_side_.resize(nc_);
		for (unsigned int j = 0; j < nc_; j++)
			Y_gas_side_(j) = omega_gas(j);

		P_gas_side_ = P_gas;
		T_gas_side_ = T_gas;

		// Internal side
		if (planar_symmetry_ == false)
		{
			for (unsigned int j = 0; j < nc_; j++)
				Y_[0](j) = omega_gas(j);

			P_(0) = P_gas;
			T_(0) = T_gas;
		}

		// External side
		{
			for (unsigned int j = 0; j < nc_; j++)
				Y_[grid_.Ni()](j) = omega_gas(j);

			P_(grid_.Ni()) = P_gas;
			T_(grid_.Ni()) = T_gas;
		}
	}

	void Reactor1D::SetInitialConditions(const double T_gas, const double P_gas, const Eigen::VectorXd& omega_gas, const Eigen::VectorXd& Gamma0, const Eigen::VectorXd& Z0)
	{
		for (unsigned int i = 0; i < np_; i++)
			for (unsigned int j = 0; j < nc_; j++)
				Y_[i](j) = omega_gas(j);

		for (unsigned int i = 0; i < np_; i++)
			for (unsigned int j = 0; j < surf_nc_; j++)
				Z_[i](j) = Z0(j);

		for (unsigned int i = 0; i < np_; i++)
			for (unsigned int j = 0; j < surf_np_; j++)
			{
				Gamma_[i](j) = Gamma0(j);
				GammaFromEqn_[i](j) = Gamma0(j);
			}

		P_.setConstant(P_gas);
		T_.setConstant(T_gas);
		epsilon_.setConstant(porousMedium_.porosity());
	}

	void Reactor1D::Properties()
	{
		for (unsigned int i = 0; i < np_; i++)
		{
			// Thermodynamics
			{
				thermodynamicsMap_.SetPressure(P_(i));
				thermodynamicsMap_.SetTemperature(T_(i));

				for (unsigned int j = 0; j < nc_; j++)
					aux_Y[j + 1] = Y_[i](j);
				thermodynamicsMap_.MoleFractions_From_MassFractions(aux_X, mw_(i), aux_Y);
				for (unsigned int j = 0; j < nc_; j++)
					X_[i](j) = aux_X[j + 1];

				// Concentrations [kmol/m3]
				const double cTot = P_(i) / PhysicalConstants::R_J_kmol / T_(i); // [kmol/m3]
				Product(cTot, aux_X, &aux_C);

				// Mixture density
				rho_gas_(i) = cTot*mw_(i);	// [kg/m3]
			}

			// Transport properties (porous medium)
			{
				porousMedium_.SetPorosity(epsilon_(i));
				porousMedium_.SetTemperature(T_(i));
				porousMedium_.SetPressure(P_(i));

				// Properties
				Sv_(i) = porousMedium_.Sv();
				rp_(i) = porousMedium_.rp();
				permeability_(i) = porousMedium_.permeability();
				rho_bulk_(i) = porousMedium_.density_bulk(rho_graphite_);
				eta_bulk_(i) = porousMedium_.eta_bulk();
				eta_knudsen_(i) = porousMedium_.eta_knudsen();
				eta_viscous_(i) = porousMedium_.eta_viscous();

				// Mixture diffusion coefficients
				{
					for (unsigned int j = 0; j < nc_; j++)
						aux_eigen(j) = aux_X[j + 1];
					porousMedium_.EffectiveDiffusionCoefficients(aux_eigen);
					for (unsigned int j = 0; j < nc_; j++)
						gamma_star_[i](j) = porousMedium_.gamma_effective()(j);
				}
			}

			// Kinetics
			if (detailed_heterogeneous_kinetics_ == false)
			{
				heterogeneousMechanism_.SetTemperature(T_(i));
				heterogeneousMechanism_.SetPressure(P_(i));

				// Homogeneous phase
				if (heterogeneousMechanism_.homogeneous_reactions() == true)
				{
					kineticsMap_.SetTemperature(T_(i));
					kineticsMap_.SetPressure(P_(i));
					kineticsMap_.ReactionRates(aux_C);
					kineticsMap_.FormationRates(&aux_R);
					ElementByElementProduct(aux_R, thermodynamicsMap_.MW(), &aux_R); // [kg/m3/s]
					aux_R.CopyTo(omega_homogeneous_from_homogeneous_[i].data());
				}

				// Heterogeneous phase (global)
				if (heterogeneousMechanism_.heterogeneous_reactions() == true)
				{
					heterogeneousMechanism_.SetTemperature(T_(i));
					heterogeneousMechanism_.SetPressure(P_(i));

					for (unsigned int j = 0; j < nc_; j++)
						aux_eigen(j) = aux_C[j + 1];
					heterogeneousMechanism_.FormationRates(Sv_(i), aux_eigen);
					for (unsigned int j = 0; j < nc_; j++)
						omega_homogeneous_from_heterogeneous_[i](j) = heterogeneousMechanism_.Rgas()(j)*thermodynamicsMap_.MW()[j + 1];						// [kg/m3/s]

					omega_deposition_per_unit_area_(i) = heterogeneousMechanism_.r_deposition_per_unit_area()*heterogeneousMechanism_.mw_carbon();			// [kg/m2/s]
					omega_deposition_per_unit_volume_(i) = heterogeneousMechanism_.r_deposition_per_unit_volume()*heterogeneousMechanism_.mw_carbon();		// [kg/m3/s]

					omega_loss_per_unit_volume_(i) = 0.;
					for (unsigned int j = 0; j < nc_; j++)
						omega_loss_per_unit_volume_(i) += heterogeneousMechanism_.Rgas()(j)*thermodynamicsMap_.MW()[j + 1];					// [kg/m3/s]
				}
			}
			else
			{
				heterogeneousDetailedMechanism_.SetTemperature(T_(i));
				heterogeneousDetailedMechanism_.SetPressure(P_(i));

				// Homogeneous phase
				if (heterogeneousDetailedMechanism_.homogeneous_reactions() == true)
				{
					kineticsMap_.SetTemperature(T_(i));
					kineticsMap_.SetPressure(P_(i));
					kineticsMap_.ReactionRates(aux_C);
					kineticsMap_.FormationRates(&aux_R);
					ElementByElementProduct(aux_R, thermodynamicsMap_.MW(), &aux_R); // [kg/m3/s]
					aux_R.CopyTo(omega_homogeneous_from_homogeneous_[i].data());
				}

				// Heterogeneous phase (detailed)
				if (heterogeneousDetailedMechanism_.heterogeneous_reactions() == true )
				{
					heterogeneousDetailedMechanism_.SetTemperature(T_(i));
					heterogeneousDetailedMechanism_.SetPressure(P_(i));

					for (unsigned int j = 0; j < nc_; j++)
						eigen_C_(j) = aux_C[j + 1];

					for (unsigned int j = 0; j < surf_nc_; j++)
						eigen_Z_(j) = Z_[i](j);

					for (unsigned int j = 0; j < bulk_nc_; j++)
						eigen_a_(j) = 1.;

					for (unsigned int j = 0; j < surf_np_; j++)
						eigen_gamma_(j) = Gamma_[i](j);

					heterogeneousDetailedMechanism_.FormationRates(Sv_(i), eigen_C_, eigen_Z_, eigen_a_, eigen_gamma_);

					for (unsigned int j = 0; j < nc_; j++)
						omega_homogeneous_from_heterogeneous_[i](j) = heterogeneousDetailedMechanism_.Rgas()(j)*thermodynamicsMap_.MW()[j + 1];				// [kg/m3/s]

					for (unsigned int j = 0; j < surf_nc_; j++)
						omega_heterogeneous_from_heterogeneous_[i](j) = heterogeneousDetailedMechanism_.Rsurface()(j);	// [kmol/m2/s]

					omega_deposition_per_unit_area_(i) = heterogeneousDetailedMechanism_.r_deposition_per_unit_area()*heterogeneousDetailedMechanism_.mw_carbon();			// [kg/m2/s]
					omega_deposition_per_unit_volume_(i) = heterogeneousDetailedMechanism_.r_deposition_per_unit_volume()*heterogeneousDetailedMechanism_.mw_carbon();		// [kg/m3/s]

					omega_loss_per_unit_volume_(i) = 0.;
					for (unsigned int j = 0; j < nc_; j++)
						omega_loss_per_unit_volume_(i) += heterogeneousDetailedMechanism_.Rgas()(j)*thermodynamicsSurfaceMap_.MW()[j + 1];					// [kg/m3/s]
				}
			}
		}
	}

	void Reactor1D::DiffusionFluxes()
	{
		// Effective diffusion fluxes
		{
			Eigen::VectorXd dummy(np_); dummy.setZero();
			grid_.Derivative(OpenSMOKE::DERIVATIVE_1ST_BACKWARD, dummy, Y_, &dY_over_dx_);
			
			for (int i = 0; i < grid_.Ni(); i++)
			{
				j_star_[i].setZero();
				for (unsigned int j = 0; j < nc_; j++)
					j_star_[i](j) = -0.50 * (rho_gas_(i)*gamma_star_[i](j) + rho_gas_(i + 1)*gamma_star_[i + 1](j)) *dY_over_dx_[i + 1](j);
			}
		}

		// Correction diffusion velocity
		{
			jc_star_.setZero();
			for (int i = 0; i < grid_.Ni(); i++)
			{
				for (unsigned int j = 0; j < nc_; j++)
					jc_star_(i) -= j_star_[i](j);
			}
		}

		// Total diffusion velocity
		{
			for (int i = 0; i < grid_.Ni(); i++)
			{
				for (unsigned int j = 0; j < nc_; j++)
					j_star_[i](j) += 0.50*(Y_[i](j) + Y_[i + 1](j))*jc_star_(i);
			}
		}
	}

	void Reactor1D::SpatialDerivatives()
	{
		Eigen::VectorXd dummy(np_);

		// First-order derivatives
		grid_.Derivative(derivative_type_mass_fractions_, dummy, Y_, &dY_over_dx_);
		grid_.Derivative(derivative_type_effective_diffusivity_, dummy, gamma_star_, &dgamma_star_over_dx_);
		grid_.Derivative(derivative_type_bulk_density_, dummy, rho_gas_, &drho_gas_over_dx_);

		// Second order derivatives
		grid_.SecondDerivative(Y_, &d2Y_over_dx2_);
	}

	void Reactor1D::SubEquations_MassFractions()
	{
		// Internal side
		if (planar_symmetry_ == true)
		{
			for (unsigned int j = 0; j < nc_; j++)
				dY_over_dt_[0](j) = Y_[0](j) - Y_[1](j);
		}
		else
		{
			for (unsigned int j = 0; j < nc_; j++)
				dY_over_dt_[0](j) = Y_[0](j) - Y_gas_side_(j);
		}

		// Internal points
		for (int i = 1; i < grid_.Ni(); i++)
		{
			for (unsigned int j = 0; j < nc_; j++)
			{
				// Explicit derivatives
				if (planar_symmetry_ == true)
				{
					dY_over_dt_[i](j) = 	gamma_star_[i](j)*rho_gas_(i)*d2Y_over_dx2_[i](j) +
								(gamma_star_[i](j)*drho_gas_over_dx_(i) + dgamma_star_over_dx_[i](j)*rho_gas_(i))*dY_over_dx_[i](j);
				}
				else
				{
					dY_over_dt_[i](j) = 	gamma_star_[i](j)*rho_gas_(i)*d2Y_over_dx2_[i](j) +
								(gamma_star_[i](j)*drho_gas_over_dx_(i) + dgamma_star_over_dx_[i](j)*rho_gas_(i))*dY_over_dx_[i](j) +
								gamma_star_[i](j)*rho_gas_(i)*(Y_[i+1](j)-Y_[i-1](j))/(grid_.x()[i+1]-grid_.x()[i-1]) / grid_.x()[i];
				}

				dY_over_dt_[i](j) += epsilon_(i)*omega_homogeneous_from_homogeneous_[i](j) + omega_homogeneous_from_heterogeneous_[i](j) - Y_[i](j)*omega_loss_per_unit_volume_(i);

				dY_over_dt_[i](j) /= (rho_gas_(i)*epsilon_(i));
			}
		}

		// External side
		for (unsigned int j = 0; j < nc_; j++)
			dY_over_dt_[grid_.Ni()](j) = Y_[grid_.Ni()](j) - Y_gas_side_(j);
	}

	void Reactor1D::SubEquations_Porosity()
	{
		// Internal points
		for (unsigned int i = 0; i < np_; i++)
		{
			const double coefficient = 1000.;
			const double smoothing_coefficient = 0.50*(std::tanh(coefficient*(epsilon_(i)-porousMedium_.epsilon_threshold()))+1.);
			depsilon_over_dt_(i) = -omega_deposition_per_unit_volume_(i) / rho_graphite_*smoothing_coefficient;
		}
	}

	void Reactor1D::SubEquations_SurfaceSpeciesFractions()
	{
		for (int i = 0; i < grid_.Np(); i++)
		{
			for (unsigned int j = 0; j < surf_np_; j++)
			{
				if (site_non_conservation_[j] == true)
					dGamma_over_dt_[i](j) = heterogeneousDetailedMechanism_.Rphases()(j);
				else
					dGamma_over_dt_[i](j) = 0.;
			}

			for (unsigned int j = 0; j < surf_nc_; j++)
			{
				const unsigned int index_phase = thermodynamicsSurfaceMap_.vector_site_phases_belonging()[j];
				dZ_over_dt_[i](j) = (thermodynamicsSurfaceMap_.vector_occupancies_site_species()[j] * omega_heterogeneous_from_heterogeneous_[i](j)
					- Z_[i](j)*dGamma_over_dt_[i](index_phase))
					/ Gamma_[i](index_phase);
			}

			if (dae_formulation_ == true)
				dZ_over_dt_[i](dae_species_index_) = 1. - Z_[i].sum();
		}
	}

	int Reactor1D::OdeEquations(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y, OpenSMOKE::OpenSMOKEVectorDouble& dy)
	{
		OpenSMOKE::OpenSMOKEVectorDouble Gamma(surf_np_);
		OpenSMOKE::OpenSMOKEVectorDouble Z(surf_nc_);
		OpenSMOKE::OpenSMOKEVectorDouble dGamma_over_dt(surf_np_);
		OpenSMOKE::OpenSMOKEVectorDouble dZ_over_dt(surf_nc_);
		OpenSMOKE::OpenSMOKEVectorDouble RfromSurface(nc_);
		OpenSMOKE::OpenSMOKEVectorDouble Rsurface(surf_nc_);
		OpenSMOKE::OpenSMOKEVectorDouble Rbulk(bulk_nc_);
		OpenSMOKE::OpenSMOKEVectorDouble RsurfacePhases(surf_np_);

		unsigned int k = 1;
		for (unsigned int j = 0; j < surf_np_; j++)
			Gamma[j + 1] = y[k++];
		for (unsigned int j = 0; j < surf_nc_; j++)
			Z[j + 1] = y[k++];

		// Molar fractions
		for (unsigned int j = 0; j < nc_; j++)
			aux_Y[j + 1] = Y_[i_current](j);
		thermodynamicsMap_.MoleFractions_From_MassFractions(aux_X, mw_(i_current), aux_Y);

		// Calculates the volume and the concentrations of species
		const double cTot = P_(i_current) / (PhysicalConstants::R_J_kmol * T_(i_current));
		Product(cTot, aux_X, &aux_C);

		// Calculates heterogeneous kinetics
		{
			thermodynamicsSurfaceMap_.SetTemperature(T_(i_current));
			thermodynamicsSurfaceMap_.SetPressure(P_(i_current));
			kineticsSurfaceMap_.SetTemperature(T_(i_current));
			kineticsSurfaceMap_.SetPressure(P_(i_current));

			kineticsSurfaceMap_.KineticConstants();

			OpenSMOKE::OpenSMOKEVectorDouble aux_a(bulk_nc_);
			for (unsigned int j = 0; j < bulk_nc_; j++)
				aux_a[j + 1] = 1.;

			kineticsSurfaceMap_.ReactionRates(aux_C, Z, aux_a, Gamma);
			kineticsSurfaceMap_.FormationRates(&RfromSurface, &Rsurface, &Rbulk, &RsurfacePhases);
		}

		// Equations
		{
			unsigned int k = 1;

			// Phases
			for (unsigned int j = 0; j < surf_np_; ++j)
			{
				if (site_non_conservation_[j] == true)
					dGamma_over_dt[j + 1] = RsurfacePhases[j + 1];
				else
					dGamma_over_dt[j + 1] = 0.;
			}

			// Heterogeneous species
			for (unsigned int j = 0; j < surf_nc_; ++j)
			{
				const unsigned int index_phase = thermodynamicsSurfaceMap_.vector_site_phases_belonging()[j];
				dZ_over_dt[j + 1] = (thermodynamicsSurfaceMap_.vector_occupancies_site_species()[j] * Rsurface[j + 1] -
					Z[j + 1] * dGamma_over_dt[index_phase + 1]) / Gamma[index_phase + 1];
			}
		}

		// Recover unknowns
		{
			unsigned int k = 1;
			for (unsigned int j = 0; j < surf_np_; j++)
				dy[k++] = dGamma_over_dt[j + 1];
			for (unsigned int j = 0; j < surf_nc_; j++)
				dy[k++] = dZ_over_dt[j + 1];
		}

		return 0;
	}

	int Reactor1D::OdePrint(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y)
	{
		// Video output
		if (count_ode_video_%n_steps_video_ == 1)
		{
			if (count_dae_video_ % (n_steps_video_ * 1000) == 1)
			{
				std::cout << std::endl;
				std::cout << std::setw(18) << std::left << "Time[s]";
				std::cout << std::endl;
			}

			std::cout << std::scientific << std::setw(18) << std::setprecision(6) << std::left << t;
			std::cout << std::endl;
		}

		count_ode_video_++;

		return 0;
	}

	void Reactor1D::Recover_Unknowns(const double* y)
	{
		unsigned int count = 0;
		for (unsigned int i = 0; i < np_; i++)
		{
			// Species
			for (unsigned int j = 0; j < nc_; j++)
				Y_[i](j) = y[count++];

			// Porosity
			epsilon_(i) = y[count++];

			// Surface densities
			for (unsigned int j = 0; j < surf_np_; j++)
				GammaFromEqn_[i](j) = y[count++];

			// Surface species
			for (unsigned int j = 0; j < surf_nc_; j++)
				Z_[i](j) = y[count++];

			// Non conservation of sites
			for (unsigned int j = 0; j < surf_np_; j++)
				if (site_non_conservation_[j] == true)
					Gamma_[i](j) = GammaFromEqn_[i](j);
		}
	}


	void Reactor1D::Recover_Residuals(double* dy)
	{
		unsigned int count = 0;
		for (unsigned int i = 0; i < np_; i++)
		{
			// Species
			for (unsigned int j = 0; j < nc_; j++)
				dy[count++] = dY_over_dt_[i](j);
			
			// Porosity
			dy[count++] = depsilon_over_dt_(i);

			// Surface densities
			for (unsigned int j = 0; j < surf_np_; j++)
				dy[count++] = dGamma_over_dt_[i](j);

			// Surface species
			for (unsigned int j = 0; j < surf_nc_; j++)
				dy[count++] = dZ_over_dt_[i](j);
		}
	}

	void Reactor1D::AlgebraicDifferentialVector(double* v)
	{
		int count = 0;
		for (unsigned int i = 0; i < id_equations_.size(); i++)
		{
			if (id_equations_[i] == true)  v[count++] = 1.;
			if (id_equations_[i] == false) v[count++] = 0.;
		}
	}

	void Reactor1D::UnknownsVector(double* v)
	{
		unsigned int count = 0;
		for (unsigned int i = 0; i < np_; i++)
		{
			for (unsigned int j = 0; j < nc_; j++)
				v[count++] = Y_[i](j);

			v[count++] = epsilon_[i];

			for (unsigned int j = 0; j < surf_np_; j++)
				v[count++] = GammaFromEqn_[i](j);

			for (unsigned int j = 0; j < surf_nc_; j++)
				v[count++] = Z_[i](j);
		}
	}

	void Reactor1D::CorrectedUnknownsVector(double* v)
	{
		unsigned int count = 0;
		for (unsigned int i = 0; i < np_; i++)
		{
			for (unsigned int j = 0; j < nc_; j++)
				Y_[i](j) = v[count++];

			epsilon_(i) = v[count++];

			for (unsigned int j = 0; j < surf_np_; j++)
				GammaFromEqn_[i](j) = v[count++];

			for (unsigned int j = 0; j < surf_nc_; j++)
				Z_[i](j) = v[count++];

			// Non conservation of sites
			for (unsigned int j = 0; j < surf_np_; j++)
				if (site_non_conservation_[j] == true)
					Gamma_[i](j) = GammaFromEqn_[i](j);
		}
	}

	void Reactor1D::MinimumUnknownsVector(double* v)
	{
		const double zero = 0.;

		unsigned int count = 0;
		for (unsigned int i = 0; i < np_; i++)
		{
			for (unsigned int j = 0; j < nc_; j++)
				v[count++] = zero;

			v[count++] = zero;

			for (unsigned int j = 0; j < surf_np_; j++)
				v[count++] = zero;

			for (unsigned int j = 0; j < surf_nc_; j++)
				v[count++] = zero;
		}
	}

	void Reactor1D::MaximumUnknownsVector(double* v)
	{
		const double one = 1.;

		unsigned int count = 0;
		for (unsigned int i = 0; i < np_; i++)
		{
			for (unsigned int j = 0; j < nc_; j++)
				v[count++] = one;

			v[count++] = one;

			for (unsigned int j = 0; j < surf_np_; j++)
				v[count++] = one;

			for (unsigned int j = 0; j < surf_nc_; j++)
				v[count++] = one;
		}
	}

	void Reactor1D::Equations(const double t, const double* y, double* dy)
	{
		// Recover unknowns
		Recover_Unknowns(y);

		// Properties
		Properties();

		// Spatial derivatives
		SpatialDerivatives();

		// Fluxes
		// TODO: they are not used because the closure is done by directly imposing the summ of mass fractions equal to 1
		// DiffusionFluxes();

		// Equations
		SubEquations_MassFractions();
		SubEquations_Porosity();
		SubEquations_SurfaceSpeciesFractions();

		// Recover residuals
		Recover_Residuals(dy);
	}

	int Reactor1D::SolveFromScratch(DaeSMOKE::DaeSolver_Parameters& dae_parameters, OdeSMOKE::OdeSolver_Parameters& ode_parameters)
	{
		std::ofstream fMonitoring;
		fMonitoring.open((output_folder_ / "log.out").string().c_str(), std::ios::out);
		fMonitoring.setf(std::ios::scientific);

		if (detailed_heterogeneous_kinetics_ == true)
		{
			// Solve independent ODE systems
			for (unsigned int i = 0; i < np_; i++)
			{
				i_current = i;
				unsigned int NE_ODE = surf_np_ + surf_nc_;

				std::cout << std::endl;
				std::cout << "Calculating initial surface species composition for point " << i + 1 << "/" << np_ << std::endl;
				std::cout << "--------------------------------------------------------------------------" << std::endl;

				// Initial conditions
				OpenSMOKE::OpenSMOKEVectorDouble yOde0(NE_ODE);
				OpenSMOKE::OpenSMOKEVectorDouble yOdef(NE_ODE);
				unsigned int k = 1;
				for (unsigned int j = 0; j < surf_np_; j++)
					yOde0[k++] = Gamma_[i](j);
				for (unsigned int j = 0; j < surf_nc_; j++)
					yOde0[k++] = Z_[i](j);

				// Print intial conditions
				{
					OpenSMOKE::OpenSMOKEVectorDouble dy0(yOde0.Size());
					OdeEquations(0., yOde0, dy0);
					OdePrint(0., yOde0);
				}

				{
					// Min and max values
					Eigen::VectorXd yMin(NE_ODE); for (unsigned int j = 0; j < NE_ODE; j++) yMin(j) = 0.;
					Eigen::VectorXd yMax(NE_ODE); for (unsigned int j = 0; j < NE_ODE; j++) yMax(j) = 1.;

					// Initial conditions
					Eigen::VectorXd y0_eigen(yOde0.Size());
					yOde0.CopyTo(y0_eigen.data());

					// Final f
					Eigen::VectorXd yf_eigen(y0_eigen.size());

					// Create the solver
					typedef OdeSMOKE::KernelDense<OpenSMOKE::ODESystem_OpenSMOKE_Reactor1D> denseOde;
					typedef OdeSMOKE::MethodGear<denseOde> methodGear;
					OdeSMOKE::MultiValueSolver<methodGear> ode_solver;
					ode_solver.SetReactor1D(this);

					// Set initial conditions
					ode_solver.SetInitialConditions(0., y0_eigen);

					// Set linear algebra options
					//ode_solver.SetLinearAlgebraSolver(ode_parameters.);
					//ode_solver.SetFullPivoting(ode_parameters.full_pivoting());

					// Set relative and absolute tolerances
					ode_solver.SetAbsoluteTolerances(ode_parameters.absolute_tolerance());
					ode_solver.SetRelativeTolerances(ode_parameters.relative_tolerance());

					// Set minimum and maximum values
					ode_solver.SetMinimumValues(yMin);
					ode_solver.SetMaximumValues(yMax);

					// Set maximum number of steps
					if (ode_parameters.maximum_number_of_steps() > 0)
						ode_solver.SetMaximumNumberOfSteps(ode_parameters.maximum_number_of_steps());

					// Set maximum integration order
					if (ode_parameters.maximum_order() > 0)
						ode_solver.SetMaximumOrder(ode_parameters.maximum_order());

					// Set maximum step size allowed
					if (ode_parameters.maximum_step() > 0)
						ode_solver.SetMaximumStepSize(ode_parameters.maximum_step());

					// Set minimum step size allowed
					if (ode_parameters.minimum_step() > 0)
						ode_solver.SetMinimumStepSize(ode_parameters.minimum_step());

					// Set initial step size
					if (ode_parameters.initial_step() > 0)
						ode_solver.SetFirstStepSize(ode_parameters.initial_step());

					if (ode_parameters.minimum_yp() > 0)
						ode_solver.SetStopConditionMaximumYPrimeNorm1(ode_parameters.minimum_yp());

					// Solve the system
					double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
					OdeSMOKE::OdeStatus status = ode_solver.Solve(ode_end_time_);
					double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();

					// Check the solution
					if (status > 0)
					{
						ode_solver.Solution(yf_eigen);
						yOdef.CopyFrom(yf_eigen.data());

						unsigned int k = 0;
						for (unsigned int j = 0; j < surf_np_; j++)
							Gamma_[i](j) = yf_eigen(k++);
						for (unsigned int j = 0; j < surf_nc_; j++)
							Z_[i](j) = yf_eigen(k++);

						// Normalize
						const double sum_local_Z = Z_[i].sum();
						for (unsigned int j = 0; j < surf_nc_; j++)
							Z_[i](j) /= sum_local_Z;
					}
				}
			}
		}

		std::cout << "--------------------------------------------------------------------------" << std::endl;
		std::cout << " Solving the whole DAE system (fully coupled)                             " << std::endl;
		std::cout << "--------------------------------------------------------------------------" << std::endl;

		// Loop
		unsigned int number_intervals = static_cast<unsigned int>(time_total_ / dae_time_interval_);
		for (unsigned int k = 1; k <= number_intervals; k++)
		{
			const double t0 = (k - 1)*dae_time_interval_;
			const double tf = t0 + dae_time_interval_;

			// Solve
			int flag = Solve(dae_parameters, t0, tf);
			if (flag < 0)
				return flag;

			// Write current solution
			std::stringstream current_index; current_index << k;
			std::stringstream hours; hours << tf / 3600.;
			std::string solution_file = "Solution." + hours.str() + ".out";
			std::string diffusion_coefficients_file = "DiffusionCoefficients." + hours.str() + ".out";
			std::string heterogeneous_rates_file = "HeterogeneousRates." + hours.str() + ".out";
			std::string homogeneous_rates_file = "HomogeneousRates." + hours.str() + ".out";
			std::string ropa_file = "ROPA." + hours.str() + ".out";
			PrintSolution(tf, (output_folder_ / solution_file).string().c_str());
			PrintDiffusionCoefficients(tf, (output_diffusion_folder_ / diffusion_coefficients_file).string().c_str());
			PrintHeterogeneousRates(tf, (output_heterogeneous_folder_ / heterogeneous_rates_file).string().c_str());
			PrintHomogeneousRates(tf, (output_homogeneous_folder_ / homogeneous_rates_file).string().c_str());

			if (ropa_analysis_ == true)
				PrintROPA(tf, (output_ropa_folder_ / ropa_file).string().c_str());
		}

		return true;
	}

	int Reactor1D::Solve(DaeSMOKE::DaeSolver_Parameters& dae_parameters, const double t0, const double tf)
	{
		int flag = DaeSMOKE::Solve_Band_OpenSMOKEppDae<Reactor1D, OpenSMOKE_Reactor1D_DaeSystem>(this, dae_parameters, t0, tf);
		return flag;
	}

	void Reactor1D::PrintSolution(const double t, const std::string name_file)
	{
		std::ofstream fOutput(name_file.c_str(), std::ios::out);

		{
			unsigned int count = 1;
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "time[s]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "x[mm]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "y[mm]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "T[K]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "P[Pa]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "eps[-]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "rhoBulk[kg/m3]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "Sv[1/m]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "rp[micron]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "K[m2]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "TauBulk[-]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "TauKnudsen[-]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "TauViscous[-]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "rhoG[kg/m3]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "mwG[kg/kmol]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "sumOmega[-]", count);


			for (unsigned int j = 0; j < thermodynamicsMap_.elements().size(); j++)
				OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, thermodynamicsMap_.elements()[j] + "_x", count);
			for (unsigned int j = 0; j < nc_; j++)
				OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, thermodynamicsMap_.NamesOfSpecies()[j] + "_x", count);
			for (unsigned int j = 0; j < nc_; j++)
				OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, thermodynamicsMap_.NamesOfSpecies()[j] + "_w", count);

			fOutput << std::endl;
		}

		for (unsigned int i = 0; i < np_; i++)
		{
			OpenSMOKE::OpenSMOKEVectorDouble yy(nc_);
			OpenSMOKE::OpenSMOKEVectorDouble xx(nc_);

			for (unsigned int j = 0; j < nc_; j++)
				yy[j + 1] = Y_[i](j);

			thermodynamicsMap_.SetPressure(P_(i));
			thermodynamicsMap_.SetTemperature(T_(i));
			thermodynamicsMap_.MoleFractions_From_MassFractions(xx, mw_(i), yy);

			fOutput << std::setprecision(9) << std::setw(20) << t;
			fOutput << std::setprecision(9) << std::setw(20) << grid_.x()[i] * 1000.;
			fOutput << std::setprecision(9) << std::setw(20) << 0.;
			fOutput << std::setprecision(9) << std::setw(20) << T_(i);
			fOutput << std::setprecision(9) << std::setw(20) << P_(i);
			
			fOutput << std::setprecision(9) << std::setw(20) << epsilon_(i);
			fOutput << std::setprecision(9) << std::setw(20) << rho_bulk_(i);
			fOutput << std::setprecision(9) << std::setw(20) << Sv_(i);
			fOutput << std::setprecision(9) << std::setw(20) << rp_(i)*1.e6;
			fOutput << std::setprecision(9) << std::setw(20) << permeability_(i);
			fOutput << std::setprecision(9) << std::setw(20) << eta_bulk_(i);
			fOutput << std::setprecision(9) << std::setw(20) << eta_knudsen_(i);
			fOutput << std::setprecision(9) << std::setw(20) << eta_viscous_(i);

			fOutput << std::setprecision(9) << std::setw(20) << rho_gas_(i);
			fOutput << std::setprecision(9) << std::setw(20) << mw_(i);
			fOutput << std::setprecision(9) << std::setw(20) << Y_[i].sum();

			// Elements
			{
				for (unsigned int j = 0; j < thermodynamicsMap_.elements().size(); j++)
				{
					double sum = 0.;
					for (unsigned int k = 0; k < nc_; k++)
						sum += thermodynamicsMap_.atomic_composition()(k, j) * xx(k + 1);
					fOutput << std::setprecision(9) << std::setw(20) << sum;
				}
			}

			// Species (mole fractions and mass fractions)
			{
				for (unsigned int j = 0; j < nc_; j++)
					fOutput << std::setprecision(9) << std::setw(20) << xx(j + 1);
				for (unsigned int j = 0; j < nc_; j++)
					fOutput << std::setprecision(9) << std::setw(20) << yy(j + 1);
			}

			fOutput << std::endl;
		}

		fOutput.close();
	}

	void Reactor1D::PrintDiffusionCoefficients(const double t, const std::string name_file)
	{
		std::ofstream fOutput(name_file.c_str(), std::ios::out);

		{
			unsigned int count = 1;
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "time[s]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "x[mm]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "T[K]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "P[Pa]", count);

			for (unsigned int j = 0; j < nc_; j++)
			{
				OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, thermodynamicsMap_.NamesOfSpecies()[j] + "_ToE", count);
				OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, thermodynamicsMap_.NamesOfSpecies()[j] + "_FiE", count);
				OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, thermodynamicsMap_.NamesOfSpecies()[j] + "_KnE", count);
				OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, thermodynamicsMap_.NamesOfSpecies()[j] + "_Fi", count);
				OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, thermodynamicsMap_.NamesOfSpecies()[j] + "_Kn", count);
			}

			fOutput << std::endl;
		}

		for (unsigned int i = 0; i < np_; i++)
		{
			// Calculate the diffusion coefficients
			{
				porousMedium_.SetPorosity(epsilon_(i));
				porousMedium_.SetTemperature(T_(i));
				porousMedium_.SetPressure(P_(i));

				thermodynamicsMap_.SetPressure(P_(i));
				thermodynamicsMap_.SetTemperature(T_(i));
				aux_Y.CopyFrom(Y_[i].data());
				thermodynamicsMap_.MoleFractions_From_MassFractions(aux_X, mw_(i), aux_Y);
				aux_X.CopyTo(X_[i].data());
				for (unsigned int j = 0; j < nc_; j++)
					aux_eigen(j) = aux_X[j + 1];
				porousMedium_.EffectiveDiffusionCoefficients(aux_eigen);
			}

			fOutput << std::setprecision(9) << std::setw(20) << t;
			fOutput << std::setprecision(9) << std::setw(20) << grid_.x()[i] * 1000.;
			fOutput << std::setprecision(9) << std::setw(20) << T_(i);
			fOutput << std::setprecision(9) << std::setw(20) << P_(i);

			// Diffusion coefficients
			{
				for (unsigned int j = 0; j < nc_; j++)
				{
					fOutput << std::setprecision(9) << std::setw(20) << porousMedium_.gamma_effective()(j);
					fOutput << std::setprecision(9) << std::setw(20) << porousMedium_.gamma_fick_effective()(j);
					fOutput << std::setprecision(9) << std::setw(20) << porousMedium_.gamma_knudsen_effective()(j);
					fOutput << std::setprecision(9) << std::setw(20) << porousMedium_.gamma_fick()(j);
					fOutput << std::setprecision(9) << std::setw(20) << porousMedium_.gamma_knudsen()(j);
				}
			}

			fOutput << std::endl;
		}

		fOutput.close();
	}

	void Reactor1D::PrintHomogeneousRates(const double t, const std::string name_file)
	{
		std::ofstream fOutput(name_file.c_str(), std::ios::out);

		{
			unsigned int count = 1;
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "time[s]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "x[mm]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "T[K]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "P[Pa]", count);

			for (unsigned int j = 0; j < nc_; j++)
				OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, thermodynamicsMap_.NamesOfSpecies()[j] + "[kmol/m3/s]", count);

			fOutput << std::endl;
		}

		for (unsigned int i = 0; i < np_; i++)
		{
			// Formation rates
			{
				thermodynamicsMap_.SetTemperature(T_(i));
				thermodynamicsMap_.SetPressure(P_(i));
				aux_Y.CopyFrom(Y_[i].data());
				thermodynamicsMap_.MoleFractions_From_MassFractions(aux_X, mw_(i), aux_Y);
				aux_X.CopyTo(X_[i].data());
				const double cTot = P_(i) / PhysicalConstants::R_J_kmol / T_(i); // [kmol/m3]
				Product(cTot, aux_X, &aux_C);
				kineticsMap_.SetTemperature(T_(i));
				kineticsMap_.SetPressure(P_(i));
				kineticsMap_.ReactionRates(aux_C);
				kineticsMap_.FormationRates(&aux_R);
			}

			fOutput << std::setprecision(9) << std::setw(20) << t;
			fOutput << std::setprecision(9) << std::setw(20) << grid_.x()[i] * 1000.;
			fOutput << std::setprecision(9) << std::setw(20) << T_(i);
			fOutput << std::setprecision(9) << std::setw(20) << P_(i);

			for (unsigned int j = 0; j < nc_; j++)
					fOutput << std::setprecision(9) << std::setw(20) << aux_R[j+1];

			fOutput << std::endl;
		}

		fOutput.close();
	}

	void Reactor1D::PrintHeterogeneousRates(const double t, const std::string name_file)
	{
		if (detailed_heterogeneous_kinetics_ == true)
			PrintDetailedHeterogeneousRates(t, name_file);
		else
			PrintGlobalHeterogeneousRates(t, name_file);
	}

	void Reactor1D::PrintGlobalHeterogeneousRates(const double t, const std::string name_file)
	{
		unsigned int width = 25;

		std::ofstream fOutput(name_file.c_str(), std::ios::out);

		{
			unsigned int count = 1;
			OpenSMOKE::PrintTagOnASCIILabel(width, fOutput, "time[s]", count);
			OpenSMOKE::PrintTagOnASCIILabel(width, fOutput, "x[mm]", count);
			OpenSMOKE::PrintTagOnASCIILabel(width, fOutput, "T[K]", count);
			OpenSMOKE::PrintTagOnASCIILabel(width, fOutput, "P[Pa]", count);

			// Total deposition rate
			OpenSMOKE::PrintTagOnASCIILabel(width, fOutput, "rDepo[kmol/m2/s]", count);

			// Reaction rates 
			for (int j = 0; j < heterogeneousMechanism_.r().size(); j++)
			{
				std::stringstream number; number << j + 1;
				OpenSMOKE::PrintTagOnASCIILabel(width, fOutput, "r" + number.str() + "[kmol/m2/s]", count);
			}

			// Hydrogen inhibition factors 
			{
				OpenSMOKE::PrintTagOnASCIILabel(width, fOutput, "I_CH4[-]", count);
				OpenSMOKE::PrintTagOnASCIILabel(width, fOutput, "I_C2H4[-]", count);
				OpenSMOKE::PrintTagOnASCIILabel(width, fOutput, "I_C2H2[-]", count);
				OpenSMOKE::PrintTagOnASCIILabel(width, fOutput, "I_C6H6[-]", count);
			}

			// Heterogeneous formation rates [kmol/m3/s] 
			{
				OpenSMOKE::PrintTagOnASCIILabel(width, fOutput, "R_CH4[kmol/m3/s]", count);
				OpenSMOKE::PrintTagOnASCIILabel(width, fOutput, "R_C2H4[kmol/m3/s]", count);
				OpenSMOKE::PrintTagOnASCIILabel(width, fOutput, "R_C2H2[kmol/m3/s]", count);
				OpenSMOKE::PrintTagOnASCIILabel(width, fOutput, "R_C6H6[kmol/m3/s]", count);
				OpenSMOKE::PrintTagOnASCIILabel(width, fOutput, "R_H2[kmol/m3/s]", count);
			}

			fOutput << std::endl;
		}

		for (unsigned int i = 0; i < np_; i++)
		{
			// Calculate the reaction and formation rates of heterogeneous reactions
			{
				porousMedium_.SetPorosity(epsilon_(i));
				porousMedium_.SetTemperature(T_(i));
				porousMedium_.SetPressure(P_(i));

				thermodynamicsMap_.SetPressure(P_(i));
				thermodynamicsMap_.SetTemperature(T_(i));
				aux_Y.CopyFrom(Y_[i].data());
				thermodynamicsMap_.MoleFractions_From_MassFractions(aux_X, mw_(i), aux_Y);
				aux_X.CopyTo(X_[i].data());
				const double cTot = P_(i) / PhysicalConstants::R_J_kmol / T_(i); // [kmol/m3]
				Product(cTot, aux_X, &aux_C);

				heterogeneousMechanism_.SetTemperature(T_(i));
				heterogeneousMechanism_.SetPressure(P_(i));
				for (unsigned int j = 0; j < nc_; j++)
					aux_eigen(j) = aux_C[j + 1];
				heterogeneousMechanism_.FormationRates(porousMedium_.Sv(), aux_eigen);
			}

			fOutput << std::setprecision(9) << std::setw(width) << t;
			fOutput << std::setprecision(9) << std::setw(width) << grid_.x()[i] * 1000.;
			fOutput << std::setprecision(9) << std::setw(width) << T_(i);
			fOutput << std::setprecision(9) << std::setw(width) << P_(i);

			// Reaction rates
			{
				fOutput << std::setprecision(9) << std::setw(width) << heterogeneousMechanism_.r_deposition_per_unit_area();
				for (int j = 0; j < heterogeneousMechanism_.r().size(); j++)
					fOutput << std::setprecision(9) << std::setw(width) << heterogeneousMechanism_.r()(j);
			}

			// Hydrogen inhibition factors
			{
				fOutput << std::setprecision(9) << std::setw(width) << heterogeneousMechanism_.I_CH4();
				fOutput << std::setprecision(9) << std::setw(width) << heterogeneousMechanism_.I_C2H4();
				fOutput << std::setprecision(9) << std::setw(width) << heterogeneousMechanism_.I_C2H2();
				fOutput << std::setprecision(9) << std::setw(width) << heterogeneousMechanism_.I_C6H6();

			}

			// Heterogeneous formation rates
			{
				fOutput << std::setprecision(9) << std::setw(width) << heterogeneousMechanism_.Rgas()(heterogeneousMechanism_.index_CH4());
				fOutput << std::setprecision(9) << std::setw(width) << heterogeneousMechanism_.Rgas()(heterogeneousMechanism_.index_C2H4());
				fOutput << std::setprecision(9) << std::setw(width) << heterogeneousMechanism_.Rgas()(heterogeneousMechanism_.index_C2H2());
				fOutput << std::setprecision(9) << std::setw(width) << heterogeneousMechanism_.Rgas()(heterogeneousMechanism_.index_C6H6());
				fOutput << std::setprecision(9) << std::setw(width) << heterogeneousMechanism_.Rgas()(heterogeneousMechanism_.index_H2());
			}

			fOutput << std::endl;
		}

		fOutput.close();
	}

	void Reactor1D::PrintDetailedHeterogeneousRates(const double t, const std::string name_file)
	{
		unsigned int width = 25;

		std::ofstream fOutput(name_file.c_str(), std::ios::out);

		{
			unsigned int count = 1;
			OpenSMOKE::PrintTagOnASCIILabel(width, fOutput, "time[s]", count);
			OpenSMOKE::PrintTagOnASCIILabel(width, fOutput, "x[mm]", count);
			OpenSMOKE::PrintTagOnASCIILabel(width, fOutput, "T[K]", count);
			OpenSMOKE::PrintTagOnASCIILabel(width, fOutput, "P[Pa]", count);

			// Total deposition rate
			OpenSMOKE::PrintTagOnASCIILabel(width, fOutput, "rDepo[kmol/m2/s]", count);

			// Reaction rates 
		//	for (int j = 0; j < heterogeneousDetailedMechanism_.r().size(); j++)
		//	{
		//		std::stringstream number; number << j + 1;
		//		OpenSMOKE::PrintTagOnASCIILabel(width, fOutput, "r" + number.str() + "[kmol/m2/s]", count);
		//	}

			// Formation rates of homogeneous species from heterogeneous reactions [kmol/m3/s] 
			for (unsigned int i = 0; i < nc_; i++)
			{
				std::string label = "RV_" + thermodynamicsMap_.NamesOfSpecies()[i];
				OpenSMOKE::PrintTagOnASCIILabel(width, fOutput, label, count);
			}

			// Formation rates of surface species from heterogeneous reactions [kmol/m2/s] 
			for (unsigned int i = 0; i < surf_nc_; i++)
			{
				std::string label = "RS_" + thermodynamicsSurfaceMap_.NamesOfSpecies()[nc_+i];
				OpenSMOKE::PrintTagOnASCIILabel(width, fOutput, label, count);
			}

			fOutput << std::endl;
		}

		for (unsigned int i = 0; i < np_; i++)
		{
			// Calculate the reaction and formation rates of heterogeneous reactions
			{
				porousMedium_.SetPorosity(epsilon_(i));
				porousMedium_.SetTemperature(T_(i));
				porousMedium_.SetPressure(P_(i));

				// Homogeneous contributions
				{
					thermodynamicsMap_.SetPressure(P_(i));
					thermodynamicsMap_.SetTemperature(T_(i));
					aux_Y.CopyFrom(Y_[i].data());
					thermodynamicsMap_.MoleFractions_From_MassFractions(aux_X, mw_(i), aux_Y);
					aux_X.CopyTo(X_[i].data());
					const double cTot = P_(i) / PhysicalConstants::R_J_kmol / T_(i); // [kmol/m3]
					Product(cTot, aux_X, &aux_C);
				}

				// Heterogeneous mechanism
				{
					heterogeneousDetailedMechanism_.SetTemperature(T_(i));
					heterogeneousDetailedMechanism_.SetPressure(P_(i));

					for (unsigned int j = 0; j < nc_; j++)
						eigen_C_(j) = aux_C[j + 1];

					for (unsigned int j = 0; j < surf_nc_; j++)
						eigen_Z_(j) = Z_[i](j);

					for (unsigned int j = 0; j < bulk_nc_; j++)
						eigen_a_(j) = 1.;

					for (unsigned int j = 0; j < surf_np_; j++)
						eigen_gamma_(j) = Gamma_[i](j);

					heterogeneousDetailedMechanism_.FormationRates(Sv_(i), eigen_C_, eigen_Z_, eigen_a_, eigen_gamma_);
				}
			}

			fOutput << std::setprecision(9) << std::setw(width) << t;
			fOutput << std::setprecision(9) << std::setw(width) << grid_.x()[i] * 1000.;
			fOutput << std::setprecision(9) << std::setw(width) << T_(i);
			fOutput << std::setprecision(9) << std::setw(width) << P_(i);

			// Reaction rates
			{
				fOutput << std::setprecision(9) << std::setw(width) << heterogeneousDetailedMechanism_.r_deposition_per_unit_area();
		//		for (int j = 0; j < heterogeneousDetailedMechanism_.r().size(); j++)
		//			fOutput << std::setprecision(9) << std::setw(width) << heterogeneousDetailedMechanism_.r()(j);
			}

			// Formation rates of homogeneous species from heterogeneous reactions[kmol/m3/s]
			for (unsigned int i = 0; i < nc_; i++)
				fOutput << std::setprecision(9) << std::setw(width) << heterogeneousDetailedMechanism_.Rgas()(i);

			// Formation rates of surface species from heterogeneous reactions[kmol/m2/s]
			for (unsigned int i = 0; i < surf_nc_; i++)
				fOutput << std::setprecision(9) << std::setw(width) << heterogeneousDetailedMechanism_.Rsurface()(i);

			fOutput << std::endl;
		}

		fOutput.close();
	}

	void Reactor1D::PrintROPA(const double t, const std::string name_file)
	{
		// Open ROPA file
		std::ofstream fROPA(name_file.c_str(), std::ios::out);

		// Write head in file
		ropa_->WriteHead(fROPA, "Reactor1D");

		// ROPA at mouth section
		{
			OpenSMOKE::OpenSMOKEVectorDouble x(nc_);
			OpenSMOKE::OpenSMOKEVectorDouble omega(nc_);
			OpenSMOKE::OpenSMOKEVectorDouble c(nc_);
			OpenSMOKE::OpenSMOKEVectorDouble z(surf_nc_);
			OpenSMOKE::OpenSMOKEVectorDouble a(bulk_nc_);
			OpenSMOKE::OpenSMOKEVectorDouble gamma(surf_np_);

			std::vector<double> list_points(3);
			list_points[0] = 0;
			list_points[1] = np_/2;
			list_points[2] = np_-1;

			for (unsigned int i=0;i<list_points.size();i++)
			{
				const unsigned int point = static_cast<unsigned int>(list_points[i]);

				fROPA << std::endl;
				fROPA << "-------------------------------------------------------------------" << std::endl;
				fROPA << "ROPA at point: " << point << " - coordinate: " << grid_.x()[point] << std::endl;
				fROPA << "-------------------------------------------------------------------" << std::endl;
				fROPA << std::endl;

				// Molar fractions
				double mw;
				for (unsigned int j = 0; j < nc_; j++)
					omega[j + 1] = Y_[point](j);
				thermodynamicsMap_.MoleFractions_From_MassFractions(x, mw, omega);

				// Concentrations
				const double cTot = rho_gas_(point) / mw_(point);
				OpenSMOKE::Product(cTot, x, &c);

				// Surface species
				for (unsigned int j = 0; j < surf_nc_; j++)
					z[j + 1] = Z_[point](j);

				// Bulk activities
				for (unsigned int j = 0; j < bulk_nc_; j++)
					a[j + 1] = 1.;

				// Surface densities
				for (unsigned int j = 0; j < surf_np_; j++)
					gamma[j + 1] = Gamma_[point](j);

				ropa_->Analyze(fROPA, 0, 0., T_(point), P_(point), c, omega, z, a, gamma);
			}
		}
	}

	void Reactor1D::Print(const double t, const double* y)
	{
		if (count_dae_video_%n_steps_video_ == 1)
		{
			if (count_dae_video_ % (n_steps_video_ * 50) == 1)
			{
				std::cout << std::endl;
				std::cout << std::left << std::setw(14) << "Time[s]";		// [s]
				std::cout << std::left << std::setw(14) << "Time[h]";		// [h]
				std::cout << std::left << std::setw(14) << "Rho(SP)[kg/m3]";	// [kg/m3]
				std::cout << std::left << std::setw(14) << "Rho(GS)[kg/m3]";	// [kg/m3]
				std::cout << std::left << std::setw(14) << "eps(SP)[-]";	// [-]
				std::cout << std::left << std::setw(14) << "eps(GS)[-]";	// [-]
				std::cout << std::left << std::setw(14) << "eps(min)[-]";	// [-]
				std::cout << std::left << std::setw(14) << "Rdep[mu/h]";	// [micron/h]
				std::cout << std::left << std::setw(16) << "Rdep[kg/m3/h]";	// [kg/m3/h]

				if (detailed_heterogeneous_kinetics_ == true)
				{
					for (unsigned int i = 0; i < surf_np_; i++)
						std::cout << std::left << std::setw(16) << "Gamma[kmol/m2]";
					std::cout << std::left << std::setw(14) << "Error";
				}

				std::cout << std::endl;
			}

			const double r_deposition_per_unit_area_mean = AreaAveraged(omega_deposition_per_unit_area_);		// [kg/m2/s]
			const double r_deposition_per_unit_volume_mean = AreaAveraged(omega_deposition_per_unit_volume_);	// [kg/m3/s]

			std::cout << std::left << std::setw(14) << std::scientific << std::setprecision(6) << t;									// [s]
			std::cout << std::left << std::setw(14) << std::scientific << std::setprecision(6) << t / 3600.;							// [h]
			std::cout << std::left << std::setw(14) << std::fixed << std::setprecision(6) << rho_bulk_(0);		// [kg/m3]
			std::cout << std::left << std::setw(14) << std::fixed << std::setprecision(6) << rho_bulk_(np_-1);	// [kg/m3]
			std::cout << std::left << std::setw(14) << std::fixed << std::setprecision(6) << epsilon_(0);		// [-]
			std::cout << std::left << std::setw(14) << std::fixed << std::setprecision(6) << epsilon_(np_-1);	// [-]
			std::cout << std::left << std::setw(14) << std::fixed << std::setprecision(6) << epsilon_.minCoeff();	// [-]
			std::cout << std::left << std::setw(14) << std::fixed << std::setprecision(6) << r_deposition_per_unit_area_mean / rho_graphite_*1e6*3600.;	// [micron/h]
			std::cout << std::left << std::setw(16) << std::fixed << std::setprecision(6) << r_deposition_per_unit_volume_mean *3600.;					// [kg/m3/h]

			if (detailed_heterogeneous_kinetics_ == true)
			{
				for (unsigned int i = 0; i < surf_np_; i++)
					std::cout << std::left << std::setw(16) << std::scientific << Gamma_[0](i);

				Eigen::VectorXd error_z_sum(np_);
				for (unsigned int i = 0; i < np_; i++)
					error_z_sum(i) = std::fabs(Z_[i].sum() - 1.);

				std::cout << std::left << std::setw(14) << std::setprecision(2) << std::scientific << error_z_sum.maxCoeff();
			}

			std::cout << std::endl;
		}

		if (count_file_ == n_steps_file_)
		{
			const int width = 20;

			const double T_mean = AreaAveraged(T_);
			const double T_std = AreaStandardDeviation(T_mean, T_);

			const double P_mean = AreaAveraged(P_);
			const double P_std = AreaStandardDeviation(P_mean, P_);

			const double epsilon_mean = AreaAveraged(epsilon_);
			const double epsilon_std = AreaStandardDeviation(epsilon_mean, epsilon_);

			const double rho_bulk_mean = AreaAveraged(rho_bulk_);
			const double rho_bulk_std = AreaStandardDeviation(rho_bulk_mean, rho_bulk_);

			const double Sv_mean = AreaAveraged(Sv_);
			const double Sv_std = AreaStandardDeviation(Sv_mean, Sv_);

			const double rp_mean = AreaAveraged(rp_);
			const double rp_std = AreaStandardDeviation(rp_mean, rp_);

			const double permeability_mean = AreaAveraged(permeability_);
			const double permeability_std = AreaStandardDeviation(permeability_mean, permeability_);

			const double eta_bulk_mean = AreaAveraged(eta_bulk_);
			const double eta_bulk_std = AreaStandardDeviation(eta_bulk_mean, eta_bulk_);

			const double eta_knudsen_mean = AreaAveraged(eta_knudsen_);
			const double eta_knudsen_std = AreaStandardDeviation(eta_knudsen_mean, eta_knudsen_);

			const double eta_viscous_mean = AreaAveraged(eta_viscous_);
			const double eta_viscous_std = AreaStandardDeviation(eta_viscous_mean, eta_viscous_);

			const double r_deposition_per_unit_area_mean = AreaAveraged(omega_deposition_per_unit_area_);
			const double r_deposition_per_unit_area_std = AreaStandardDeviation(r_deposition_per_unit_area_mean, omega_deposition_per_unit_area_);

			const double r_deposition_per_unit_volume_mean = AreaAveraged(omega_deposition_per_unit_volume_);
			const double r_deposition_per_unit_volume_std = AreaStandardDeviation(r_deposition_per_unit_volume_mean, omega_deposition_per_unit_volume_);

			fMonitoring_ << std::left << std::setprecision(9) << std::setw(width) << t / 3600.;	// [h]
			fMonitoring_ << std::left << std::setprecision(9) << std::setw(width) << t;			// [s]

			fMonitoring_ << std::left << std::setw(width) << std::fixed << std::setprecision(4) << T_mean;		// [K]
			fMonitoring_ << std::left << std::setw(width) << std::fixed << std::setprecision(4) << T_std;		// [K]

			fMonitoring_ << std::left << std::setw(width) << std::fixed << std::setprecision(4) << P_mean;		// [Pa]
			fMonitoring_ << std::left << std::setw(width) << std::fixed << std::setprecision(4) << P_std;		// [Pa]

			fMonitoring_ << std::left << std::setw(width) << std::fixed << std::setprecision(4) << epsilon_mean;
			fMonitoring_ << std::left << std::setw(width) << std::fixed << std::setprecision(4) << epsilon_std;

			fMonitoring_ << std::left << std::setw(width) << std::fixed << std::setprecision(4) << rho_bulk_mean;		// [kg/m3]
			fMonitoring_ << std::left << std::setw(width) << std::fixed << std::setprecision(4) << rho_bulk_std;		// [kg/m3]

			fMonitoring_ << std::left << std::setw(width) << std::fixed << std::setprecision(4) << Sv_mean;	// [1/m]
			fMonitoring_ << std::left << std::setw(width) << std::fixed << std::setprecision(4) << Sv_std;		// [1/m]

			fMonitoring_ << std::left << std::setw(width) << std::fixed << std::setprecision(4) << rp_mean*1e6;	// [micron]
			fMonitoring_ << std::left << std::setw(width) << std::fixed << std::setprecision(4) << rp_std*1e6;		// [micron]

			fMonitoring_ << std::left << std::setw(width) << std::fixed << std::setprecision(4) << permeability_mean;	// [m2]
			fMonitoring_ << std::left << std::setw(width) << std::fixed << std::setprecision(4) << permeability_std;	// [m2]

			fMonitoring_ << std::left << std::setw(width) << std::fixed << std::setprecision(4) << eta_bulk_mean;
			fMonitoring_ << std::left << std::setw(width) << std::fixed << std::setprecision(4) << eta_bulk_std;

			fMonitoring_ << std::left << std::setw(width) << std::fixed << std::setprecision(4) << eta_knudsen_mean;
			fMonitoring_ << std::left << std::setw(width) << std::fixed << std::setprecision(4) << eta_knudsen_std;

			fMonitoring_ << std::left << std::setw(width) << std::fixed << std::setprecision(4) << eta_viscous_mean;
			fMonitoring_ << std::left << std::setw(width) << std::fixed << std::setprecision(4) << eta_viscous_std;

			fMonitoring_ << std::left << std::setw(width) << std::scientific << std::setprecision(6) << r_deposition_per_unit_area_mean / rho_graphite_*1000.;		// [mm/s]
			fMonitoring_ << std::left << std::setw(width) << std::scientific << std::setprecision(6) << r_deposition_per_unit_area_std / rho_graphite_*1000.;	// [mm/s]

			fMonitoring_ << std::left << std::setw(width) << std::scientific << std::setprecision(6) << r_deposition_per_unit_area_mean *1000.*3600.;	// [g/m2/h]
			fMonitoring_ << std::left << std::setw(width) << std::scientific << std::setprecision(6) << r_deposition_per_unit_area_std *1000.*3600.;		// [g/m2/h]

			fMonitoring_ << std::left << std::setw(width) << std::scientific << std::setprecision(6) << r_deposition_per_unit_volume_mean *1000.*3600.;	// [g/m3/h]
			fMonitoring_ << std::left << std::setw(width) << std::scientific << std::setprecision(6) << r_deposition_per_unit_volume_std *1000.*3600.;	// [g/m3/h]

			fMonitoring_ << std::endl;

			count_file_ = 0;
		}

		count_dae_video_++;
		count_file_++;
	}

	void Reactor1D::SparsityPattern(std::vector<unsigned int>& rows, std::vector<unsigned int>& cols)
	{
		const unsigned int dimblock = BlockDimensions();
		const unsigned int nblocks = NumberOfEquations() / dimblock;
		const unsigned int nonzeros = (dimblock*dimblock*(3 * nblocks - 2));

		rows.resize(nonzeros);
		cols.resize(nonzeros);

		int count = 0;

		// First block (point)
		{
			for (unsigned int ki = 1; ki <= dimblock; ki++)
			for (unsigned int kj = 1; kj <= 2 * dimblock; kj++)
			{
				rows[count] = ki;
				cols[count] = kj;
				count++;
			}
		}

		// Internal blocks (points)
		for (unsigned int block = 2; block <= nblocks - 1; block++)
		{
			const unsigned int i = (block - 1)*dimblock;
			const unsigned int j = (block - 2)*dimblock;

			for (unsigned int ki = 1; ki <= dimblock; ki++)
			for (unsigned int kj = 1; kj <= 3 * dimblock; kj++)
			{
				rows[count] = i + ki;
				cols[count] = j + kj;
				count++;
			}
		}

		// Last block (point)
		{
			const unsigned int i = (nblocks - 1)*dimblock;
			const unsigned int j = (nblocks - 2)*dimblock;

			for (unsigned int ki = 1; ki <= dimblock; ki++)
			for (unsigned int kj = 1; kj <= 2 * dimblock; kj++)
			{
				rows[count] = i + ki;
				cols[count] = j + kj;
				count++;
			}
		}
	}

	double Reactor1D::AreaAveraged(const Eigen::VectorXd& v)
	{
		if (planar_symmetry_ == true)
		{
			double sum = 0.;

			// Internal
			for (unsigned int i = 1; i < np_ - 1; i++)
			{
				const double area = grid_.dxc()(i)*0.50;
				sum += v(i)*area;
			}

			// West (zero gradient)
			{
				const double area = grid_.dxe()(0)*0.50;
				sum += v(0)*area;
			}

			// East (gas side)
			{
				const double area = grid_.dxw()(np_ - 1)*0.50;
				sum += v(np_ - 1)*area;
			}


			const double area_total = grid_.L();
			return sum / area_total;
		}
		else
		{
			double sum = 0.;

			// Internal
			for (unsigned int i = 1; i < np_ - 1; i++)
			{
				const double ri = (grid_.x()(i) + grid_.x()(i - 1)) / 2.;
				const double re = (grid_.x()(i) + grid_.x()(i + 1)) / 2.;
				const double area = boost::math::constants::pi<double>()*(re*re - ri*ri);
				sum += v(i)*area;
			}

			// West (zero gradient)
			{
				const double ri = grid_.x()(0);
				const double re = (grid_.x()(0) + grid_.x()(1)) / 2.;
				const double area = boost::math::constants::pi<double>()*(re*re - ri*ri);
				sum += v(0)*area;
			}

			// East (gas side)
			{
				const double ri = (grid_.x()(np_-1) + grid_.x()(np_ - 2)) / 2.;
				const double re = grid_.x()(np_-1);
				const double area = boost::math::constants::pi<double>()*(re*re - ri*ri);
				sum += v(np_ - 1)*area;
			}

			const double Ri = grid_.x()(0);
			const double Re = grid_.x()(np_ - 1);
			const double area_total = boost::math::constants::pi<double>()*(Re*Re - Ri*Ri);
			return sum / area_total;
		}
	}

	double Reactor1D::AreaStandardDeviation(const double mean, const Eigen::VectorXd& v)
	{
		if (planar_symmetry_ == true)
		{
			double sum_std = 0.;

			// Internal
			for (unsigned int i = 1; i < np_ - 1; i++)
			{
				const double area = grid_.dxc()(i)*0.50;
				sum_std += boost::math::pow<2>(v(i) - mean)*area;
			}

			// West (zero gradient)
			{
				const double area = grid_.dxe()(0)*0.50;
				sum_std += boost::math::pow<2>(v(0) - mean)*area;
			}

			// East (gas side)
			{
				const double area = grid_.dxw()(np_ - 1)*0.50;
				sum_std += boost::math::pow<2>(v(np_ - 1) - mean)*area;
			}

			const double coefficient = double(np_ - 1) / double(np_);
			const double area_total = grid_.L();
			return std::sqrt(sum_std / area_total / coefficient);
		}
		else
		{
			double sum_std = 0.;

			// Internal
			for (unsigned int i = 1; i < np_ - 1; i++)
			{
				const double ri = (grid_.x()(i) + grid_.x()(i - 1)) / 2.;
				const double re = (grid_.x()(i) + grid_.x()(i + 1)) / 2.;
				const double area = boost::math::constants::pi<double>()*(re*re - ri*ri);
				sum_std += boost::math::pow<2>(v(np_ - 1) - mean)*area;
			}

			// West (zero gradient)
			{
				const double ri = grid_.x()(0);
				const double re = (grid_.x()(0) + grid_.x()(1)) / 2.;
				const double area = boost::math::constants::pi<double>()*(re*re - ri*ri);
				sum_std += boost::math::pow<2>(v(np_ - 1) - mean)*area;
			}

			// East (gas side)
			{
				const double ri = (grid_.x()(np_ - 1) + grid_.x()(np_ - 2)) / 2.;
				const double re = grid_.x()(np_ - 1);
				const double area = boost::math::constants::pi<double>()*(re*re - ri*ri);
				sum_std += boost::math::pow<2>(v(np_ - 1) - mean)*area;
			}

			const double Ri = grid_.x()(0);
			const double Re = grid_.x()(np_ - 1);
			const double area_total = boost::math::constants::pi<double>()*(Re*Re - Ri*Ri);
			const double coefficient = double(np_ - 1) / double(np_);
			return std::sqrt(sum_std / area_total / coefficient);
		}
	}

	void Reactor1D::PrintLabelMonitoringFile()
	{
		unsigned int count = 1;
		const unsigned int width = 20;
		OpenSMOKE::PrintTagOnASCIILabel(width, fMonitoring_, "time[h]", count);
		OpenSMOKE::PrintTagOnASCIILabel(width, fMonitoring_, "time[s]", count);
		OpenSMOKE::PrintTagOnASCIILabel(width, fMonitoring_, "T[K]", count);
		OpenSMOKE::PrintTagOnASCIILabel(width, fMonitoring_, "Tstd[K]", count);
		OpenSMOKE::PrintTagOnASCIILabel(width, fMonitoring_, "Pa[Pa]", count);
		OpenSMOKE::PrintTagOnASCIILabel(width, fMonitoring_, "Pastd[Pa]", count);
		OpenSMOKE::PrintTagOnASCIILabel(width, fMonitoring_, "eps[-]", count);
		OpenSMOKE::PrintTagOnASCIILabel(width, fMonitoring_, "epsstd[-]", count);
		OpenSMOKE::PrintTagOnASCIILabel(width, fMonitoring_, "rhoB[kg/m3]", count);
		OpenSMOKE::PrintTagOnASCIILabel(width, fMonitoring_, "rhoBstd[kg/m3]", count);
		OpenSMOKE::PrintTagOnASCIILabel(width, fMonitoring_, "Sv[1/]", count);
		OpenSMOKE::PrintTagOnASCIILabel(width, fMonitoring_, "Svstd[1/m]", count);
		OpenSMOKE::PrintTagOnASCIILabel(width, fMonitoring_, "rp[micron]", count);
		OpenSMOKE::PrintTagOnASCIILabel(width, fMonitoring_, "rpstd[micron]", count);
		OpenSMOKE::PrintTagOnASCIILabel(width, fMonitoring_, "K[m2]", count);
		OpenSMOKE::PrintTagOnASCIILabel(width, fMonitoring_, "Kstd[m2]", count);
		OpenSMOKE::PrintTagOnASCIILabel(width, fMonitoring_, "etaBulk[-]", count);
		OpenSMOKE::PrintTagOnASCIILabel(width, fMonitoring_, "etaBulkstd[K]", count);
		OpenSMOKE::PrintTagOnASCIILabel(width, fMonitoring_, "etaKn[-]", count);
		OpenSMOKE::PrintTagOnASCIILabel(width, fMonitoring_, "etaKnstd[K]", count);
		OpenSMOKE::PrintTagOnASCIILabel(width, fMonitoring_, "etaVisc[-]", count);
		OpenSMOKE::PrintTagOnASCIILabel(width, fMonitoring_, "etaViscstd[K]", count);
		OpenSMOKE::PrintTagOnASCIILabel(width, fMonitoring_, "rDep[mm/s]", count);
		OpenSMOKE::PrintTagOnASCIILabel(width, fMonitoring_, "rDepstd[mm/s]", count);
		OpenSMOKE::PrintTagOnASCIILabel(width, fMonitoring_, "rDep[g/m2/h]", count);
		OpenSMOKE::PrintTagOnASCIILabel(width, fMonitoring_, "rDepstd[g/m2/h]", count);
		OpenSMOKE::PrintTagOnASCIILabel(width, fMonitoring_, "rDep[g/m3/h]", count);
		OpenSMOKE::PrintTagOnASCIILabel(width, fMonitoring_, "rDepstd[g/m3/h]", count);

		fMonitoring_ << std::endl;
	}
}
