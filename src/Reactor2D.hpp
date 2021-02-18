/*-----------------------------------------------------------------------*\
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
#include "Interface_Reactor2D_OpenSMOKEppDae.h"

// Pointer
CVI::Reactor2D* reactor2d;

// Interfaces (Sundials)
#if OPENSMOKE_USE_SUNDIALS == 1
#include "sundials_header.h"
#include "Interface_Reactor2D_Ida.h"
#endif

// Interfaces (BzzMath)
#if OPENSMOKE_USE_BZZMATH == 1
#include "BzzMath.hpp"
#include "Interface_Reactor2D_BzzDae.h"
#endif

// Interfaces (DASPK)
#if OPENSMOKE_USE_DASPK == 1
#include "Interface_Reactor2D_Daspk.h"
#endif

// Boost C++
#include <boost/math/special_functions/pow.hpp>

namespace CVI
{
	Reactor2D::Reactor2D(	OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap,
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
							const boost::filesystem::path output_folder) :

	thermodynamicsMap_(thermodynamicsMap),
	kineticsMap_(kineticsMap),
	transportMap_(transportMap),
	thermodynamicsSurfaceMap_(thermodynamicsSurfaceMap),
	kineticsSurfaceMap_(kineticsSurfaceMap),
	porousMedium_(porousMedium),
	porosityDefect_(porosityDefect),
	heterogeneousMechanism_(heterogeneousMechanism),
	heterogeneousDetailedMechanism_(heterogeneousDetailedMechanism),
	grid_x_(grid_x),
	grid_y_(grid_y),
	plugFlowReactor_(plugFlowReactor),
	detailed_heterogeneous_kinetics_(detailed_heterogeneous_kinetics),
	site_non_conservation_(site_non_conservation),
	output_folder_(output_folder)
	{
		t_old_ = 0.;
		time_smoothing_ = 10.;

		gaseous_phase_ = GASEOUS_PHASE_FROM_PLUG_FLOW;
		equations_set_ = EQUATIONS_SET_COMPLETE;

		n_steps_video_ = 10;
		count_dae_video_ = 1;
		count_ode_video_ = 1;
		n_steps_file_ = 3;
		n_steps_update_plug_flow_ = 30;
		count_file_ = n_steps_file_;
		count_tecplot_ = 0;
		count_update_plug_flow_ = 0;
		ropa_analysis_ = false;

		time_total_ = 48.*3600.;
		dae_time_interval_ = 3600.;
		tecplot_time_interval_ = 3600.;
		ode_end_time_ = 1.;
		time_profiles_ = false;

		time_starting_point_ = 0.;
		start_from_backup_ = false;

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
			surface_dae_species_index_ = thermodynamicsSurfaceMap_.IndexOfSpecies(surface_dae_species) - (nc_ + 1);
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
		
		nx_ = grid_x_.Np();
		ny_ = grid_y_.Np();
		np_ = nx_*ny_;
		ne_ = block_*np_;
		band_size_ = block_*(nx_+1)-1;

		planar_symmetry_ = true;
		hole_ = true;

		gas_dae_species_index_ = thermodynamicsMap_.IndexOfSpecies(gas_dae_species)-1;

		vx_ = 0.;
		vy_ = 0.;

		output_tecplot_folder_ = output_folder_ / "Tecplot";
		OpenSMOKE::CreateDirectory(output_tecplot_folder_);

		output_plug_flow_folder_ = output_folder_ / "PlugFlow";
		OpenSMOKE::CreateDirectory(output_plug_flow_folder_);

		output_matlab_folder_ = output_folder_ / "Matlab";
		OpenSMOKE::CreateDirectory(output_matlab_folder_);

		output_diffusion_folder_ = output_folder_ / "DiffusionCoefficients";
		OpenSMOKE::CreateDirectory(output_diffusion_folder_);

		output_heterogeneous_folder_ = output_folder_ / "HeterogeneousReactions";
		OpenSMOKE::CreateDirectory(output_heterogeneous_folder_);

		output_homogeneous_folder_ = output_folder_ / "HomogeneousReactions";
		OpenSMOKE::CreateDirectory(output_homogeneous_folder_);

		output_disks_source_terms_folder_ = output_folder_ / "DisksSourceTerms";
		OpenSMOKE::CreateDirectory(output_disks_source_terms_folder_);

		output_ropa_folder_ = output_folder_ / "ROPA";
		OpenSMOKE::CreateDirectory(output_ropa_folder_);

		output_backup_folder_ = output_folder_ / "Backup";
		OpenSMOKE::CreateDirectory(output_backup_folder_);

		fMonitoring_.open((output_folder_ / "monitor.out").string().c_str(), std::ios::out);
		fMonitoring_.setf(std::ios::scientific);
		PrintLabelMonitoringFile();

		MemoryAllocation();
		SetAlgebraicDifferentialEquations();

		list_points_south_.resize(nx_);
		list_points_north_.resize(nx_);
		list_points_east_.resize(ny_);
		list_points_west_.resize(ny_);

		for (unsigned int i = 0; i < nx_; i++)
			list_points_south_(i) = i;

		for (unsigned int i = 0; i < nx_; i++)
			list_points_north_(i) = nx_*(ny_-1) + i;

		for (unsigned int i = 0; i < ny_; i++)
			list_points_east_(i) = (i+1)*nx_-1;

		for (unsigned int i = 0; i < ny_; i++)
			list_points_west_(i) = i*nx_;
	}

	void Reactor2D::MemoryAllocation()
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

		// Time derivatives: temperature [K/s]
		dT_over_dt_.resize(np_);

		// Time derivatives: mass fractions [1/s]
		dY_over_dt_.resize(np_);
		for (unsigned int i = 0; i < np_; i++)
			dY_over_dt_[i].resize(nc_);

		// Time derivatives: porosity [1/s]
		depsilon_over_dt_.resize(np_);


		// Gas side Mass fractions [-]
		Y_gas_north_side_.resize(nx_);
		Y_gas_south_side_.resize(nx_);
		for (unsigned int i = 0; i < nx_; i++)
		{
			Y_gas_north_side_[i].resize(nc_);
			Y_gas_south_side_[i].resize(nc_);
		}

		Y_gas_east_side_.resize(ny_);
		Y_gas_west_side_.resize(ny_);
		for (unsigned int i = 0; i < ny_; i++)
		{
			Y_gas_east_side_[i].resize(nc_);
			Y_gas_west_side_[i].resize(nc_);
		}

		// Gas side temperature and pressure
		T_gas_north_side_.resize(nx_);
		P_gas_north_side_.resize(nx_);
		T_gas_south_side_.resize(nx_);
		P_gas_south_side_.resize(nx_);
		T_gas_east_side_.resize(ny_);
		P_gas_east_side_.resize(ny_);
		T_gas_west_side_.resize(ny_);
		P_gas_west_side_.resize(ny_);

		// Increment of bulk density due to the single reactions
		delta_rhobulk_.resize(np_);
		delta_rhobulk_due_to_single_reaction_.resize(heterogeneousMechanism_.r().size());
		delta_rhobulk_due_to_single_reaction_over_rhobulk_.resize(heterogeneousMechanism_.r().size());
		for (int i = 0; i < heterogeneousMechanism_.r().size(); i++)
		{
			delta_rhobulk_due_to_single_reaction_[i].resize(np_);
			delta_rhobulk_due_to_single_reaction_over_rhobulk_[i].resize(np_);
			delta_rhobulk_due_to_single_reaction_[i].setZero();
			delta_rhobulk_due_to_single_reaction_over_rhobulk_[i].setZero();
		}

		// Reset to zero
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

		// Total mass produced/consumed (kg)
		homogeneous_total_mass_source_.resize(nc_);
		heterogeneous_total_mass_source_.resize(nc_);
		total_mass_exchanged_.resize(nc_);
		homogeneous_total_mass_source_.setZero();
		heterogeneous_total_mass_source_.setZero();
		total_mass_exchanged_.setZero();
	}

	void Reactor2D::SetSurfaceOnTheFlyROPA(OpenSMOKE::SurfaceOnTheFlyROPA* ropa)
	{
		ropa_ = ropa;
		ropa_analysis_ = true;
	}

	void Reactor2D::SetPlanarSymmetry(const bool flag)
	{
		planar_symmetry_ = flag;
		if (std::fabs(grid_x_.x()(0)) <= 1e-12)
			hole_ = false;
	}

	void Reactor2D::SetSiteNonConservation(std::vector<bool>& site_non_conservation)
	{
		site_non_conservation_ = site_non_conservation;
	}

	void Reactor2D::SetAlgebraicDifferentialEquations()
	{
		id_equations_.resize(ne_);

		unsigned int count = 0;

		if (equations_set_ == EQUATIONS_SET_COMPLETE)
		{
			// South side
			for (unsigned int i = 0; i < nx_; i++)
			{
				// Species
				for (unsigned int j = 0; j < nc_; j++)
					id_equations_[count++] = false;

				// Porosity
				id_equations_[count++] = true;

				// Surface species
				for (unsigned int j = 0; j < surf_np_; j++)
					id_equations_[count++] = true;

				if (dae_formulation_ == false)
				{
					for (unsigned int j = 0; j < surf_nc_; j++)
						id_equations_[count++] = true;
				}
				else
				{
					for (unsigned int j = 0; j < surface_dae_species_index_; j++)
						id_equations_[count++] = true;
					id_equations_[count++] = false;
					for (unsigned int j = surface_dae_species_index_ + 1; j < surf_nc_; j++)
						id_equations_[count++] = true;
				}
			}


			// Internal strips
			for (unsigned int j = 1; j < ny_ - 1; j++)
			{
				// West side
				{
					// Species
					for (unsigned int j = 0; j < nc_; j++)
						id_equations_[count++] = false;

					// Porosity
					id_equations_[count++] = true;

					// Surface species
					for (unsigned int j = 0; j < surf_np_; j++)
						id_equations_[count++] = true;

					if (dae_formulation_ == false)
					{
						for (unsigned int j = 0; j < surf_nc_; j++)
							id_equations_[count++] = true;
					}
					else
					{
						for (unsigned int j = 0; j < surface_dae_species_index_; j++)
							id_equations_[count++] = true;
						id_equations_[count++] = false;
						for (unsigned int j = surface_dae_species_index_ + 1; j < surf_nc_; j++)
							id_equations_[count++] = true;
					}
				}

				// Internal point
				for (unsigned int i = 1; i < nx_ - 1; i++)
				{
					// Species
					for (unsigned int j = 0; j < gas_dae_species_index_; j++)
						id_equations_[count++] = true;
					id_equations_[count++] = false;
					for (unsigned int j = gas_dae_species_index_ + 1; j < nc_; j++)
						id_equations_[count++] = true;

					// Porosity
					id_equations_[count++] = true;

					// Surface species
					for (unsigned int j = 0; j < surf_np_; j++)
						id_equations_[count++] = true;

					if (dae_formulation_ == false)
					{
						for (unsigned int j = 0; j < surf_nc_; j++)
							id_equations_[count++] = true;
					}
					else
					{
						for (unsigned int j = 0; j < surface_dae_species_index_; j++)
							id_equations_[count++] = true;
						id_equations_[count++] = false;
						for (unsigned int j = surface_dae_species_index_ + 1; j < surf_nc_; j++)
							id_equations_[count++] = true;
					}
				}

				// East side (gas side)
				{
					// Species
					for (unsigned int j = 0; j < nc_; j++)
						id_equations_[count++] = false;

					// Porosity
					id_equations_[count++] = true;

					// Surface species
					for (unsigned int j = 0; j < surf_np_; j++)
						id_equations_[count++] = true;

					if (dae_formulation_ == false)
					{
						for (unsigned int j = 0; j < surf_nc_; j++)
							id_equations_[count++] = true;
					}
					else
					{
						for (unsigned int j = 0; j < surface_dae_species_index_; j++)
							id_equations_[count++] = true;
						id_equations_[count++] = false;
						for (unsigned int j = surface_dae_species_index_ + 1; j < surf_nc_; j++)
							id_equations_[count++] = true;
					}
				}
			}

			// North side
			for (unsigned int i = 0; i < nx_; i++)
			{
				// Species
				for (unsigned int j = 0; j < nc_; j++)
					id_equations_[count++] = false;

				// Porosity
				id_equations_[count++] = true;

				// Surface species
				for (unsigned int j = 0; j < surf_np_; j++)
					id_equations_[count++] = true;

				if (dae_formulation_ == false)
				{
					for (unsigned int j = 0; j < surf_nc_; j++)
						id_equations_[count++] = true;
				}
				else
				{
					for (unsigned int j = 0; j < surface_dae_species_index_; j++)
						id_equations_[count++] = true;
					id_equations_[count++] = false;
					for (unsigned int j = surface_dae_species_index_ + 1; j < surf_nc_; j++)
						id_equations_[count++] = true;
				}
			}
		}
		else if (equations_set_ == EQUATIONS_SET_ONLYTEMPERATURE)
		{
			// South side
			for (unsigned int i = 0; i < nx_; i++)
				id_equations_[count++] = false;

			// Internal strips
			for (unsigned int j = 1; j < ny_ - 1; j++)
			{
				// West side
				id_equations_[count++] = false;

				// Internal point
				for (unsigned int i = 1; i < nx_ - 1; i++)
					id_equations_[count++] = true;

				// East side (gas side)
				id_equations_[count++] = false;
			}

			// North side
			for (unsigned int i = 0; i < nx_; i++)
				id_equations_[count++] = false;
		}

		// Data needed by Sundials Ida
		{
			unsigned int n_algebraic = std::count_if(id_equations_.begin(), id_equations_.end(), std::bind2nd(std::equal_to<bool>(), false));
			unsigned int n_differential = std::count_if(id_equations_.begin(), id_equations_.end(), std::bind2nd(std::equal_to<bool>(), true));

			algebraic_equations_.resize(n_algebraic);
			differential_equations_.resize(n_differential);

			int count_differential = 0;
			int count_algebraic = 0;
			for (unsigned int i = 0; i < id_equations_.size(); i++)
			{
				if (id_equations_[i] == true)  differential_equations_(count_differential++) = i;
				if (id_equations_[i] == false) algebraic_equations_(count_algebraic++) = i;
			}
		}
	}

	void Reactor2D::SetGasSide(const double T_gas, const double P_gas, const std::vector<Eigen::VectorXd>& omega_gas)
	{
		time_profiles_ = false;

		gaseous_phase_ = GASEOUS_PHASE_FROM_PLUG_FLOW;

		// Sides
		{
			// South side
			for (unsigned int i = 0; i < nx_; i++)
			{
				for (unsigned int j = 0; j < nc_; j++)
					Y_gas_south_side_[i](j) = omega_gas[i](j);

				P_gas_south_side_(i) = P_gas;
				T_gas_south_side_(i) = T_gas;
			}

			// East side
			for (unsigned int i = 0; i < ny_; i++)
			{
				for (unsigned int j = 0; j < nc_; j++)
					Y_gas_east_side_[i](j) = omega_gas[nx_ + i](j);

				P_gas_east_side_(i) = P_gas;
				T_gas_east_side_(i) = T_gas;
			}

			// North side
			for (unsigned int i = 0; i < nx_; i++)
			{
				for (unsigned int j = 0; j < nc_; j++)
					Y_gas_north_side_[i](j) = omega_gas[2 * nx_ + ny_ - 1 - i](j);

				P_gas_north_side_(i) = P_gas;
				T_gas_north_side_(i) = T_gas;
			}
		}

		// Force consistency
		{
			// Set initial fields consistent along the east side
			if (plugFlowReactor_.geometric_pattern() == CVI::PlugFlowReactorCoupled::GeometricPattern::ONE_SIDE)
			{
				for (unsigned int i = 0; i < ny_; i++)
				{
					for (unsigned int j = 0; j < nc_; j++)
						Y_[list_points_east_(i)](j) = omega_gas[nx_ + i](j);

					P_(list_points_east_(i)) = P_gas;
					T_(list_points_east_(i)) = T_gas;
				}
			}
			// Set initial fields consistent along the south, east, and north sides
			else if (plugFlowReactor_.geometric_pattern() == CVI::PlugFlowReactorCoupled::GeometricPattern::THREE_SIDES)
			{
				for (unsigned int i = 0; i < nx_; i++)
				{
					for (unsigned int j = 0; j < nc_; j++)
						Y_[list_points_south_(i)](j) = omega_gas[i](j);

					P_(list_points_south_(i)) = P_gas;
					T_(list_points_south_(i)) = T_gas;
				}

				for (unsigned int i = 0; i < ny_; i++)
				{
					for (unsigned int j = 0; j < nc_; j++)
						Y_[list_points_east_(i)](j) = omega_gas[nx_ + i](j);

					P_(list_points_east_(i)) = P_gas;
					T_(list_points_east_(i)) = T_gas;
				}

				for (unsigned int i = 0; i < nx_; i++)
				{
					for (unsigned int j = 0; j < nc_; j++)
						Y_[list_points_north_(i)](j) = omega_gas[2 * nx_ + ny_ - 1 - i](j);

					P_(list_points_north_(i)) = P_gas;
					T_(list_points_north_(i)) = T_gas;
				}
			}
		}
	}

	void Reactor2D::SetGasSide(const double T_gas, const double P_gas, const CVI::DiskFromCFD& disk_from_cfd)
	{
		time_profiles_ = false;

		gaseous_phase_ = GASEOUS_PHASE_FROM_CFD;

		// Set initial fields consistent along the west side
		if (hole_ == false)
		{
			for (unsigned int i = 0; i < ny_; i++)
			{
				for (unsigned int j = 0; j < nc_; j++)
					Y_[list_points_west_(i)](j) = disk_from_cfd.east_mass_fractions()[i][j];
				P_(list_points_west_(i)) = P_gas;
				T_(list_points_west_(i)) = disk_from_cfd.east_temperature()[i];

				for (unsigned int j = 0; j < nc_; j++)
					Y_gas_west_side_[i](j) = disk_from_cfd.east_mass_fractions()[i][j];
				P_gas_west_side_(i) = P_gas;
				T_gas_west_side_(i) = disk_from_cfd.east_temperature()[i];
			}
		}
		else
		{
			for (unsigned int i = 0; i < ny_; i++)
			{
				for (unsigned int j = 0; j < nc_; j++)
					Y_[list_points_west_(i)](j) = disk_from_cfd.west_mass_fractions()[i][j];
				P_(list_points_west_(i)) = P_gas;
				T_(list_points_west_(i)) = disk_from_cfd.west_temperature()[i];

				for (unsigned int j = 0; j < nc_; j++)
					Y_gas_west_side_[i](j) = disk_from_cfd.west_mass_fractions()[i][j];
				P_gas_west_side_(i) = P_gas;
				T_gas_west_side_(i) = disk_from_cfd.west_temperature()[i];
			}
		}

		// Set initial fields consistent along the east side
		{
			for (unsigned int i = 0; i < ny_; i++)
			{
				for (unsigned int j = 0; j < nc_; j++)
					Y_[list_points_east_(i)](j) = disk_from_cfd.east_mass_fractions()[i][j];
				P_(list_points_east_(i)) = P_gas;
				T_(list_points_east_(i)) = disk_from_cfd.east_temperature()[i];

				for (unsigned int j = 0; j < nc_; j++)
					Y_gas_east_side_[i](j) = disk_from_cfd.east_mass_fractions()[i][j];
				P_gas_east_side_(i) = P_gas;
				T_gas_east_side_(i) = disk_from_cfd.east_temperature()[i];
			}
		}

		// Set initial fields consistent along the south side
		{
			for (unsigned int i = 0; i < nx_; i++)
			{
				for (unsigned int j = 0; j < nc_; j++)
					Y_[list_points_south_(i)](j) = disk_from_cfd.south_mass_fractions()[i][j];
				P_(list_points_south_(i)) = P_gas;
				T_(list_points_south_(i)) = disk_from_cfd.south_temperature()[i];

				for (unsigned int j = 0; j < nc_; j++)
					Y_gas_south_side_[i](j) = disk_from_cfd.south_mass_fractions()[i][j];
				P_gas_south_side_(i) = P_gas;
				T_gas_south_side_(i) = disk_from_cfd.south_temperature()[i];
			}
		}

		// Set initial fields consistent along the north side
		{
			for (unsigned int i = 0; i < nx_; i++)
			{
				for (unsigned int j = 0; j < nc_; j++)
					Y_[list_points_north_(i)](j) = disk_from_cfd.north_mass_fractions()[i][j];
				P_(list_points_north_(i)) = P_gas;
				T_(list_points_north_(i)) = disk_from_cfd.north_temperature()[i];

				for (unsigned int j = 0; j < nc_; j++)
					Y_gas_north_side_[i](j) = disk_from_cfd.north_mass_fractions()[i][j];
				P_gas_north_side_(i) = P_gas;
				T_gas_north_side_(i) = disk_from_cfd.north_temperature()[i];
			}
		}
	}

	void Reactor2D::SetGasSide(OpenSMOKE::FixedProfile* profile_temperature, const double P_gas, const CVI::DiskFromCFD& disk_from_cfd)
	{
		gaseous_phase_ = GASEOUS_PHASE_FROM_CFD;

		time_profiles_ = true;

		profile_temperature_ = new OpenSMOKE::FixedProfile(	profile_temperature->x().size(),
															profile_temperature->x().data(),
															profile_temperature->y().data());

		// Update boundary conditions
		UpdateTemperatureBoundaryConditions(0.);

		// Set initial fields consistent along the west side
		if (hole_ == false)
		{
			for (unsigned int i = 0; i < ny_; i++)
			{
				for (unsigned int j = 0; j < nc_; j++)
					Y_[list_points_west_(i)](j) = disk_from_cfd.east_mass_fractions()[i][j];
				P_(list_points_west_(i)) = P_gas;
				T_(list_points_west_(i)) = T_gas_east_side_(i);

				for (unsigned int j = 0; j < nc_; j++)
					Y_gas_west_side_[i](j) = disk_from_cfd.east_mass_fractions()[i][j];
				P_gas_west_side_(i) = P_gas;
			}
		}
		else
		{
			for (unsigned int i = 0; i < ny_; i++)
			{
				for (unsigned int j = 0; j < nc_; j++)
					Y_[list_points_west_(i)](j) = disk_from_cfd.west_mass_fractions()[i][j];
				P_(list_points_west_(i)) = P_gas;
				T_(list_points_west_(i)) = T_gas_west_side_(i);

				for (unsigned int j = 0; j < nc_; j++)
					Y_gas_west_side_[i](j) = disk_from_cfd.west_mass_fractions()[i][j];
				P_gas_west_side_(i) = P_gas;
			}
		}

		// Set initial fields consistent along the east side
		{
			for (unsigned int i = 0; i < ny_; i++)
			{
				for (unsigned int j = 0; j < nc_; j++)
					Y_[list_points_east_(i)](j) = disk_from_cfd.east_mass_fractions()[i][j];
				P_(list_points_east_(i)) = P_gas;
				T_(list_points_east_(i)) = T_gas_east_side_(i);

				for (unsigned int j = 0; j < nc_; j++)
					Y_gas_east_side_[i](j) = disk_from_cfd.east_mass_fractions()[i][j];
				P_gas_east_side_(i) = P_gas;
			}
		}

		// Set initial fields consistent along the south side
		{
			for (unsigned int i = 0; i < nx_; i++)
			{
				for (unsigned int j = 0; j < nc_; j++)
					Y_[list_points_south_(i)](j) = disk_from_cfd.south_mass_fractions()[i][j];
				P_(list_points_south_(i)) = P_gas;
				T_(list_points_south_(i)) = T_gas_south_side_(i);

				for (unsigned int j = 0; j < nc_; j++)
					Y_gas_south_side_[i](j) = disk_from_cfd.south_mass_fractions()[i][j];
				P_gas_south_side_(i) = P_gas;
			}
		}

		// Set initial fields consistent along the north side
		{
			for (unsigned int i = 0; i < nx_; i++)
			{
				for (unsigned int j = 0; j < nc_; j++)
					Y_[list_points_north_(i)](j) = disk_from_cfd.north_mass_fractions()[i][j];
				P_(list_points_north_(i)) = P_gas;
				T_(list_points_north_(i)) = T_gas_north_side_(i);

				for (unsigned int j = 0; j < nc_; j++)
					Y_gas_north_side_[i](j) = disk_from_cfd.north_mass_fractions()[i][j];
				P_gas_north_side_(i) = P_gas;
			}
		}
	}

	void Reactor2D::UpdateTemperatureBoundaryConditions(const double time)
	{
		T_gas_west_side_.setConstant(profile_temperature_->Interpolate(time));
		T_gas_east_side_.setConstant(profile_temperature_->Interpolate(time));
		T_gas_north_side_.setConstant(profile_temperature_->Interpolate(time));
		T_gas_south_side_.setConstant(profile_temperature_->Interpolate(time));

		// Update mass fractions of species (in case of PFR) TODO
		//Y_gas_side_.resize(nc_);
		//for (unsigned int j = 0; j < nc_; j++)
		//	Y_gas_side_(j) = profiles_omega_[j]->Interpolate(time);
	}

	void Reactor2D::UpdateTemperatureField(const double time)
	{
		T_.setConstant(profile_temperature_->Interpolate(time));
	}

	void Reactor2D::SetInitialConditions(const boost::filesystem::path path_to_backup_file, const double T_gas, const double P_gas, const Eigen::VectorXd& omega_gas, const Eigen::VectorXd& Gamma0, const Eigen::VectorXd& Z0)
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

		epsilon_.setConstant(porousMedium_.porosity());
		P_.setConstant(P_gas);
		//if (gaseous_phase_ != GASEOUS_PHASE_FROM_CFD)
			T_.setConstant(T_gas);

		// Set porosity defect
		if (porosityDefect_.is_active() == true)
		{
			for (unsigned int k = 1; k < ny_ - 1; k++)
				for (unsigned int i = 1; i < nx_ - 1; i++)
				{
					const int center = k*nx_ + i;
					epsilon_(center) = porosityDefect_.set_porosity(grid_x_.x()(i), grid_y_.x()(k), epsilon_(center));
				}
		}

		if (path_to_backup_file.empty() == false)
		{
			SetInitialConditionsFromBackupFile(path_to_backup_file);
			start_from_backup_ = true;
		}
	}

	void Reactor2D::SetTimeTotal(const double time_total)
	{
		time_total_ = time_total;
	}

	void Reactor2D::SetDaeTimeInterval(const double time_interval)
	{
		dae_time_interval_ = time_interval;
	}

	void Reactor2D::SetOdeEndTime(const double time_interval)
	{
		ode_end_time_ = time_interval;
	}

	void Reactor2D::SetTecplotTimeInterval(const double time_interval)
	{
		tecplot_time_interval_ = time_interval;
	}

	void Reactor2D::SetStepsVideo(const int steps_video)
	{
		n_steps_video_ = steps_video;
		count_dae_video_ = 1;
		count_ode_video_ = 1;
	}

	void Reactor2D::SetStepsFile(const int steps_file)
	{
		n_steps_file_ = steps_file;
		count_file_ = n_steps_file_;
	}

	void Reactor2D::SetStepsUpdatePlugFlow(const int steps_update_plug_flow)
	{
		n_steps_update_plug_flow_ = steps_update_plug_flow;
		count_update_plug_flow_ = n_steps_update_plug_flow_;
	}

	void Reactor2D::SetUniformVelocity(const double vx, const double vy)
	{
		vx_ = vx;
		vy_ = vy;
	}

	void Reactor2D::SetOutputFile(const std::string filename)
	{
		output_disk_file_name_ = filename;
	}

	void Reactor2D::SetDerivativeMassFractions(const OpenSMOKE::derivative_type value)
	{
		derivative_type_mass_fractions_ = value;
	}

	void Reactor2D::SetDerivativeEffectiveDiffusivity(const OpenSMOKE::derivative_type value)
	{
		derivative_type_effective_diffusivity_ = value;
	}

	void Reactor2D::SetDerivativeBulkDensity(const OpenSMOKE::derivative_type value)
	{
		derivative_type_bulk_density_ = value;
	}

	void Reactor2D::Properties()
	{
		for (unsigned int i = 0; i < np_; i++)
		{
			// Thermodynamics
			{
				thermodynamicsMap_.SetPressure(P_(i));
				thermodynamicsMap_.SetTemperature(T_(i));

				aux_Y.CopyFrom(Y_[i].data());
				thermodynamicsMap_.MoleFractions_From_MassFractions(aux_X.GetHandle(), mw_(i), aux_Y.GetHandle());
				aux_X.CopyTo(X_[i].data());

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
					porousMedium_.EffectiveDiffusionCoefficients(X_[i]);
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
					kineticsMap_.ReactionRates(aux_C.GetHandle());
					kineticsMap_.FormationRates(aux_R.GetHandle());
					OpenSMOKE::ElementByElementProduct(aux_R.Size(), aux_R.GetHandle(), thermodynamicsMap_.MWs().data(), aux_R.GetHandle()); // [kg/m3/s]
					aux_R.CopyTo(omega_homogeneous_from_homogeneous_[i].data());
				}

				// Heterogeneous phase
				if (heterogeneousMechanism_.heterogeneous_reactions() == true)
				{
					heterogeneousMechanism_.SetTemperature(T_(i));
					heterogeneousMechanism_.SetPressure(P_(i));

					aux_C.CopyTo(aux_eigen.data());
					heterogeneousMechanism_.FormationRates(Sv_(i), aux_eigen);
					for (unsigned int j = 0; j < nc_; j++)
						omega_homogeneous_from_heterogeneous_[i](j) = heterogeneousMechanism_.Rgas()(j)*thermodynamicsMap_.MW(j);								// [kg/m3/s]
					
					omega_deposition_per_unit_area_(i) = heterogeneousMechanism_.r_deposition_per_unit_area()*heterogeneousMechanism_.mw_carbon();	// [kg/m2/s]
					omega_deposition_per_unit_volume_(i) = heterogeneousMechanism_.r_deposition_per_unit_volume()*heterogeneousMechanism_.mw_carbon();		// [kg/m3/s]

					omega_loss_per_unit_volume_(i) = 0.;
					for (unsigned int j = 0; j < nc_; j++)
						omega_loss_per_unit_volume_(i) += heterogeneousMechanism_.Rgas()(j)*thermodynamicsMap_.MW(j);					// [kg/m3/s]
				}
			}
			else
			{
				const double coefficient = porousMedium_.epsilon_smoothing_coefficient();
				const double smoothing_coefficient = 0.50*(std::tanh(coefficient*(epsilon_(i) - porousMedium_.epsilon_threshold())) + 1.);

				heterogeneousDetailedMechanism_.SetTemperature(T_(i));
				heterogeneousDetailedMechanism_.SetPressure(P_(i));

				// Homogeneous phase
				if (heterogeneousDetailedMechanism_.homogeneous_reactions() == true)
				{
					kineticsMap_.SetTemperature(T_(i));
					kineticsMap_.SetPressure(P_(i));
					kineticsMap_.ReactionRates(aux_C.GetHandle());
					kineticsMap_.FormationRates(aux_R.GetHandle());
					OpenSMOKE::ElementByElementProduct(aux_R.Size(), aux_R.GetHandle(), thermodynamicsMap_.MWs().data(), aux_R.GetHandle()); // [kg/m3/s]
					aux_R.CopyTo(omega_homogeneous_from_homogeneous_[i].data());

					// Smoothing
					omega_homogeneous_from_homogeneous_[i] *= smoothing_coefficient;
				}

				// Heterogeneous phase (detailed)
				if (heterogeneousDetailedMechanism_.heterogeneous_reactions() == true)
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
						omega_homogeneous_from_heterogeneous_[i](j) = heterogeneousDetailedMechanism_.Rgas()(j)*thermodynamicsMap_.MW(j);				// [kg/m3/s]

					for (unsigned int j = 0; j < surf_nc_; j++)
						omega_heterogeneous_from_heterogeneous_[i](j) = heterogeneousDetailedMechanism_.Rsurface()(j);	// [kmol/m2/s]

					omega_deposition_per_unit_area_(i) = heterogeneousDetailedMechanism_.r_deposition_per_unit_area()*heterogeneousDetailedMechanism_.mw_carbon();			// [kg/m2/s]
					omega_deposition_per_unit_volume_(i) = heterogeneousDetailedMechanism_.r_deposition_per_unit_volume()*heterogeneousDetailedMechanism_.mw_carbon();		// [kg/m3/s]

					omega_loss_per_unit_volume_(i) = 0.;
					for (unsigned int j = 0; j < nc_; j++)
						omega_loss_per_unit_volume_(i) += heterogeneousDetailedMechanism_.Rgas()(j)*thermodynamicsSurfaceMap_.MW(j);					// [kg/m3/s]

					// Smoothing
					omega_homogeneous_from_heterogeneous_[i] *= smoothing_coefficient;
					omega_heterogeneous_from_heterogeneous_[i] *= smoothing_coefficient;
					omega_deposition_per_unit_area_(i) *= smoothing_coefficient;
					omega_deposition_per_unit_volume_(i) *= smoothing_coefficient;
					omega_loss_per_unit_volume_(i) *= smoothing_coefficient;
				}
			}
		}
	}

	void Reactor2D::SubEquations_MassFractions_BoundaryConditions_WestSide(const double t)
	{
		if (gaseous_phase_ == GASEOUS_PHASE_FROM_PLUG_FLOW)
		{
			for (unsigned int i = 0; i < ny_; i++)
			{
				const int point = list_points_west_(i);
				for (unsigned int j = 0; j < nc_; j++)
					dY_over_dt_[point](j) = Y_[point](j) - Y_[point + 1](j);
			}
		}
		else if (gaseous_phase_ == GASEOUS_PHASE_FROM_CFD)
		{
			if (hole_ == false)
			{
				for (unsigned int i = 0; i < ny_; i++)
				{
					const int point = list_points_west_(i);
					for (unsigned int j = 0; j < nc_; j++)
						dY_over_dt_[point](j) = Y_[point](j) - Y_[point + 1](j);
				}
			}
			else
			{
				for (unsigned int i = 0; i < ny_; i++)
				{
					const int point = list_points_west_(i);
					for (unsigned int j = 0; j < nc_; j++)
						dY_over_dt_[point](j) = Y_[point](j) - Y_gas_west_side_[i](j);
				}
			}
		}
	}

	void Reactor2D::SubEquations_MassFractions_BoundaryConditions_EastSide(const double t)
	{
		if (gaseous_phase_ == GASEOUS_PHASE_FROM_PLUG_FLOW)
		{
			if (plugFlowReactor_.internal_boundary_layer_correction() == true && t > time_smoothing_)
			{
				for (unsigned int i = 0; i < ny_; i++)
				{
					const int point = list_points_east_(i);

					const double kc = plugFlowReactor_.mass_transfer_coefficient(T_(point), P_(point), rho_gas_(point), grid_y_.x()(i) + plugFlowReactor_.inert_length());
					const double gamma_over_dx = gamma_star_[point](heterogeneousMechanism_.index_CH4()) / grid_x_.dxw()(nx_ - 1);

					for (unsigned int j = 0; j < nc_; j++)
						dY_over_dt_[point](j) = Y_[point](j) - (gamma_over_dx*Y_[point - 1](j) + kc*Y_gas_east_side_[i](j)) / (gamma_over_dx + kc);
				}
			}
			else
			{
				for (unsigned int i = 0; i < ny_; i++)
				{
					const int point = list_points_east_(i);
					for (unsigned int j = 0; j < nc_; j++)
						dY_over_dt_[point](j) = Y_[point](j) - Y_gas_east_side_[i](j);
				}
			}
		}
		else if (gaseous_phase_ == GASEOUS_PHASE_FROM_CFD)
		{
			for (unsigned int i = 0; i < ny_; i++)
			{
				const int point = list_points_east_(i);
				for (unsigned int j = 0; j < nc_; j++)
					dY_over_dt_[point](j) = Y_[point](j) - Y_gas_east_side_[i](j);
			}
		}
	}

	void Reactor2D::SubEquations_MassFractions_BoundaryConditions_NorthSide(const double t)
	{
		if (gaseous_phase_ == GASEOUS_PHASE_FROM_PLUG_FLOW)
		{
			if (plugFlowReactor_.geometric_pattern() == CVI::PlugFlowReactorCoupled::ONE_SIDE)
			{
				for (unsigned int i = 0; i < nx_; i++)
				{
					const int point = list_points_north_(i);
					for (unsigned int j = 0; j < nc_; j++)
						dY_over_dt_[point](j) = Y_[point](j) - Y_[point - nx_](j);
				}
			}
			else if (plugFlowReactor_.geometric_pattern() == CVI::PlugFlowReactorCoupled::THREE_SIDES)
			{
				for (unsigned int i = 0; i < nx_; i++)
				{
					const int point = list_points_north_(i);
					for (unsigned int j = 0; j < nc_; j++)
						dY_over_dt_[point](j) = Y_[point](j) - Y_gas_north_side_[i](j);
				}
			}
		}
		else if (gaseous_phase_ == GASEOUS_PHASE_FROM_CFD)
		{
			for (unsigned int i = 0; i < nx_; i++)
			{
				const int point = list_points_north_(i);
				for (unsigned int j = 0; j < nc_; j++)
					dY_over_dt_[point](j) = Y_[point](j) - Y_gas_north_side_[i](j);
			}
		}
	}

	void Reactor2D::SubEquations_MassFractions_BoundaryConditions_SouthSide(const double t)
	{
		if (gaseous_phase_ == GASEOUS_PHASE_FROM_PLUG_FLOW)
		{
			if (plugFlowReactor_.geometric_pattern() == CVI::PlugFlowReactorCoupled::ONE_SIDE)
			{
				for (unsigned int i = 0; i < nx_; i++)
				{
					const int point = list_points_south_(i);
					for (unsigned int j = 0; j < nc_; j++)
						dY_over_dt_[point](j) = Y_[point](j) - Y_[point + nx_](j);
				}
			}
			else if (plugFlowReactor_.geometric_pattern() == CVI::PlugFlowReactorCoupled::THREE_SIDES)
			{
				for (unsigned int i = 0; i < nx_; i++)
				{
					const int point = list_points_south_(i);
					for (unsigned int j = 0; j < nc_; j++)
						dY_over_dt_[point](j) = Y_[point](j) - Y_gas_south_side_[i](j);
				}
			}
		}
		else if (gaseous_phase_ == GASEOUS_PHASE_FROM_CFD)
		{
			for (unsigned int i = 0; i < nx_; i++)
			{
				const int point = list_points_south_(i);
				for (unsigned int j = 0; j < nc_; j++)
					dY_over_dt_[point](j) = Y_[point](j) - Y_gas_south_side_[i](j);
			}
		}
	}

	void Reactor2D::SubEquations_MassFractions_BoundaryConditions(const double t)
	{
		// South (zero gradient)
		SubEquations_MassFractions_BoundaryConditions_SouthSide(t);

		// North (zero gradient)
		SubEquations_MassFractions_BoundaryConditions_NorthSide(t);

		// West (zero gradient)
		SubEquations_MassFractions_BoundaryConditions_WestSide(t);

		// East (gas side)
		SubEquations_MassFractions_BoundaryConditions_EastSide(t);
	}

	void Reactor2D::SubEquations_MassFractions(const double t)
	{
		// Boundary conditions (according to the geometric pattern)
		SubEquations_MassFractions_BoundaryConditions(t);

		// Internal
		for (unsigned int k = 1; k < ny_ - 1; k++)
			for (unsigned int i = 1; i < nx_ - 1; i++)
			{
				const int center = k*nx_ + i;
				const int east  = center + 1;
				const int west  = center - 1;
				const int north = center + nx_;
				const int south = center - nx_;

				for (unsigned int j = 0; j < nc_; j++)
				{
					const double c_east = 0.50* (gamma_star_[east](j)*rho_gas_[east] + gamma_star_[center](j)*rho_gas_[center]);
					const double c_west = 0.50* (gamma_star_[west](j)*rho_gas_[west] + gamma_star_[center](j)*rho_gas_[center]);
					const double c_north = 0.50* (gamma_star_[north](j)*rho_gas_[north] + gamma_star_[center](j)*rho_gas_[center]);
					const double c_south = 0.50* (gamma_star_[south](j)*rho_gas_[south] + gamma_star_[center](j)*rho_gas_[center]);

					      double diffusion_x = (c_east*(Y_[east](j) - Y_[center](j)) / grid_x_.dxe()(i)-c_west*(Y_[center](j) - Y_[west](j)) / grid_x_.dxw()(i)) / grid_x_.dxc_over_2()(i);
					const double diffusion_y = (c_north*(Y_[north](j) - Y_[center](j)) / grid_y_.dxe()(k)-c_south*(Y_[center](j) - Y_[south](j)) / grid_y_.dxw()(k)) / grid_y_.dxc_over_2()(k);
					
					      double convection_x = rho_gas_[center] * vx_ * (Y_[center](j) - Y_[west](j)) / grid_x_.dxw()(i);
					const double convection_y = rho_gas_[center] * vy_ * (Y_[center](j) - Y_[south](j)) / grid_y_.dxw()(k);
					
					const double homogeneous_reactions = epsilon_(center)*omega_homogeneous_from_homogeneous_[center](j);
					const double heterogeneous_reactions = omega_homogeneous_from_heterogeneous_[center](j) + Y_[center](j)*omega_deposition_per_unit_volume_(center);



					if (planar_symmetry_ == false)
					{
						const double dY_over_dr = (Y_[east](j) - Y_[west](j)) / (grid_x_.x()[i + 1] - grid_x_.x()[i - 1]);
						const double diffusion_radial = gamma_star_[center](j)*rho_gas_[center] * dY_over_dr / grid_x_.x()[i];
						diffusion_x += diffusion_radial;

						// Correction because of the cylindrical symmetry
						const double vr = vx_*grid_x_.x()[0]/grid_x_.x()[i];
						convection_x = rho_gas_[center] * vr * (Y_[center](j) - Y_[west](j)) / grid_x_.dxw()(i);
					}

					const double diffusion = diffusion_x + diffusion_y;
					const double convection = convection_x + convection_y;

					dY_over_dt_[center](j) = diffusion + homogeneous_reactions + heterogeneous_reactions - convection;
					dY_over_dt_[center](j) /= (rho_gas_(center)*epsilon_(center));
				}

				dY_over_dt_[center](gas_dae_species_index_) = 1.-Y_[center].sum();
			}
	}

	void Reactor2D::SubEquations_Porosity()
	{
		for (unsigned int i = 0; i < np_; i++)
			depsilon_over_dt_(i) = -omega_deposition_per_unit_volume_(i) / rho_graphite_;
	}

	void Reactor2D::SubEquations_SurfaceSpeciesFractions()
	{
		for (unsigned int i = 0; i < np_; i++)
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
				dZ_over_dt_[i](surface_dae_species_index_) = 1. - Z_[i].sum();
		}
	}

	void Reactor2D::SubEquations_Temperature_BoundaryConditions(const double t)
	{
		// South (Dirichlet condition)
		SubEquations_Temperature_BoundaryConditions_SouthSide(t);

		// North (Dirichlet condition)
		SubEquations_Temperature_BoundaryConditions_NorthSide(t);

		// West (Dirichlet condition)
		SubEquations_Temperature_BoundaryConditions_WestSide(t);

		// East (Dirichlet condition)
		SubEquations_Temperature_BoundaryConditions_EastSide(t);
	}

	void Reactor2D::SubEquations_Temperature_BoundaryConditions_WestSide(const double t)
	{
		if (gaseous_phase_ == GASEOUS_PHASE_FROM_CFD)
		{
			if (hole_ == false)
			{
				for (unsigned int i = 0; i < ny_; i++)
				{
					const int point = list_points_west_(i);
					dT_over_dt_(point) = T_(point) - T_(point+1);
				}
			}
			else
			{
				for (unsigned int i = 0; i < ny_; i++)
				{
					const int point = list_points_west_(i);
					dT_over_dt_(point) = T_(point) - T_gas_west_side_(i);
				}
			}
		}
	}

	void Reactor2D::SubEquations_Temperature_BoundaryConditions_EastSide(const double t)
	{
		if (gaseous_phase_ == GASEOUS_PHASE_FROM_CFD)
		{
			for (unsigned int i = 0; i < ny_; i++)
			{
				const int point = list_points_east_(i);
				dT_over_dt_(point) = T_(point) - T_gas_east_side_(i);
			}
		}
	}

	void Reactor2D::SubEquations_Temperature_BoundaryConditions_NorthSide(const double t)
	{
		if (gaseous_phase_ == GASEOUS_PHASE_FROM_CFD)
		{
			for (unsigned int i = 0; i < nx_; i++)
			{
				const int point = list_points_north_(i);
				dT_over_dt_(point) = T_(point) - T_gas_north_side_(i);
			}
		}
	}

	void Reactor2D::SubEquations_Temperature_BoundaryConditions_SouthSide(const double t)
	{
		if (gaseous_phase_ == GASEOUS_PHASE_FROM_CFD)
		{
			for (unsigned int i = 0; i < nx_; i++)
			{
				const int point = list_points_south_(i);
				dT_over_dt_(point) = T_(point) - T_gas_south_side_(i);
			}
		}
	}

	void Reactor2D::SubEquations_Temperature(const double t)
	{
		// Properties
		Eigen::VectorXd rho_times_cp(np_);
		Eigen::VectorXd lambda(np_);

		for (unsigned int i = 0; i < np_; i++)
		{
			porousMedium_.SetPorosity(epsilon_(i));
			porousMedium_.SetTemperature(T_(i));
			porousMedium_.SetPressure(P_(i));

			// Properties
			rho_times_cp(i) = porousMedium_.cp_times_rho(rho_gas_(i), rho_graphite_);
			lambda(i) = porousMedium_.lambda();
		}

		// Boundary conditions (according to the geometric pattern)
		SubEquations_Temperature_BoundaryConditions(t);

		// Internal
		for (unsigned int k = 1; k < ny_ - 1; k++)
			for (unsigned int i = 1; i < nx_ - 1; i++)
			{
				const int center = k*nx_ + i;
				const int east   = center + 1;
				const int west   = center - 1;
				const int north  = center + nx_;
				const int south  = center - nx_;

				{
					const double kappa_east = 0.50* (lambda(east) + lambda(center));
					const double kappa_west = 0.50* (lambda(west) + lambda(center));
					const double kappa_north = 0.50* (lambda(north) + lambda(center));
					const double kappa_south = 0.50* (lambda(south) + lambda(center));

					double conduction_x = (kappa_east*(T_(east) - T_(center)) / grid_x_.dxe()(i) - kappa_west*(T_(center) - T_(west)) / grid_x_.dxw()(i)) / grid_x_.dxc_over_2()(i);
					const double conduction_y = (kappa_north*(T_(north) - T_(center)) / grid_y_.dxe()(k) - kappa_south*(T_(center) - T_(south)) / grid_y_.dxw()(k)) / grid_y_.dxc_over_2()(k);

					if (planar_symmetry_ == false)
					{
						const double dT_over_dr = (T_(east) - T_(west)) / (grid_x_.x()[i + 1] - grid_x_.x()[i - 1]);
						const double conduction_radial = lambda(center) * dT_over_dr / grid_x_.x()[i];
						conduction_x += conduction_radial;
					}

					const double conduction = conduction_x + conduction_y;

					dT_over_dt_(center) = conduction;
					dT_over_dt_(center) /= rho_times_cp(center);
				}
			}
	}

	int Reactor2D::OdeEquations(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y, OpenSMOKE::OpenSMOKEVectorDouble& dy)
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
		thermodynamicsMap_.MoleFractions_From_MassFractions(aux_X.GetHandle(), mw_(i_current), aux_Y.GetHandle());

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

			kineticsSurfaceMap_.ReactionRates(aux_C.GetHandle(), Z.GetHandle(), aux_a.GetHandle(), Gamma.GetHandle());
			kineticsSurfaceMap_.FormationRates(RfromSurface.GetHandle(), Rsurface.GetHandle(), Rbulk.GetHandle(), RsurfacePhases.GetHandle());
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

	int Reactor2D::OdePrint(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y)
	{
		// Video output
		const bool verbose_ode = true;
		if (count_ode_video_%n_steps_video_ == 1 && verbose_ode == true)
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

	void Reactor2D::Recover_Unknowns(const double* y)
	{
		if (equations_set_ == EQUATIONS_SET_COMPLETE)
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
		else if (equations_set_ == EQUATIONS_SET_ONLYTEMPERATURE)
		{
			for (unsigned int i = 0; i < np_; i++)
				T_(i) = y[i];
		}
	}

	void Reactor2D::Recover_Residuals(double* dy)
	{
		if (equations_set_ == EQUATIONS_SET_COMPLETE)
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
		else if (equations_set_ == EQUATIONS_SET_ONLYTEMPERATURE)
		{
			for (unsigned int i = 0; i < np_; i++)
				dy[i] = dT_over_dt_(i);
		}
	}

	void Reactor2D::AlgebraicDifferentialVector(double* v)
	{
		int count = 0;
		for (unsigned int i = 0; i < ne_; i++)
		{
			if (id_equations_[i] == true)  v[count++] = 1.;
			if (id_equations_[i] == false) v[count++] = 0.;
		}
	}

	void Reactor2D::UnknownsVector(double* v)
	{
		if (equations_set_ == EQUATIONS_SET_COMPLETE)
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
		else if (equations_set_ == EQUATIONS_SET_ONLYTEMPERATURE)
		{
			for (unsigned int i = 0; i < np_; i++)
					v[i] = T_(i);
		}
	}

	void Reactor2D::CorrectedUnknownsVector(double* v)
	{
		if (equations_set_ == EQUATIONS_SET_COMPLETE)
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

				// Normalize
				//const double sum = Y_[i].sum();
				//for (unsigned int j = 0; j < nc_; j++)
				//	Y_[i](j) /= sum;
			}
		}
		else if (equations_set_ == EQUATIONS_SET_ONLYTEMPERATURE)
		{
			for (unsigned int i = 0; i < np_; i++)
				T_(i) = v[i];
		}

		Properties();
	}

	void Reactor2D::MinimumUnknownsVector(double* v)
	{
		const double zero = 0.;
		const double minimum_temperature = 200.;

		if (equations_set_ == EQUATIONS_SET_COMPLETE)
		{
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
		else if (equations_set_ == EQUATIONS_SET_ONLYTEMPERATURE)
		{
			for (unsigned int i = 0; i < np_; i++)
				v[i] = minimum_temperature;
		}
	}

	void Reactor2D::MaximumUnknownsVector(double* v)
	{
		const double one = 1.;
		const double maximum_temperature = 2500.;

		if (equations_set_ == EQUATIONS_SET_COMPLETE)
		{
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
		else if (equations_set_ == EQUATIONS_SET_ONLYTEMPERATURE)
		{
			for (unsigned int i = 0; i < np_; i++)
				v[i] = maximum_temperature;
		}
	}

	void Reactor2D::Equations(const double t, const double* y, double* dy)
	{
		if (equations_set_ == EQUATIONS_SET_COMPLETE)
			EquationsComplete(t, y, dy);
		else if (equations_set_ == EQUATIONS_SET_ONLYTEMPERATURE)
			EquationsOnlyTemperature(t, y, dy);
	}

	void Reactor2D::EquationsComplete(const double t, const double* y, double* dy)
	{
		// Recover unknowns
		Recover_Unknowns(y);

		// Update boundary conditions
		if (time_profiles_ == true)
		{
			UpdateTemperatureBoundaryConditions(t);
			UpdateTemperatureField(t);
		}

		// Properties
		Properties();

		// Equations
		SubEquations_MassFractions(t);
		SubEquations_Porosity();
		SubEquations_SurfaceSpeciesFractions();

		// Recover residuals
		Recover_Residuals(dy);
	}

	void Reactor2D::EquationsOnlyTemperature(const double t, const double* y, double* dy)
	{
		// Recover unknowns
		Recover_Unknowns(y);

		// Update boundary conditions
		if (time_profiles_ == true)
			std::cout << "Solving equations of temperature only" << std::endl;

		// Properties
		Properties();

		// Equations
		SubEquations_Temperature(t);

		// Recover residuals
		Recover_Residuals(dy);
	}

	int Reactor2D::SolveFromScratch(DaeSMOKE::DaeSolver_Parameters& dae_parameters, OdeSMOKE::OdeSolver_Parameters& ode_parameters)
	{
		// If the boundary conditions are taken from the CFD,
		// the temperature is solved
		if (gaseous_phase_ == GASEOUS_PHASE_FROM_CFD)
		{
			equations_set_ = EQUATIONS_SET_ONLYTEMPERATURE;
			ne_ = np_;
			unsigned int previous_block = block_;
			block_ = 1;
			band_size_ = (nx_ + 1) - 1;
			SetAlgebraicDifferentialEquations();

			std::cout << std::endl;
			std::cout << "--------------------------------------------------------------------------" << std::endl;
			std::cout << "Solving the temperature equation..." << std::endl;
			std::cout << "--------------------------------------------------------------------------" << std::endl;

			// Solve
			Properties();
			int flag = Solve(dae_parameters, 0, 10000.);
			
			equations_set_ = EQUATIONS_SET_COMPLETE;
			block_ = previous_block;
			ne_ = block_*np_;
			band_size_ = block_*(nx_ + 1) - 1;
			count_file_ = n_steps_file_;
			SetAlgebraicDifferentialEquations();
		}

		// Print initial solution
		Properties();
		PrintTecplot(0., (output_tecplot_folder_ / "Solution.tec.0").string().c_str());
		
		if (detailed_heterogeneous_kinetics_ == true && start_from_backup_ == true)
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
					typedef OdeSMOKE::KernelDense<OpenSMOKE::ODESystem_OpenSMOKE_Reactor2D> denseOde;
					typedef OdeSMOKE::MethodGear<denseOde> methodGear;
					OdeSMOKE::MultiValueSolver<methodGear> ode_solver;
					ode_solver.SetReactor2D(this);

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
		unsigned int number_intervals = static_cast<unsigned int>(time_total_/dae_time_interval_);
		for (unsigned int k = 1; k <= number_intervals; k++)
		{
			const double t0 = time_starting_point_ + (k - 1)*dae_time_interval_;
			const double tf = t0 + dae_time_interval_;

			// Reset to zero
			t_old_ = t0;
			homogeneous_total_mass_source_.setZero();
			heterogeneous_total_mass_source_.setZero();
			total_mass_exchanged_.setZero();

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
			std::string integral_homogeneous_rates_file = output_disk_file_name_ + ".source." + hours.str() + ".xml";

			PrintSolution(tf, (output_folder_ / solution_file).string().c_str());
			PrintDiffusionCoefficients(tf, (output_diffusion_folder_ / diffusion_coefficients_file).string().c_str());
			PrintHeterogeneousRates(tf, (output_heterogeneous_folder_ / heterogeneous_rates_file).string().c_str());
			PrintHomogeneousRates(tf, (output_homogeneous_folder_ / homogeneous_rates_file).string().c_str());
			PrintIntegralHomogeneousRates(tf, (output_disks_source_terms_folder_ / integral_homogeneous_rates_file).string().c_str());

			if (ropa_analysis_ == true)
			{
				std::cout << "Writing ROPA file..." << std::endl;
				std::string ropa_file = "ROPA." + hours.str() + ".out";
				PrintROPA(tf, (output_ropa_folder_ / ropa_file).string().c_str());
			}

			if (detailed_heterogeneous_kinetics_ == true)
			{
				std::cout << "Writing Backup file..." << std::endl;
				std::string backup_file = "Solution." + hours.str() + ".xml";
				PrintXMLFile( (output_backup_folder_ / backup_file).string(), tf);
			}
		}

		return true;
	}

	int Reactor2D::Solve(DaeSMOKE::DaeSolver_Parameters& dae_parameters, const double t0, const double tf)
	{
		int flag = 0;
		
		if (dae_parameters.type() == DaeSMOKE::DaeSolver_Parameters::DAE_INTEGRATOR_OPENSMOKEPP)
		{
			if (dae_parameters.sparse_linear_algebra() == false)
				flag = DaeSMOKE::Solve_Band_OpenSMOKEppDae<Reactor2D, OpenSMOKE_Reactor2D_DaeSystem>(this, dae_parameters, t0, tf);
			else
				flag = DaeSMOKE::Solve_Sparse_OpenSMOKEppDae<Reactor2D, OpenSMOKE_Reactor2D_DaeSystem>(this, dae_parameters, t0, tf);
		}
		#if OPENSMOKE_USE_BZZMATH == 1
		else if (dae_parameters.type() == DaeSMOKE::DaeSolver_Parameters::DAE_INTEGRATOR_BZZDAE)
		{
			if(dae_parameters.jacobian_structure() == OpenSMOKE::JACOBIAN_STRUCTURE_BAND)
				flag = DaeSMOKE::Solve_Band_BzzDae<Reactor2D, OpenSMOKE_Reactor2D_BzzDaeSystem>(this, dae_object_, dae_parameters, t0, tf);
			else if (dae_parameters.jacobian_structure() == OpenSMOKE::JACOBIAN_STRUCTURE_TRIDIAGONAL_BLOCK)
				flag = DaeSMOKE::Solve_TridiagonalBlock_BzzDae_Reactor2D<Reactor2D, OpenSMOKE_Reactor2D_BzzDaeSystem>(this, dae_object_, dae_parameters, t0, tf);
		}
		#endif
		#if OPENSMOKE_USE_SUNDIALS == 1
		else if(dae_parameters.type() == DaeSMOKE::DaeSolver_Parameters::DAE_INTEGRATOR_IDA)
			flag = DaeSMOKE::Solve_Band_Ida<Reactor2D>(this, dae_parameters, t0, tf);
		#endif
		#if OPENSMOKE_USE_DASPK == 1
		else if (dae_parameters.type() == DaeSMOKE::DaeSolver_Parameters::DAE_INTEGRATOR_DASPK)
			flag = DaeSMOKE::Solve_Band_Daspk<Reactor2D>(this, dae_parameters, t0, tf);
		#endif
		
		return flag;
	}

	void Reactor2D::PrintTecplot(const double t, const std::string name_file)
	{
		if (detailed_heterogeneous_kinetics_ == false)
			PrintTecplotGlobalKinetics(t, name_file);
		else
			PrintTecplotDetailedKinetics(t, name_file);
	}

	void Reactor2D::PrintTecplotGlobalKinetics(const double t, const std::string name_file)
	{
		OpenSMOKE::OpenSMOKEVectorDouble yy(nc_);
		OpenSMOKE::OpenSMOKEVectorDouble xx(nc_);

		std::ofstream fOutput(name_file.c_str(), std::ios::out);

		// Tecplot file
		{
			fOutput << "Title=Solution" << std::endl;
			fOutput << "Variables=";
			fOutput << "\"x[mm]\"" << ", ";
			fOutput << "\"y[mm]\"" << ", ";
			fOutput << "\"time[h]\"" << ", ";
			fOutput << "\"T[K]\"" << ", ";
			fOutput << "\"P[Pa]\"" << ", ";
			fOutput << "\"eps[-]\"" << ", ";
			fOutput << "\"rhoBulk[kg/m3]\"" << ", ";
			fOutput << "\"Sv[1/m]\"" << ", ";
			fOutput << "\"rpore[micron]\"" << ", ";
			fOutput << "\"K[m2]\"" << ", ";
			fOutput << "\"TauBulk[-]\"" << ", ";
			fOutput << "\"TauKnudsen[-]\"" << ", ";
			fOutput << "\"TauViscous[-]\"" << ", ";
			fOutput << "\"rhoGas[kg/m3]\"" << ", ";
			fOutput << "\"mwGas[kg/kmol]\"" << ", ";

			// Species (mass fractions)
			for (unsigned int j = 0; j < nc_; j++)
				fOutput << "\"" << thermodynamicsMap_.NamesOfSpecies()[j] << "\", ";
			
			// Heterogeneous reaction rates
			{
				fOutput << "\"rDep[kg/m2/s]\"" << ", ";
				fOutput << "\"rDep[kg/m3/s]\"" << ", ";
				fOutput << "\"rDep[m/s]\"" << ", ";
				for (int j = 0; j < heterogeneousMechanism_.r().size(); j++)
				{
					fOutput << "\"" << ("rDep" + heterogeneousMechanism_.tags()[j] + "[kg/m2/s]") << "\", ";
					fOutput << "\"" << ("rDep" + heterogeneousMechanism_.tags()[j] + "[m/s]") << "\", ";
				}
			}

			// Heterogeneous formation rates [kmol/m3/s] 
			{
				fOutput << "\"RhetCH4[kmml/m3/s]\"" << ", ";
				fOutput << "\"RhetC2H4[kml/m3/s]\"" << ", ";
				fOutput << "\"RhetC2H2[kml/m3/s]\"" << ", ";
				fOutput << "\"RhetC6H6[kml/m3/s]\"" << ", ";
				fOutput << "\"RhetH2[kml/m3/s]\"" << ", ";
			}
			
			// Hydrogen inhibition factors
			{
				fOutput << "\"InhCH4\"" << ", ";
				fOutput << "\"InhC2H4\"" << ", ";
				fOutput << "\"InhC2H2\"" << ", ";
				fOutput << "\"InhC6H6\"" << ", ";
			}
			
			// Heterogeneous reaction rates: contributions to the bulk density
			{
				fOutput << "\"dRhoBulk[kg/m3]\"" << ", ";
				for (int j = 0; j < heterogeneousMechanism_.r().size(); j++)
				{
					fOutput << "\"" << ("dRhoBulk" + heterogeneousMechanism_.tags()[j] + "[kg/m3]") << "\", ";
					fOutput << "\"" << ("dRhoBulk" + heterogeneousMechanism_.tags()[j]) << "\", ";
				}
			}
			
			// Finalize
			fOutput << std::endl;

			// Titles
			fOutput << "Zone I=" << nx_ << ", J=" << ny_ << ", F=POINT" << std::endl;
		}

		for (unsigned int k = 0; k < ny_; k++)
			for (unsigned int i = 0; i < nx_; i++)
			{
				const int point = k*nx_ + i;

				// Porous medium
				porousMedium_.SetPorosity(epsilon_(point));
				porousMedium_.SetTemperature(T_(point));
				porousMedium_.SetPressure(P_(point));

				// Thermodynamics
				thermodynamicsMap_.SetPressure(P_(point));
				thermodynamicsMap_.SetTemperature(T_(point));

				// Concentrations
				aux_Y.CopyFrom(Y_[point].data());
				thermodynamicsMap_.MoleFractions_From_MassFractions(aux_X.GetHandle(), mw_(point), aux_Y.GetHandle());
				aux_X.CopyTo(X_[point].data());
				const double cTot = P_(point) / PhysicalConstants::R_J_kmol / T_(point); // [kmol/m3]
				Product(cTot, aux_X, &aux_C);

				// Heterogeneous reactions
				heterogeneousMechanism_.SetTemperature(T_(point));
				heterogeneousMechanism_.SetPressure(P_(point));
				for (unsigned int j = 0; j < nc_; j++)
					aux_eigen(j) = aux_C[j + 1];
				heterogeneousMechanism_.FormationRates(porousMedium_.Sv(), aux_eigen);


				// Write grid
				fOutput << std::setprecision(9) << std::setw(20) << grid_x_.x()(i);
				fOutput << std::setprecision(9) << std::setw(20) << grid_y_.x()(k);

				// Write rlevant variables
				fOutput << std::setprecision(9) << std::setw(20) << t/3600.;
				fOutput << std::setprecision(9) << std::setw(20) << T_(point);
				fOutput << std::setprecision(9) << std::setw(20) << P_(point);

				// Write solid phase properties
				fOutput << std::setprecision(9) << std::setw(20) << epsilon_(point);
				fOutput << std::setprecision(9) << std::setw(20) << rho_bulk_(point);
				fOutput << std::setprecision(9) << std::setw(20) << Sv_(point);
				fOutput << std::setprecision(9) << std::setw(20) << rp_(point)*1.e6;
				fOutput << std::setprecision(9) << std::setw(20) << permeability_(point);
				fOutput << std::setprecision(9) << std::setw(20) << eta_bulk_(point);
				fOutput << std::setprecision(9) << std::setw(20) << eta_knudsen_(point);
				fOutput << std::setprecision(9) << std::setw(20) << eta_viscous_(point);

				// Write gas phase properties
				fOutput << std::setprecision(9) << std::setw(20) << rho_gas_(point);
				fOutput << std::setprecision(9) << std::setw(20) << mw_(point);

				// Write species mass fractions
				for (unsigned int j = 0; j < nc_; j++)
					fOutput << std::setprecision(9) << std::setw(20) << Y_[point](j);

				// Heterogeneous reaction rates
				{
					fOutput << std::setprecision(9) << std::setw(20) << heterogeneousMechanism_.r_deposition_per_unit_area()*heterogeneousMechanism_.mw_carbon();		// [kg/m2/s]
					fOutput << std::setprecision(9) << std::setw(20) << heterogeneousMechanism_.r_deposition_per_unit_volume()*heterogeneousMechanism_.mw_carbon();		// [kg/m3/s]
					fOutput << std::setprecision(9) << std::setw(20) << heterogeneousMechanism_.r_deposition_per_unit_volume()*heterogeneousMechanism_.mw_carbon() /
						(porousMedium_.Sv() / rho_graphite_);																		// [m/s]
					for (int j = 0; j < heterogeneousMechanism_.r().size(); j++)
					{
						fOutput << std::setprecision(9) << std::setw(20) << heterogeneousMechanism_.r()(j)*heterogeneousMechanism_.mw_carbon();									// [kg/m2/s]
						fOutput << std::setprecision(9) << std::setw(20) << heterogeneousMechanism_.r()(j)*heterogeneousMechanism_.mw_carbon() / rho_graphite_;	// [m/s]
					}
				}

				// Heterogeneous formation rates
				{
					fOutput << std::setprecision(9) << std::setw(20) << -heterogeneousMechanism_.Rgas()(heterogeneousMechanism_.index_CH4());
					fOutput << std::setprecision(9) << std::setw(20) << -heterogeneousMechanism_.Rgas()(heterogeneousMechanism_.index_C2H4());
					fOutput << std::setprecision(9) << std::setw(20) << -heterogeneousMechanism_.Rgas()(heterogeneousMechanism_.index_C2H2());
					fOutput << std::setprecision(9) << std::setw(20) << -heterogeneousMechanism_.Rgas()(heterogeneousMechanism_.index_C6H6());
					fOutput << std::setprecision(9) << std::setw(20) << -heterogeneousMechanism_.Rgas()(heterogeneousMechanism_.index_H2());
				}

				// Hydrogen inhibition factors
				{
					fOutput << std::setprecision(9) << std::setw(20) << heterogeneousMechanism_.I_CH4();
					fOutput << std::setprecision(9) << std::setw(20) << heterogeneousMechanism_.I_C2H4();
					fOutput << std::setprecision(9) << std::setw(20) << heterogeneousMechanism_.I_C2H2();
					fOutput << std::setprecision(9) << std::setw(20) << heterogeneousMechanism_.I_C6H6();
				}

				// Heterogeneous reaction rates: contributions to the bulk density
				{
					fOutput << std::setprecision(9) << std::setw(20) << delta_rhobulk_(point);	// [kg/m3]
					for (int j = 0; j < heterogeneousMechanism_.r().size(); j++)
					{
						fOutput << std::setprecision(9) << std::setw(20) << delta_rhobulk_due_to_single_reaction_[j](point);				// [kg/m3]
						fOutput << std::setprecision(9) << std::setw(20) << delta_rhobulk_due_to_single_reaction_over_rhobulk_[j](point);	// [-]
					}
				}

				// Finalize
				fOutput << std::endl;
			}

		fOutput.close();
	}

	void Reactor2D::PrintTecplotDetailedKinetics(const double t, const std::string name_file)
	{
		OpenSMOKE::OpenSMOKEVectorDouble yy(nc_);
		OpenSMOKE::OpenSMOKEVectorDouble xx(nc_);

		std::ofstream fOutput(name_file.c_str(), std::ios::out);

		// Tecplot file
		{
			fOutput << "Title=Solution" << std::endl;
			fOutput << "Variables=";
			fOutput << "\"x[mm]\"" << ", ";
			fOutput << "\"y[mm]\"" << ", ";
			fOutput << "\"time[h]\"" << ", ";
			fOutput << "\"T[K]\"" << ", ";
			fOutput << "\"P[Pa]\"" << ", ";
			fOutput << "\"eps[-]\"" << ", ";
			fOutput << "\"rhoBulk[kg/m3]\"" << ", ";
			fOutput << "\"Sv[1/m]\"" << ", ";
			fOutput << "\"rpore[micron]\"" << ", ";
			fOutput << "\"K[m2]\"" << ", ";
			fOutput << "\"TauBulk[-]\"" << ", ";
			fOutput << "\"TauKnudsen[-]\"" << ", ";
			fOutput << "\"TauViscous[-]\"" << ", ";
			fOutput << "\"rhoGas[kg/m3]\"" << ", ";
			fOutput << "\"mwGas[kg/kmol]\"" << ", ";

			// Species (mass fractions)
			for (unsigned int j = 0; j < nc_; j++)
				fOutput << "\"" << thermodynamicsMap_.NamesOfSpecies()[j] << "\", ";

			// Surface densities
			for (unsigned int j = 0; j < surf_np_; j++)
				fOutput << "\"" << "Gamma[kmol/m2]" << j << "\", ";

			// Surface species fractions
			for (unsigned int j = 0; j < surf_nc_; j++)
				fOutput << "\"" << thermodynamicsSurfaceMap_.NamesOfSpecies()[nc_ + j] << "\", ";

			// Heterogeneous reaction rates
			{
				fOutput << "\"rDep[kg/m2/s]\"" << ", ";
				fOutput << "\"rDep[kg/m3/s]\"" << ", ";
				fOutput << "\"rDep[m/s]\"" << ", ";
			}

			// Heterogeneous formation rates [kmol/m3/s] 
			{
				fOutput << "\"RhetCH4[kml/m3/s]\"" << ", ";
				fOutput << "\"RhetC2H4[kml/m3/s]\"" << ", ";
				fOutput << "\"RhetC2H2[kml/m3/s]\"" << ", ";
				fOutput << "\"RhetC6H6[kml/m3/s]\"" << ", ";
				fOutput << "\"RhetH2[kml/m3/s]\"" << ", ";
			}

			// Finalize
			fOutput << std::endl;

			// Titles
			fOutput << "Zone I=" << nx_ << ", J=" << ny_ << ", F=POINT" << std::endl;
		}

		for (unsigned int k = 0; k < ny_; k++)
			for (unsigned int i = 0; i < nx_; i++)
			{
				const int point = k*nx_ + i;

				// Calculate the reaction and formation rates of heterogeneous reactions
				{
					porousMedium_.SetPorosity(epsilon_(point));
					porousMedium_.SetTemperature(T_(point));
					porousMedium_.SetPressure(P_(point));

					// Homogeneous contributions
					{
						thermodynamicsMap_.SetPressure(P_(point));
						thermodynamicsMap_.SetTemperature(T_(point));
						aux_Y.CopyFrom(Y_[point].data());
						thermodynamicsMap_.MoleFractions_From_MassFractions(aux_X.GetHandle(), mw_(point), aux_Y.GetHandle());
						aux_X.CopyTo(X_[point].data());
						const double cTot = P_(point) / PhysicalConstants::R_J_kmol / T_(point); // [kmol/m3]
						Product(cTot, aux_X, &aux_C);
					}

					// Heterogeneous mechanism
					{
						heterogeneousDetailedMechanism_.SetTemperature(T_(point));
						heterogeneousDetailedMechanism_.SetPressure(P_(point));

						for (unsigned int j = 0; j < nc_; j++)
							eigen_C_(j) = aux_C[j + 1];

						for (unsigned int j = 0; j < surf_nc_; j++)
							eigen_Z_(j) = Z_[point](j);

						for (unsigned int j = 0; j < bulk_nc_; j++)
							eigen_a_(j) = 1.;

						for (unsigned int j = 0; j < surf_np_; j++)
							eigen_gamma_(j) = Gamma_[point](j);

						heterogeneousDetailedMechanism_.FormationRates(Sv_(point), eigen_C_, eigen_Z_, eigen_a_, eigen_gamma_);
					}
				}

				// Write grid
				fOutput << std::setprecision(9) << std::setw(20) << grid_x_.x()(i);
				fOutput << std::setprecision(9) << std::setw(20) << grid_y_.x()(k);

				// Write rlevant variables
				fOutput << std::setprecision(9) << std::setw(20) << t / 3600.;
				fOutput << std::setprecision(9) << std::setw(20) << T_(point);
				fOutput << std::setprecision(9) << std::setw(20) << P_(point);

				// Write solid phase properties
				fOutput << std::setprecision(9) << std::setw(20) << epsilon_(point);
				fOutput << std::setprecision(9) << std::setw(20) << rho_bulk_(point);
				fOutput << std::setprecision(9) << std::setw(20) << Sv_(point);
				fOutput << std::setprecision(9) << std::setw(20) << rp_(point)*1.e6;
				fOutput << std::setprecision(9) << std::setw(20) << permeability_(point);
				fOutput << std::setprecision(9) << std::setw(20) << eta_bulk_(point);
				fOutput << std::setprecision(9) << std::setw(20) << eta_knudsen_(point);
				fOutput << std::setprecision(9) << std::setw(20) << eta_viscous_(point);

				// Write gas phase properties
				fOutput << std::setprecision(9) << std::setw(20) << rho_gas_(point);
				fOutput << std::setprecision(9) << std::setw(20) << mw_(point);

				// Write species mass fractions
				for (unsigned int j = 0; j < nc_; j++)
					fOutput << std::setprecision(9) << std::setw(20) << Y_[point](j);

				// Write surface species fractions
				for (unsigned int j = 0; j < surf_np_; j++)
					fOutput << std::setprecision(9) << std::setw(20) << Gamma_[point](j);

				// Write surface species fractions
				for (unsigned int j = 0; j < surf_nc_; j++)
					fOutput << std::setprecision(9) << std::setw(20) << Z_[point](j);

				// Heterogeneous reaction rates
				{
					fOutput << std::setprecision(9) << std::setw(20) << heterogeneousDetailedMechanism_.r_deposition_per_unit_area()*heterogeneousDetailedMechanism_.mw_carbon();		// [kg/m2/s]
					fOutput << std::setprecision(9) << std::setw(20) << heterogeneousDetailedMechanism_.r_deposition_per_unit_volume()*heterogeneousDetailedMechanism_.mw_carbon();		// [kg/m3/s]
					fOutput << std::setprecision(9) << std::setw(20) << heterogeneousDetailedMechanism_.r_deposition_per_unit_volume()*heterogeneousDetailedMechanism_.mw_carbon() /
																							(porousMedium_.Sv() / rho_graphite_);														// [m/s]																		// [m/s]
				}

				// Heterogeneous formation rates
				{
					fOutput << std::setprecision(9) << std::setw(20) << -heterogeneousDetailedMechanism_.Rgas()(thermodynamicsSurfaceMap_.IndexOfSpecies("CH4"));
					fOutput << std::setprecision(9) << std::setw(20) << -heterogeneousDetailedMechanism_.Rgas()(thermodynamicsSurfaceMap_.IndexOfSpecies("C2H4"));
					fOutput << std::setprecision(9) << std::setw(20) << -heterogeneousDetailedMechanism_.Rgas()(thermodynamicsSurfaceMap_.IndexOfSpecies("C2H2"));
					fOutput << std::setprecision(9) << std::setw(20) << -heterogeneousDetailedMechanism_.Rgas()(thermodynamicsSurfaceMap_.IndexOfSpecies("C6H6"));
					fOutput << std::setprecision(9) << std::setw(20) << -heterogeneousDetailedMechanism_.Rgas()(thermodynamicsSurfaceMap_.IndexOfSpecies("H2"));
				}

				// Finalize
				fOutput << std::endl;
			}

		fOutput.close();
	}

	void Reactor2D::PrintSolution(const double t, const std::string name_file)
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
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "rhoB[kg/m3]", count);
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
			const unsigned int iy = int(i/nx_);
			const unsigned int ix = i - iy*nx_;

			OpenSMOKE::OpenSMOKEVectorDouble yy(nc_);
			OpenSMOKE::OpenSMOKEVectorDouble xx(nc_);

			for (unsigned int j = 0; j < nc_; j++)
				yy[j + 1] = Y_[i](j);

			thermodynamicsMap_.SetPressure(P_(i));
			thermodynamicsMap_.SetTemperature(T_(i));
			thermodynamicsMap_.MoleFractions_From_MassFractions(xx.GetHandle(), mw_(i), yy.GetHandle());

			fOutput << std::setprecision(9) << std::setw(20) << t;
			fOutput << std::setprecision(9) << std::setw(20) << grid_x_.x()[ix] * 1000.;
			fOutput << std::setprecision(9) << std::setw(20) << grid_y_.x()[iy] * 1000.;
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

	void Reactor2D::PrintDiffusionCoefficients(const double t, const std::string name_file)
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
				thermodynamicsMap_.MoleFractions_From_MassFractions(aux_X.GetHandle(), mw_(i), aux_Y.GetHandle());
				aux_X.CopyTo(X_[i].data());
				for (unsigned int j = 0; j < nc_; j++)
					aux_eigen(j) = aux_X[j + 1];
				porousMedium_.EffectiveDiffusionCoefficients(aux_eigen);
			}

			fOutput << std::setprecision(9) << std::setw(20) << t;
			fOutput << std::setprecision(9) << std::setw(20) << i;// grid_.x()[i] * 1000.;
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

	void Reactor2D::PrintHomogeneousRates(const double t, const std::string name_file)
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
				thermodynamicsMap_.MoleFractions_From_MassFractions(aux_X.GetHandle(), mw_(i), aux_Y.GetHandle());
				aux_X.CopyTo(X_[i].data());
				const double cTot = P_(i) / PhysicalConstants::R_J_kmol / T_(i); // [kmol/m3]
				Product(cTot, aux_X, &aux_C);
				kineticsMap_.SetTemperature(T_(i));
				kineticsMap_.SetPressure(P_(i));
				kineticsMap_.ReactionRates(aux_C.GetHandle());
				kineticsMap_.FormationRates(aux_R.GetHandle());
			}

			fOutput << std::setprecision(9) << std::setw(20) << t;
			fOutput << std::setprecision(9) << std::setw(20) << i;// grid_.x()[i] * 1000.;
			fOutput << std::setprecision(9) << std::setw(20) << T_(i);
			fOutput << std::setprecision(9) << std::setw(20) << P_(i);

			for (unsigned int j = 0; j < nc_; j++)
					fOutput << std::setprecision(9) << std::setw(20) << aux_R[j+1];

			fOutput << std::endl;
		}

		fOutput.close();
	}

	void Reactor2D::PrintIntegralHomogeneousRates(const double t, const std::string name_file)
	{
		Eigen::MatrixXd	omegadot_from_homogeneous_(np_, aux_Y.Size());
		Eigen::MatrixXd	omegadot_from_heterogeneous_(np_, aux_Y.Size());

		for (unsigned int i = 0; i < np_; i++)
		{
			// Concentrations
			thermodynamicsMap_.SetTemperature(T_(i));
			thermodynamicsMap_.SetPressure(P_(i));
			aux_Y.CopyFrom(Y_[i].data());
			thermodynamicsMap_.MoleFractions_From_MassFractions(aux_X.GetHandle(), mw_(i), aux_Y.GetHandle());
			aux_X.CopyTo(X_[i].data());
			const double cTot = P_(i) / PhysicalConstants::R_J_kmol / T_(i); // [kmol/m3]
			Product(cTot, aux_X, &aux_C);

			// Formation rates: homogeneous reactions
			if (heterogeneousDetailedMechanism_.homogeneous_reactions() == true)
			{
				kineticsMap_.SetTemperature(T_(i));
				kineticsMap_.SetPressure(P_(i));
				kineticsMap_.ReactionRates(aux_C.GetHandle());
				kineticsMap_.FormationRates(aux_R.GetHandle());
				OpenSMOKE::ElementByElementProduct(aux_R.Size(), aux_R.GetHandle(), thermodynamicsMap_.MWs().data(), aux_R.GetHandle()); // [kg/m3/s]
				
				for (unsigned int j = 0; j < nc_; j++)
					omegadot_from_homogeneous_(i, j) = epsilon_(i)*aux_R[j+1];
			}

			// Formation rates: heterogenous reactions
			if (heterogeneousDetailedMechanism_.heterogeneous_reactions() == true)
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

				const double omega_deposition_per_unit_volume = heterogeneousMechanism_.r_deposition_per_unit_volume()*heterogeneousMechanism_.mw_carbon();		// [kg/m3/s]
				for (unsigned int j = 0; j < nc_; j++)
					omegadot_from_heterogeneous_(i,j) = heterogeneousDetailedMechanism_.Rgas()(j)*thermodynamicsMap_.MW(j)
														+ Y_[i](j)*omega_deposition_per_unit_volume;															// [kg/m3/s]
			}
		}

		// Print on file
		{
			const double Ri = grid_x_.x()(0);
			const double Re = grid_x_.x()(nx_ - 1);
			const double total_volume = boost::math::constants::pi<double>()*(Re*Re - Ri * Ri)*grid_y_.L();

			std::ofstream fOutputXML(name_file.c_str(), std::ios::out);
			fOutputXML.setf(std::ios::scientific);

			fOutputXML << "<?xml version=\"1.0\" encoding=\"utf-8\"?>" << std::endl;
			fOutputXML << "<opensmoke version=\"0.1a\">" << std::endl;

			fOutputXML << "<number-species>" << std::endl;
			fOutputXML << nc_ << std::endl;
			fOutputXML << "</number-species>" << std::endl;

			fOutputXML << "<total-volume>" << std::endl;
			fOutputXML << total_volume << std::endl;
			fOutputXML << "</total-volume>" << std::endl;

			fOutputXML << "<slice-volume>" << std::endl;
			fOutputXML << total_volume*5./360. << std::endl;
			fOutputXML << "</slice-volume>" << std::endl;

			// New evaluation of source terms
			{
				// Correction coefficient
				const double cc = 1.0;

				// Time interval (s)
				const double delta_time = dae_time_interval_;
				fOutputXML << "<delta_time>" << std::endl;
				fOutputXML << delta_time << std::endl;
				fOutputXML << "</delta_time>" << std::endl;

				// Mass source terms (from time history, per unit of volume)
				std::cout << "Writing source-terms section..." << std::endl;
				{
					fOutputXML << "<source-terms>" << std::endl;
					fOutputXML << "<!--Species Hom.(kg/m3/s) Het.(kg/m3/s) Net(kg/m3/s)-->" << std::endl;
					double sum_homogeneous = 0.;
					double sum_heterogeneous = 0.;
					for (unsigned int j = 0; j < nc_; j++)
					{
						sum_homogeneous += homogeneous_total_mass_source_(j);
						sum_heterogeneous += heterogeneous_total_mass_source_(j);

						fOutputXML << std::left << std::setw(24) << thermodynamicsMap_.NamesOfSpecies()[j];

						fOutputXML << std::right << std::setprecision(9) << std::setw(20) << cc * homogeneous_total_mass_source_(j) / delta_time / total_volume;
						fOutputXML << std::right << std::setprecision(9) << std::setw(20) << cc * heterogeneous_total_mass_source_(j) / delta_time / total_volume;
						fOutputXML << std::right << std::setprecision(9) << std::setw(20) << cc * (homogeneous_total_mass_source_(j) + heterogeneous_total_mass_source_(j)) / delta_time / total_volume;

						fOutputXML << std::endl;
					}
					fOutputXML << "</source-terms>" << std::endl;

					fOutputXML << "<total-mass-source-terms>" << std::endl;
					fOutputXML << "<!--Species Hom.(kg/s) Het.(kg/s) Net(kg/s)-->" << std::endl;
					fOutputXML << std::right << std::setprecision(9) << std::setw(20) << cc * sum_homogeneous / delta_time / total_volume;
					fOutputXML << std::right << std::setprecision(9) << std::setw(20) << cc * sum_heterogeneous / delta_time / total_volume;
					fOutputXML << std::right << std::setprecision(9) << std::setw(20) << cc * (sum_homogeneous + sum_heterogeneous) / delta_time / total_volume;
					fOutputXML << std::endl;
					fOutputXML << "</total-mass-source-terms>" << std::endl;
				}

				// Mass exchanged terms (from time history, per unit of volume)
				std::cout << "Writing exchanged-terms section..." << std::endl;
				{
					fOutputXML << "<exchanged-terms>" << std::endl;
					fOutputXML << "<!--Species dummy dummy Net(kg/m3/s)-->" << std::endl;
					double sum_homogeneous = 0.;
					for (unsigned int j = 0; j < nc_; j++)
					{
						sum_homogeneous += total_mass_exchanged_(j);

						fOutputXML << std::left << std::setw(24) << thermodynamicsMap_.NamesOfSpecies()[j];

						fOutputXML << std::right << std::setprecision(9) << std::setw(20) << 0.;
						fOutputXML << std::right << std::setprecision(9) << std::setw(20) << 0.;
						fOutputXML << std::right << std::setprecision(9) << std::setw(20) << total_mass_exchanged_(j) / delta_time / total_volume;

						fOutputXML << std::endl;
					}
					fOutputXML << "</exchanged-terms>" << std::endl;
					
					fOutputXML << "<total-mass-exchanged-terms>" << std::endl;
					fOutputXML << "<!--Species dummy dummy Net(kg/s)-->" << std::endl;
					fOutputXML << std::right << std::setprecision(9) << std::setw(20) << 0.;
					fOutputXML << std::right << std::setprecision(9) << std::setw(20) << 0.;
					fOutputXML << std::right << std::setprecision(9) << std::setw(20) << sum_homogeneous / delta_time / total_volume;
					fOutputXML << std::endl;
					fOutputXML << "</total-mass-exchanged-terms>" << std::endl;
				}

				// Source terms (heterogeneous contribution only, instantaneous and from time history)
				std::cout << "Writing het-source-terms section..." << std::endl;
				{
					fOutputXML << "<het-source-terms>" << std::endl;
					fOutputXML << "<!--Species Het.Inst.(kg/m3/s) Het.Int.(kg/s) Het.Int.(kg/m3/s)-->" << std::endl;

					double sum_heterogeneous_integral = 0.;
					double sum_heterogeneous_instantaneous = 0.;
					for (unsigned int j = 0; j < nc_; j++)
					{
						sum_heterogeneous_integral += heterogeneous_total_mass_source_(j);

						const double heterogeneous_instantaneous = VolumeAveraged(omegadot_from_heterogeneous_.col(j));
						sum_heterogeneous_instantaneous += heterogeneous_instantaneous;

						fOutputXML << std::left << std::setw(24) << thermodynamicsMap_.NamesOfSpecies()[j];

						fOutputXML << std::right << std::setprecision(9) << std::setw(20) << cc * heterogeneous_instantaneous;
						fOutputXML << std::right << std::setprecision(9) << std::setw(20) << cc * heterogeneous_total_mass_source_(j) / delta_time;
						fOutputXML << std::right << std::setprecision(9) << std::setw(20) << cc * heterogeneous_total_mass_source_(j) / delta_time / total_volume;

						fOutputXML << std::endl;
					}
					fOutputXML << "</het-source-terms>" << std::endl;

					fOutputXML << "<total-het-source-terms>" << std::endl;
					fOutputXML << "<!--Species Het.Inst.(kg/m3/s) Het.Int.(kg/s) Het.Int.(kg/m3/s)-->" << std::endl;
					fOutputXML << std::right << std::setprecision(9) << std::setw(20) << cc * sum_heterogeneous_instantaneous;
					fOutputXML << std::right << std::setprecision(9) << std::setw(20) << cc * sum_heterogeneous_integral / delta_time;
					fOutputXML << std::right << std::setprecision(9) << std::setw(20) << cc * sum_heterogeneous_integral / delta_time / total_volume;
					fOutputXML << std::endl;
					fOutputXML << "</total-het-source-terms>" << std::endl;
				}

				// Source terms (homogeneous contribution only, instantaneous and from time history)
				std::cout << "Writing hom-source-terms section..." << std::endl;
				{
					fOutputXML << "<hom-source-terms>" << std::endl;
					fOutputXML << "<!--Species Hom.Inst.(kg/m3/s) Hom.Int.(kg/s) Hom.Int.(kg/m3/s)-->" << std::endl;

					double sum_homogeneous_integral = 0.;
					double sum_homogeneous_instantaneous = 0.;
					for (unsigned int j = 0; j < nc_; j++)
					{
						sum_homogeneous_integral += homogeneous_total_mass_source_(j);

						const double homogeneous_instantaneous = VolumeAveraged(omegadot_from_homogeneous_.col(j));
						sum_homogeneous_instantaneous += homogeneous_instantaneous;

						fOutputXML << std::left << std::setw(24) << thermodynamicsMap_.NamesOfSpecies()[j];

						fOutputXML << std::right << std::setprecision(9) << std::setw(20) << cc * homogeneous_instantaneous;
						fOutputXML << std::right << std::setprecision(9) << std::setw(20) << cc * homogeneous_total_mass_source_(j) / delta_time;
						fOutputXML << std::right << std::setprecision(9) << std::setw(20) << cc * homogeneous_total_mass_source_(j) / delta_time / total_volume;

						fOutputXML << std::endl;
					}
					fOutputXML << "</hom-source-terms>" << std::endl;

					fOutputXML << "<total-hom-source-terms>" << std::endl;
					fOutputXML << "<!--Species Hom.Inst.(kg/m3/s) Hom.Int.(kg/s) Hom.Int.(kg/m3/s)-->" << std::endl;
					fOutputXML << std::right << std::setprecision(9) << std::setw(20) << cc * sum_homogeneous_instantaneous;
					fOutputXML << std::right << std::setprecision(9) << std::setw(20) << cc * sum_homogeneous_integral / delta_time;
					fOutputXML << std::right << std::setprecision(9) << std::setw(20) << cc * sum_homogeneous_integral / delta_time / total_volume;
					fOutputXML << std::endl;
					fOutputXML << "</total-hom-source-terms>" << std::endl;
				}
				
				// Mass source terms (integral)
				std::cout << "Writing mass-source-terms section..." << std::endl;
				{
					fOutputXML << "<mass-source-terms>" << std::endl;
					fOutputXML << "<!--Species Hom.(kg/s) Het.(kg/s) Net(kg/s)-->" << std::endl;
					double sum_homogeneous = 0.;
					double sum_heterogeneous = 0.;
					for (unsigned int j = 0; j < nc_; j++)
					{
						sum_homogeneous += homogeneous_total_mass_source_(j);
						sum_heterogeneous += heterogeneous_total_mass_source_(j);

						fOutputXML << std::left << std::setw(24) << thermodynamicsMap_.NamesOfSpecies()[j];

						fOutputXML << std::right << std::setprecision(9) << std::setw(20) << cc*homogeneous_total_mass_source_(j) / delta_time;
						fOutputXML << std::right << std::setprecision(9) << std::setw(20) << cc*heterogeneous_total_mass_source_(j) / delta_time;
						fOutputXML << std::right << std::setprecision(9) << std::setw(20) << cc*(homogeneous_total_mass_source_(j) + heterogeneous_total_mass_source_(j)) / delta_time;

						fOutputXML << std::endl;
					}
					fOutputXML << "</mass-source-terms>" << std::endl;

					fOutputXML << "<total-mass-source-terms>" << std::endl;
					fOutputXML << "<!--Species Het.(kg/s) Het.(kg/s) Net(kg/s)-->" << std::endl;
					fOutputXML << std::right << std::setprecision(9) << std::setw(20) << cc*sum_homogeneous / delta_time;
					fOutputXML << std::right << std::setprecision(9) << std::setw(20) << cc*sum_heterogeneous / delta_time;
					fOutputXML << std::right << std::setprecision(9) << std::setw(20) << cc*(sum_homogeneous + sum_heterogeneous) / delta_time;
					fOutputXML << std::endl;
					fOutputXML << "</total-mass-source-terms>" << std::endl;
				}

				// Old evaluation of source terms
				std::cout << "Writing source-terms-old section..." << std::endl;
				{
					fOutputXML << "<source-terms-old>" << std::endl;
					fOutputXML << "<!--Species Hom.(kg/m3/s) Het.(kg/m3/s) Net(kg/m3/s)-->" << std::endl;
					double sum_homogeneous = 0.;
					double sum_heterogeneous = 0.;
					for (unsigned int j = 0; j < nc_; j++)
					{
						const double homogeneous = VolumeAveraged(omegadot_from_homogeneous_.col(j));
						const double heterogeneous = VolumeAveraged(omegadot_from_heterogeneous_.col(j));
						sum_homogeneous += homogeneous;
						sum_heterogeneous += heterogeneous;
						fOutputXML << std::left << std::setw(24) << thermodynamicsMap_.NamesOfSpecies()[j];
						fOutputXML << std::right << std::setprecision(9) << std::setw(20) << homogeneous;
						fOutputXML << std::right << std::setprecision(9) << std::setw(20) << heterogeneous;
						fOutputXML << std::right << std::setprecision(9) << std::setw(20) << homogeneous + heterogeneous;
						fOutputXML << std::endl;
					}
					fOutputXML << "</source-terms-old>" << std::endl;

					fOutputXML << "<total-source-terms-old>" << std::endl;
					fOutputXML << "<!--Species Hom.(kg/m3/s) Het.(kg/m3/s) Net(kg/m3/s)-->" << std::endl;
					fOutputXML << std::right << std::setprecision(9) << std::setw(20) << sum_homogeneous;
					fOutputXML << std::right << std::setprecision(9) << std::setw(20) << sum_heterogeneous;
					fOutputXML << std::right << std::setprecision(9) << std::setw(20) << sum_homogeneous + sum_heterogeneous;
					fOutputXML << std::endl;
					fOutputXML << "</total-source-terms-old>" << std::endl;
				}
			}

			fOutputXML << "</opensmoke>" << std::endl;
			fOutputXML.close();
		}
	}

	void Reactor2D::PrintHeterogeneousRates(const double t, const std::string name_file)
	{
		if (detailed_heterogeneous_kinetics_ == true)
			PrintDetailedHeterogeneousRates(t, name_file);
		else
			PrintGlobalHeterogeneousRates(t, name_file);
	}

	void Reactor2D::PrintGlobalHeterogeneousRates(const double t, const std::string name_file)
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
				thermodynamicsMap_.MoleFractions_From_MassFractions(aux_X.GetHandle(), mw_(i), aux_Y.GetHandle());
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
			fOutput << std::setprecision(9) << std::setw(width) << i;// grid_.x()[i] * 1000.;
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

	void Reactor2D::PrintDetailedHeterogeneousRates(const double t, const std::string name_file)
	{
		// TODO
	}

	void Reactor2D::PrintROPA(const double t, const std::string name_file)
	{
		// Open ROPA file
		std::ofstream fROPA(name_file.c_str(), std::ios::out);

		// Write head in file
		ropa_->WriteHead(fROPA, "Reactor2D");

		// ROPA at mouth section
		{
			OpenSMOKE::OpenSMOKEVectorDouble x(nc_);
			OpenSMOKE::OpenSMOKEVectorDouble omega(nc_);
			OpenSMOKE::OpenSMOKEVectorDouble c(nc_);
			OpenSMOKE::OpenSMOKEVectorDouble z(surf_nc_);
			OpenSMOKE::OpenSMOKEVectorDouble a(bulk_nc_);
			OpenSMOKE::OpenSMOKEVectorDouble gamma(surf_np_);

			std::vector<double> list_points(5);
			list_points[0] = np_/2-1;
			list_points[1] = 0;
			list_points[2] = grid_x_.Np()-1;
			list_points[3] = (np_ - grid_x_.Np()) - 1;
			list_points[4] = np_-1;

			std::vector<std::string> list_points_names(5);
			list_points_names[0] = "Center";
			list_points_names[1] = "SE Corner";
			list_points_names[2] = "SW Corner";
			list_points_names[3] = "NE Corner";
			list_points_names[4] = "NW Corner";

			for (unsigned int i = 0; i<list_points.size(); i++)
			{
				const unsigned int point = static_cast<unsigned int>(list_points[i]);

				fROPA << std::endl;
				fROPA << "-------------------------------------------------------------------" << std::endl;
				fROPA << " * ROPA at: " << list_points_names[i] << std::endl;
				fROPA << "-------------------------------------------------------------------" << std::endl;
				fROPA << std::endl;

				// Molar fractions
				double mw;
				for (unsigned int j = 0; j < nc_; j++)
					omega[j + 1] = Y_[point](j);
				thermodynamicsMap_.MoleFractions_From_MassFractions(x.GetHandle(), mw, omega.GetHandle());

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

	void Reactor2D::PrintLabelMonitoringFile()
	{
		unsigned int count = 1;
		const unsigned int width = 20;
		const unsigned int width_increased = 26;

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

		// Heterogeneous reaction rates: contributions to the bulk density
		OpenSMOKE::PrintTagOnASCIILabel(width_increased, fMonitoring_, "dRhoB[kg/m3]", count);
		for (int j = 0; j < heterogeneousMechanism_.r().size(); j++)
			OpenSMOKE::PrintTagOnASCIILabel(width_increased, fMonitoring_, "dRhoB_" + heterogeneousMechanism_.tags()[j] + "[kg/m3]", count);

		fMonitoring_ << std::endl;
	}

	void Reactor2D::Print(const double t, const double* y)
	{
		if (equations_set_ == EQUATIONS_SET_COMPLETE)
		{
			if (count_dae_video_%n_steps_video_ == 1)
			{
				if (count_dae_video_ % (n_steps_video_ * 50) == 1)
				{
					std::cout << std::endl;
					std::cout << std::left << std::setw(14) << "Time[s]";		// [s]
					std::cout << std::left << std::setw(14) << "Time[h]";		// [h]
					std::cout << std::left << std::setw(14) << "Rho(M)[kg/m3]";	// [kg/m3]
					std::cout << std::left << std::setw(14) << "Rho(Std)[kg/m3]";	// [kg/m3]
					std::cout << std::left << std::setw(14) << "Rdep[mu/h]";	// [micron/h]
					std::cout << std::left << std::setw(16) << "Rdep[kg/m3/h]";	// [kg/m3/h]
					std::cout << std::left << std::setw(16) << "epsMin[-]";		// [-]
					std::cout << std::left << std::setw(16) << "epsMax[-]";		// [-]

					if (detailed_heterogeneous_kinetics_ == true)
					{
						for (unsigned int i = 0; i < surf_np_; i++)
							std::cout << std::left << std::setw(16) << "Gamma[kmol/m2]";
						std::cout << std::left << std::setw(14) << "Error";
					}

					std::cout << std::endl;
				}


				const double rho_bulk_mean = VolumeAveraged(rho_bulk_);
				const double rho_bulk_std = VolumeStandardDeviation(rho_bulk_mean, rho_bulk_);
				const double r_deposition_per_unit_area_mean = VolumeAveraged(omega_deposition_per_unit_area_);		// [kg/m2/s]
				const double r_deposition_per_unit_volume_mean = VolumeAveraged(omega_deposition_per_unit_volume_);	// [kg/m3/s]

				std::cout << std::left << std::setw(14) << std::scientific << std::setprecision(6) << t;		// [s]
				std::cout << std::left << std::setw(14) << std::scientific << std::setprecision(6) << t / 3600.;
				std::cout << std::left << std::setw(14) << std::fixed << std::setprecision(6) << rho_bulk_mean;
				std::cout << std::left << std::setw(14) << std::fixed << std::setprecision(6) << rho_bulk_std;
				std::cout << std::left << std::setw(14) << std::fixed << std::setprecision(6) << r_deposition_per_unit_area_mean / rho_graphite_*1e6*3600.;	// [micron/h]
				std::cout << std::left << std::setw(16) << std::fixed << std::setprecision(6) << r_deposition_per_unit_volume_mean *3600.;			// [kg/m3/h]
				std::cout << std::left << std::setw(14) << std::fixed << std::setprecision(6) << epsilon_.minCoeff();
				std::cout << std::left << std::setw(16) << std::fixed << std::setprecision(6) << epsilon_.maxCoeff();

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
				const int width_increased = 26;

				const double T_mean = VolumeAveraged(T_);
				const double T_std = VolumeStandardDeviation(T_mean, T_);

				const double P_mean = VolumeAveraged(P_);
				const double P_std = VolumeStandardDeviation(P_mean, P_);

				const double epsilon_mean = VolumeAveraged(epsilon_);
				const double epsilon_std = VolumeStandardDeviation(epsilon_mean, epsilon_);

				const double rho_bulk_mean = VolumeAveraged(rho_bulk_);
				const double rho_bulk_std = VolumeStandardDeviation(rho_bulk_mean, rho_bulk_);

				const double Sv_mean = VolumeAveraged(Sv_);
				const double Sv_std = VolumeStandardDeviation(Sv_mean, Sv_);

				const double rp_mean = VolumeAveraged(rp_);
				const double rp_std = VolumeStandardDeviation(rp_mean, rp_);

				const double permeability_mean = VolumeAveraged(permeability_);
				const double permeability_std = VolumeStandardDeviation(permeability_mean, permeability_);

				const double eta_bulk_mean = VolumeAveraged(eta_bulk_);
				const double eta_bulk_std = VolumeStandardDeviation(eta_bulk_mean, eta_bulk_);

				const double eta_knudsen_mean = VolumeAveraged(eta_knudsen_);
				const double eta_knudsen_std = VolumeStandardDeviation(eta_knudsen_mean, eta_knudsen_);

				const double eta_viscous_mean = VolumeAveraged(eta_viscous_);
				const double eta_viscous_std = VolumeStandardDeviation(eta_viscous_mean, eta_viscous_);

				const double r_deposition_per_unit_area_mean = VolumeAveraged(omega_deposition_per_unit_area_);
				const double r_deposition_per_unit_area_std = VolumeStandardDeviation(r_deposition_per_unit_area_mean, omega_deposition_per_unit_area_);

				const double r_deposition_per_unit_volume_mean = VolumeAveraged(omega_deposition_per_unit_volume_);
				const double r_deposition_per_unit_volume_std = VolumeStandardDeviation(r_deposition_per_unit_volume_mean, omega_deposition_per_unit_volume_);

				const double delta_rhobulk_mean = VolumeAveraged(delta_rhobulk_);
				Eigen::VectorXd delta_rhobulk_due_to_single_reaction_mean(heterogeneousMechanism_.r().size());
				for (int j = 0; j < heterogeneousMechanism_.r().size(); j++)
					delta_rhobulk_due_to_single_reaction_mean(j) = VolumeAveraged(delta_rhobulk_due_to_single_reaction_[j]);

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

				// Heterogeneous reaction rates: contributions to the bulk density
				fMonitoring_ << std::left << std::setw(width_increased) << std::fixed << std::setprecision(4) << delta_rhobulk_mean;	// [kg/m3]
				for (int j = 0; j < heterogeneousMechanism_.r().size(); j++)
					fMonitoring_ << std::left << std::setw(width_increased) << std::fixed << std::setprecision(4) << delta_rhobulk_due_to_single_reaction_mean(j);	// [kg/m3]

				fMonitoring_ << std::endl;

				count_file_ = 0;
			}

			// Post-processing
			{
				const double delta_time = t - t_old_;
				const double coefficient = heterogeneousMechanism_.mw_carbon()*delta_time;

				for (unsigned int i = 0; i < np_; i++)
				{
					// Porous medium
					porousMedium_.SetPorosity(epsilon_(i));
					porousMedium_.SetTemperature(T_(i));
					porousMedium_.SetPressure(P_(i));

					// Thermodynamics
					thermodynamicsMap_.SetPressure(P_(i));
					thermodynamicsMap_.SetTemperature(T_(i));

					// Concentrations
					aux_Y.CopyFrom(Y_[i].data());
					thermodynamicsMap_.MoleFractions_From_MassFractions(aux_X.GetHandle(), mw_(i), aux_Y.GetHandle());
					aux_X.CopyTo(X_[i].data());
					const double cTot = P_(i) / PhysicalConstants::R_J_kmol / T_(i); // [kmol/m3]
					Product(cTot, aux_X, &aux_C);

					// Heterogeneous reactions
					heterogeneousMechanism_.SetTemperature(T_(i));
					heterogeneousMechanism_.SetPressure(P_(i));
					for (unsigned int j = 0; j < nc_; j++)
						aux_eigen(j) = aux_C[j + 1];
					heterogeneousMechanism_.FormationRates(porousMedium_.Sv(), aux_eigen);

					// Single contributions
					delta_rhobulk_(i) = 0.;
					for (int k = 0; k < heterogeneousMechanism_.r().size(); k++)
					{
						delta_rhobulk_due_to_single_reaction_[k](i) += heterogeneousMechanism_.r_deposition_per_unit_volume_per_single_reaction()(k)*coefficient;
						delta_rhobulk_(i) += delta_rhobulk_due_to_single_reaction_[k](i);
					}

					// Single contributions (normalized)
					for (int k = 0; k < heterogeneousMechanism_.r().size(); k++)
						delta_rhobulk_due_to_single_reaction_over_rhobulk_[k](i) = delta_rhobulk_due_to_single_reaction_[k](i) / (delta_rhobulk_(i) + 1.e-12);
				}
			}

			if (t >= (count_tecplot_ + 1)*tecplot_time_interval_)
			{
				count_tecplot_++;
				std::stringstream current_index; current_index << count_tecplot_;

				std::string tecplot_file = "Solution.tec." + current_index.str();
				PrintTecplot(t, (output_tecplot_folder_ / tecplot_file).string().c_str());

				if (gaseous_phase_ == GASEOUS_PHASE_FROM_PLUG_FLOW)
				{
					std::string plug_flow_file = "PlugFlow.out." + current_index.str();
					plugFlowReactor_.Print(t, (output_plug_flow_folder_ / plug_flow_file).string().c_str());
				}
			}

			if (gaseous_phase_ == GASEOUS_PHASE_FROM_PLUG_FLOW)
			{
				if (plugFlowReactor_.coupling() == true && count_update_plug_flow_ == n_steps_update_plug_flow_)
				{
					const double inlet_T = plugFlowReactor_.inlet_temperature();
					const double inlet_P = plugFlowReactor_.inlet_pressure();
					const Eigen::VectorXd inlet_omega = plugFlowReactor_.inlet_mass_fractions();

					// Plug flow ractor simulation
					std::vector<Eigen::VectorXd> Y_gas_side;
					{
						// Set initial conditions
						plugFlowReactor_.SetInitialConditions(inlet_T, inlet_P, inlet_omega);
						plugFlowReactor_.SetVerboseOutput(false);

						// Set external profile
						{
							Eigen::VectorXd csi_external(ny_ + 2);
							Eigen::MatrixXd omega_external(ny_ + 2, thermodynamicsMap_.NumberOfSpecies());

							// Set csi external
							csi_external(0) = 0.;
							csi_external(ny_ + 1) = 1000.;
							for (unsigned int i = 0; i < ny_; i++)
								csi_external(i + 1) = plugFlowReactor_.inert_length() + grid_y_.x()[i];

							// Internal points
							for (unsigned int i = 0; i < ny_; i++)
							{
								const int point = list_points_east_(i);
								for (unsigned int j = 0; j < nc_; j++)
									omega_external(i + 1, j) = Y_[point](j);
							}

							// First point
							for (unsigned int j = 0; j < nc_; j++)
								omega_external(0, j) = omega_external(1, j);

							// Last point
							for (unsigned int j = 0; j < nc_; j++)
								omega_external(ny_ + 1, j) = omega_external(ny_, j);

							// Set external profile
							plugFlowReactor_.SetExternalMassFractionsProfile(csi_external, omega_external);
						}

						// Solve the plug flow reactor
						plugFlowReactor_.Solve(plugFlowReactor_.last_residence_time());

						// Extract the plug flow history
						Eigen::VectorXd csi(plugFlowReactor_.history_csi().size());
						for (unsigned int i = 0; i < plugFlowReactor_.history_csi().size(); i++)
							csi(i) = plugFlowReactor_.history_csi()[i];

						// Assign the boundary conditions
						PlugFlowReactorCoupledProfiles* profiles = new PlugFlowReactorCoupledProfiles(csi);
						Y_gas_side.resize(2 * grid_x_.Np() + grid_y_.Np());
						for (unsigned int i = 0; i < Y_gas_side.size(); i++)
						{
							Y_gas_side[i].resize(thermodynamicsMap_.NumberOfSpecies());
							Y_gas_side[i].setZero();
						}

						// Case: only east side
						if (plugFlowReactor_.geometric_pattern() == CVI::PlugFlowReactorCoupled::ONE_SIDE)
						{
							for (int i = 0; i < grid_y_.Np(); i++)
							{
								const int point = i + grid_x_.Np();
								const double coordinate = grid_y_.x()(i) + plugFlowReactor_.inert_length();
								profiles->Interpolate(coordinate, plugFlowReactor_.history_Y(), Y_gas_side[point]);
							}
						}

						SetGasSide(inlet_T, inlet_P, Y_gas_side);
					}

					count_update_plug_flow_ = 0;
				}
			}
		}
		else if (equations_set_ == EQUATIONS_SET_ONLYTEMPERATURE)
		{
			std::cout << t << " " << T_(np_ / 2) << std::endl;
		}

		// Update the total amount of consumed/produced species
		{
			Eigen::MatrixXd	omegadot_from_homogeneous_(np_, aux_Y.Size());
			Eigen::MatrixXd	omegadot_from_heterogeneous_(np_, aux_Y.Size());

			for (unsigned int i = 0; i < np_; i++)
			{
				// Concentrations
				thermodynamicsMap_.SetTemperature(T_(i));
				thermodynamicsMap_.SetPressure(P_(i));
				aux_Y.CopyFrom(Y_[i].data());
				thermodynamicsMap_.MoleFractions_From_MassFractions(aux_X.GetHandle(), mw_(i), aux_Y.GetHandle());
				aux_X.CopyTo(X_[i].data());
				const double cTot = P_(i) / PhysicalConstants::R_J_kmol / T_(i); // [kmol/m3]
				Product(cTot, aux_X, &aux_C);

				// Formation rates: homogeneous reactions
				if (heterogeneousDetailedMechanism_.homogeneous_reactions() == true)
				{
					kineticsMap_.SetTemperature(T_(i));
					kineticsMap_.SetPressure(P_(i));
					kineticsMap_.ReactionRates(aux_C.GetHandle());
					kineticsMap_.FormationRates(aux_R.GetHandle());
					OpenSMOKE::ElementByElementProduct(aux_R.Size(), aux_R.GetHandle(), thermodynamicsMap_.MWs().data(), aux_R.GetHandle()); // [kg/m3/s]

					for (unsigned int j = 0; j < nc_; j++)
						omegadot_from_homogeneous_(i, j) = epsilon_(i) * aux_R[j + 1];
				}

				// Formation rates: heterogeneous reactions
				if (heterogeneousDetailedMechanism_.heterogeneous_reactions() == true)
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

					// Option 1
					//const double omega_deposition_per_unit_volume = heterogeneousMechanism_.r_deposition_per_unit_volume() * heterogeneousMechanism_.mw_carbon();		// [kg/m3/s]
					//for (unsigned int j = 0; j < nc_; j++)
					//	omegadot_from_heterogeneous_(i,j) = heterogeneousDetailedMechanism_.Rgas()(j) * thermodynamicsMap_.MW(j)
					//										+ Y_[i](j) * omega_deposition_per_unit_volume;						// [kg/m3/s]

					// Option 2
					const double omega_deposition_per_unit_volume = heterogeneousMechanism_.r_deposition_per_unit_volume() * heterogeneousMechanism_.mw_carbon();		// [kg/m3/s]
					for (unsigned int j = 0; j < nc_; j++)
						omegadot_from_heterogeneous_(i,j) = epsilon_(i) * heterogeneousDetailedMechanism_.Rgas()(j) * thermodynamicsMap_.MW(j);	
				}
			}

			//const double cc = std::tanh(0.05 * t / 3600.);
			const double cc = 1.0;
			for (unsigned int j = 0; j < nc_; j++)
			{
				homogeneous_total_mass_source_(j) += cc*VolumeIntegral(omegadot_from_homogeneous_.col(j)) * (t - t_old_);
				heterogeneous_total_mass_source_(j) += cc*VolumeIntegral(omegadot_from_heterogeneous_.col(j)) * (t - t_old_);
			}

			// Updating fluxes across the faces (XXX)

			{
				// Entering mass flow rates
				std::vector<double> mfr_west(nc_);
				std::vector<double> mfr_east(nc_);
				std::vector<double> mfr_north(nc_);
				std::vector<double> mfr_south(nc_);

				std::fill(mfr_west.begin(), mfr_west.end(), 0);
				std::fill(mfr_east.begin(), mfr_east.end(), 0);
				std::fill(mfr_north.begin(), mfr_north.end(), 0);
				std::fill(mfr_south.begin(), mfr_south.end(), 0);

				// Loop over all the species
				for (unsigned int j = 0; j < nc_; j++)
				{
					// West face
					for (unsigned int i = 0; i < ny_ - 1; i++)
					{
						const int point1 = list_points_west_(i);
						const int point2 = list_points_west_(i + 1);

						const double gradY1 = Y_gas_west_side_[i](j) - Y_[point1](j);
						const double gradY2 = Y_gas_west_side_[i + 1](j) - Y_[point2](j);
						const double gradY = 0.50 * (gradY1 + gradY2);
						const double rho = 0.50 * (rho_gas_(point1) + rho_gas_(point2));
						const double Gamma = 0.50 * (gamma_star_[point1](j) + gamma_star_[point2](j));

						const double r = grid_x_.x()(0);
						const double height = grid_y_.x()(i+1)- grid_y_.x()(i);
						const double area = 2 * boost::math::constants::pi<double>() * r * height;

						mfr_west[j] += rho * Gamma * gradY * area;	// [kg/s]
					}

					// East face
					for (unsigned int i = 0; i < ny_ - 1; i++)
					{
						const int point1 = list_points_east_(i);
						const int point2 = list_points_east_(i + 1);

						const double gradY1 = Y_gas_east_side_[i](j) - Y_[point1](j);
						const double gradY2 = Y_gas_east_side_[i + 1](j) - Y_[point2](j);
						const double gradY = 0.50 * (gradY1 + gradY2);
						const double rho = 0.50 * (rho_gas_(point1) + rho_gas_(point2));
						const double Gamma = 0.50 * (gamma_star_[point1](j) + gamma_star_[point2](j));

						const double r = grid_x_.x()(nx_-1);
						const double height = grid_y_.x()(i + 1) - grid_y_.x()(i);
						const double area = 2 * boost::math::constants::pi<double>() * r * height;

						mfr_east[j] += rho * Gamma * gradY * area;	// [kg/s]
					}

					// North face
					for (unsigned int i = 0; i < nx_ - 1; i++)
					{
						const int point1 = list_points_north_(i);
						const int point2 = list_points_north_(i + 1);

						const double gradY1 = Y_gas_north_side_[i](j) - Y_[point1](j);
						const double gradY2 = Y_gas_north_side_[i + 1](j) - Y_[point2](j);
						const double gradY = 0.50 * (gradY1 + gradY2);
						const double rho = 0.50 * (rho_gas_(point1) + rho_gas_(point2));
						const double Gamma = 0.50 * (gamma_star_[point1](j) + gamma_star_[point2](j));

						const double ri = grid_x_.x()(i);
						const double re = grid_x_.x()(i+1);
						const double area = boost::math::constants::pi<double>() * (re*re-ri*ri);

						mfr_north[j] += rho * Gamma * gradY * area;	// [kg/s]
					}

					// South face
					for (unsigned int i = 0; i < nx_ - 1; i++)
					{
						const int point1 = list_points_south_(i);
						const int point2 = list_points_south_(i + 1);

						const double gradY1 = Y_gas_south_side_[i](j) - Y_[point1](j);
						const double gradY2 = Y_gas_south_side_[i + 1](j) - Y_[point2](j);
						const double gradY = 0.50 * (gradY1 + gradY2);
						const double rho = 0.50 * (rho_gas_(point1) + rho_gas_(point2));
						const double Gamma = 0.50 * (gamma_star_[point1](j) + gamma_star_[point2](j));

						const double ri = grid_x_.x()(i);
						const double re = grid_x_.x()(i + 1);
						const double area = boost::math::constants::pi<double>() * (re * re - ri * ri);

						mfr_south[j] += rho * Gamma * gradY * area;	// [kg/s]
					}
				}

				// Total mass exchanged
				for (unsigned int j = 0; j < nc_; j++)
				{
					const double mfr = mfr_west[j] + mfr_east[j] + mfr_north[j] + mfr_south[j];
					total_mass_exchanged_(j) += mfr * (t - t_old_);
				}
			}
		}

		t_old_ = t;
		count_dae_video_++;
		count_file_++;
		count_update_plug_flow_++;
	}

	void Reactor2D::SparsityPattern(std::vector<unsigned int>& rows, std::vector<unsigned int>& cols)
	{
		std::vector<unsigned int> rows_single;
		std::vector<unsigned int> cols_single;
		OpenSMOKE::SparsityPatternPentadiagonal(np_, nx_, rows_single, cols_single);
		OpenSMOKE::SparsityPatternBlock(np_, block_, rows_single, cols_single, rows, cols);
	}

	void Reactor2D::PrintOnTheScreen(const std::string, const int k, double* v)
	{	
		int count = k;
		for (unsigned int j = 0; j < ny_; j++)
		{
			for (unsigned int i = 0; i < nx_; i++)
			{
				std::cout << v[count] << " ";
				count += block_;
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}

	void Reactor2D::CorrectDifferentialEquations(double* upv, double* resv)
	{
		for (int i = 0; i < differential_equations_.size(); i++)
		{
			const int k = differential_equations_[i];
			resv[k] -= upv[k];
		}
	}

	void Reactor2D::CorrectAlgebraicEquations(double* yp)
	{
		for (int i = 0; i < algebraic_equations_.size(); i++)
		{
			const int k = algebraic_equations_[i];
			yp[k] = 0.;
		}
	}

	void Reactor2D::DiagonalJacobian(const double t, double* y, double* J)
	{
		// Calculated as suggested by Buzzi (private communication)
		const double ZERO_DER = std::sqrt(OPENSMOKE_TINY_FLOAT);
		const double ETA2 = std::sqrt(OpenSMOKE::OPENSMOKE_MACH_EPS_DOUBLE);
		const double TOLR = 100. * OpenSMOKE::OPENSMOKE_MACH_EPS_FLOAT;
		const double TOLA = 1.e-10;

		double* dy_original = new double[NumberOfEquations()];
		double* dy_plus = new double[NumberOfEquations()];

		double* y_plus = new double[NumberOfEquations()];
		for (unsigned int i = 0; i < NumberOfEquations(); i++)
			y_plus[i] = y[i];

		// Fixed values
		Equations(t, y, dy_original);

		// Derivatives with respect to y[kd]
		for (unsigned int kd = 0; kd < NumberOfEquations(); kd++)
		{
			double hf = 1.e0;
			double error_weight = 1. / (TOLA + TOLR*std::fabs(y[kd]));
			double hJ = ETA2 * std::fabs(std::max(y[kd], 1. / error_weight));
			double hJf = hf / error_weight;
			hJ = std::max(hJ, hJf);
			hJ = std::max(hJ, ZERO_DER);

			// This is what is done by Buzzi
			double dy = std::min(hJ, 1.e-3 + 1e-3*std::fabs(y[kd]));
			double udy = 1. / dy;
			y_plus[kd] += dy;
			Equations(t, y_plus, dy_plus);

			//for (int j = 1; j <= y.Size(); j++)
			J[kd] = (dy_plus[kd] - dy_original[kd]) * udy;

			y_plus[kd] = y[kd];
		}

		delete dy_original;
		delete dy_plus;
		delete y_plus;
	}

	void Reactor2D::DiagonalJacobianForIDA(const double alfa, double* J)
	{
		for (int i = 0; i < differential_equations_.size(); i++)
		{
			const int k = differential_equations_[i];
			J[k] -= alfa*1.;
		}
	}

	double Reactor2D::VolumeAveraged(const Eigen::VectorXd& v)
	{
		if (planar_symmetry_ == true)
		{
			const double volume_total = grid_x_.L() * grid_y_.L();
			return VolumeIntegral(v) / volume_total;
		}
		else
		{
			const double Ri = grid_x_.x()(0);
			const double Re = grid_x_.x()(nx_ - 1);
			const double volume_total = boost::math::constants::pi<double>() * (Re * Re - Ri * Ri) * grid_y_.L();
			return VolumeIntegral(v) / volume_total;
		}
	}

	double Reactor2D::VolumeIntegral(const Eigen::VectorXd& v)
	{
		if (planar_symmetry_ == true)
		{
			double sum = 0.;

			// Internal
			for (unsigned int k = 1; k < ny_ - 1; k++)
				for (unsigned int i = 1; i < nx_ - 1; i++)
				{
					const int point = k*nx_ + i;

					const double area = (grid_x_.dxc()(i)*0.50)*(grid_y_.dxc()(k)*0.50);
					sum += v(point)*area;
				}

			// South (zero gradient)
			for (unsigned int i = 1; i < nx_ - 1; i++)
			{
				const int point = list_points_south_(i);
				const double area = (grid_x_.dxc()(i)*0.50)*(grid_y_.dxe()(0)*0.50);
				sum += v(point)*area;
			}

			// North (zero gradient)
			for (unsigned int i = 1; i < nx_ - 1; i++)
			{
				const int point = list_points_north_(i);
				const double area = (grid_x_.dxc()(i)*0.50)*(grid_y_.dxw()(ny_ - 1)*0.50);
				sum += v(point)*area;
			}

			// West (zero gradient)
			for (unsigned int i = 1; i < ny_ - 1; i++)
			{
				const int point = list_points_west_(i);
				const double area = (grid_x_.dxe()(0)*0.50)*(grid_y_.dxc()(i)*0.50);
				sum += v(point)*area;
			}

			// East (gas side)
			for (unsigned int i = 1; i < ny_ - 1; i++)
			{
				const int point = list_points_east_(i);
				const double area = (grid_x_.dxw()(nx_ - 1)*0.50)*(grid_y_.dxc()(i)*0.50);
				sum += v(point)*area;
			}

			// Corner: north/east
			{
				const int point = list_points_north_(nx_ - 1);
				const double area = (grid_x_.dxw()(nx_ - 1)*0.50)*(grid_y_.dxw()(ny_ - 1)*0.50);
				sum += v(point)*area;
			}

			// Corner: north/west
			{
				const int point = list_points_north_(0);
				const double area = (grid_x_.dxe()(0)*0.50)*(grid_y_.dxw()(ny_ - 1)*0.50);
				sum += v(point)*area;
			}

			// Corner: south/west
			{
				const int point = list_points_south_(0);
				const double area = (grid_x_.dxe()(0)*0.50)*(grid_y_.dxe()(0)*0.50);
				sum += v(point)*area;
			}

			// Corner: south/east
			{
				const int point = list_points_south_(nx_ - 1);
				const double area = (grid_x_.dxw()(nx_ - 1)*0.50)*(grid_y_.dxe()(0)*0.50);
				sum += v(point)*area;
			}

			return sum;
		}
		else
		{
			double sum = 0.;

			// Internal
			for (unsigned int k = 1; k < ny_ - 1; k++)
				for (unsigned int i = 1; i < nx_ - 1; i++)
				{
					const double ri = (grid_x_.x()(i) + grid_x_.x()(i - 1)) / 2.;
					const double re = (grid_x_.x()(i) + grid_x_.x()(i + 1)) / 2.;
					const double height = grid_y_.dxc()(k)*0.50;
					const double volume = boost::math::constants::pi<double>()*(re*re - ri*ri)*height;

					const int point = k*nx_ + i;
					sum += v(point)*volume;
				}

			// South (zero gradient)
			for (unsigned int i = 1; i < nx_ - 1; i++)
			{
				const double ri = (grid_x_.x()(i) + grid_x_.x()(i - 1)) / 2.;
				const double re = (grid_x_.x()(i) + grid_x_.x()(i + 1)) / 2.;
				const double height = grid_y_.dxe()(0)*0.50;
				const double volume = boost::math::constants::pi<double>()*(re*re - ri*ri)*height;

				const int point = list_points_south_(i);
				sum += v(point)*volume;
			}

			// North (zero gradient)
			for (unsigned int i = 1; i < nx_ - 1; i++)
			{
				const double ri = (grid_x_.x()(i) + grid_x_.x()(i - 1)) / 2.;
				const double re = (grid_x_.x()(i) + grid_x_.x()(i + 1)) / 2.;
				const double height = grid_y_.dxw()(ny_ - 1)*0.50;
				const double volume = boost::math::constants::pi<double>()*(re*re - ri*ri)*height;

				const int point = list_points_north_(i);
				sum += v(point)*volume;
			}

			// West (zero gradient)
			for (unsigned int i = 1; i < ny_ - 1; i++)
			{
				const double ri = (grid_x_.x()(nx_ - 1) + grid_x_.x()(nx_ - 2)) / 2.;
				const double re = grid_x_.x()(nx_ - 1);
				const double height = grid_y_.dxc()(i)*0.50;
				const double volume = boost::math::constants::pi<double>()*(re*re - ri*ri)*height;

				const int point = list_points_west_(i);
				sum += v(point)*volume;
			}

			// East (gas side)
			for (unsigned int i = 1; i < ny_ - 1; i++)
			{
				const double ri = grid_x_.x()(0);
				const double re = (grid_x_.x()(0) + grid_x_.x()(1)) / 2.;
				const double height = grid_y_.dxc()(i)*0.50;
				const double volume = boost::math::constants::pi<double>()*(re*re - ri*ri)*height;

				const int point = list_points_east_(i);
				sum += v(point)*volume;
			}

			// Corner: north/west
			{
				const double ri = (grid_x_.x()(nx_ - 1) + grid_x_.x()(nx_ - 2)) / 2.;
				const double re = grid_x_.x()(nx_ - 1);
				const double height = grid_y_.dxw()(ny_ - 1)*0.50;
				const double volume = boost::math::constants::pi<double>()*(re*re - ri*ri)*height;

				const int point = list_points_north_(nx_ - 1);
				sum += v(point)*volume;
			}

			// Corner: north/east
			{
				const double ri = grid_x_.x()(0);
				const double re = (grid_x_.x()(0) + grid_x_.x()(1)) / 2.;
				const double height = grid_y_.dxw()(ny_ - 1)*0.50;
				const double volume = boost::math::constants::pi<double>()*(re*re - ri*ri)*height;

				const int point = list_points_north_(0);
				sum += v(point)*volume;
			}

			// Corner: south/west
			{
				const double ri = (grid_x_.x()(nx_ - 1) + grid_x_.x()(nx_ - 2)) / 2.;
				const double re = grid_x_.x()(nx_ - 1);
				const double height = grid_y_.dxe()(0)*0.50;
				const double volume = boost::math::constants::pi<double>()*(re*re - ri*ri)*height;

				const int point = list_points_south_(0);
				sum += v(point)*volume;
			}

			// Corner: south/east
			{
				const double ri = grid_x_.x()(0);
				const double re = (grid_x_.x()(0) + grid_x_.x()(1)) / 2.;
				const double height = grid_y_.dxe()(0)*0.50;
				const double volume = boost::math::constants::pi<double>()*(re*re - ri*ri)*height;

				const int point = list_points_south_(nx_ - 1);
				sum += v(point)*volume;
			}

			return sum;
		}
	}

	double Reactor2D::VolumeStandardDeviation(const double mean, const Eigen::VectorXd& v)
	{
		if (planar_symmetry_ == true)
		{
			double sum_std = 0.;

			// Internal
			for (unsigned int k = 1; k < ny_ - 1; k++)
				for (unsigned int i = 1; i < nx_ - 1; i++)
				{
					const int point = k*nx_ + i;

					const double area = (grid_x_.dxc()(i)*0.50)*(grid_y_.dxc()(k)*0.50);
					sum_std += boost::math::pow<2>(v(point) - mean)*area;
				}

			// South (zero gradient)
			for (unsigned int i = 1; i < nx_ - 1; i++)
			{
				const int point = list_points_south_(i);
				const double area = (grid_x_.dxc()(i)*0.50)*(grid_y_.dxe()(0)*0.50);
				sum_std += boost::math::pow<2>(v(point) - mean)*area;
			}

			// North (zero gradient)
			for (unsigned int i = 1; i < nx_ - 1; i++)
			{
				const int point = list_points_north_(i);
				const double area = (grid_x_.dxc()(i)*0.50)*(grid_y_.dxw()(ny_ - 1)*0.50);
				sum_std += boost::math::pow<2>(v(point) - mean)*area;
			}

			// West (zero gradient)
			for (unsigned int i = 1; i < ny_ - 1; i++)
			{
				const int point = list_points_west_(i);
				const double area = (grid_x_.dxe()(0)*0.50)*(grid_y_.dxc()(i)*0.50);
				sum_std += boost::math::pow<2>(v(point) - mean)*area;
			}

			// East (gas side)
			for (unsigned int i = 1; i < ny_ - 1; i++)
			{
				const int point = list_points_east_(i);
				const double area = (grid_x_.dxw()(nx_ - 1)*0.50)*(grid_y_.dxc()(i)*0.50);
				sum_std += boost::math::pow<2>(v(point) - mean)*area;
			}

			// Corner: north/east
			{
				const int point = list_points_north_(nx_ - 1);
				const double area = (grid_x_.dxw()(nx_ - 1)*0.50)*(grid_y_.dxw()(ny_ - 1)*0.50);
				sum_std += boost::math::pow<2>(v(point) - mean)*area;
			}

			// Corner: north/west
			{
				const int point = list_points_north_(0);
				const double area = (grid_x_.dxe()(0)*0.50)*(grid_y_.dxw()(ny_ - 1)*0.50);
				sum_std += boost::math::pow<2>(v(point) - mean)*area;
			}

			// Corner: south/west
			{
				const int point = list_points_south_(0);
				const double area = (grid_x_.dxe()(0)*0.50)*(grid_y_.dxe()(0)*0.50);
				sum_std += boost::math::pow<2>(v(point) - mean)*area;
			}

			// Corner: south/east
			{
				const int point = list_points_south_(nx_ - 1);
				const double area = (grid_x_.dxw()(nx_ - 1)*0.50)*(grid_y_.dxe()(0)*0.50);
				sum_std += boost::math::pow<2>(v(point) - mean)*area;
			}

			const double coefficient = double(np_ - 1) / double(np_);
			const double volume_total = grid_x_.L()*grid_y_.L();
			return std::sqrt(sum_std / volume_total / coefficient);
		}
		else
		{
			double sum_std = 0.;

			// Internal
			for (unsigned int k = 1; k < ny_ - 1; k++)
				for (unsigned int i = 1; i < nx_ - 1; i++)
				{
					const double ri = (grid_x_.x()(i) + grid_x_.x()(i - 1)) / 2.;
					const double re = (grid_x_.x()(i) + grid_x_.x()(i + 1)) / 2.;
					const double height = grid_y_.dxc()(k)*0.50;
					const double volume = boost::math::constants::pi<double>()*(re*re - ri*ri)*height;

					const int point = k*nx_ + i;
					sum_std += boost::math::pow<2>(v(point) - mean)*volume;
				}

			// South (zero gradient)
			for (unsigned int i = 1; i < nx_ - 1; i++)
			{
				const double ri = (grid_x_.x()(i) + grid_x_.x()(i - 1)) / 2.;
				const double re = (grid_x_.x()(i) + grid_x_.x()(i + 1)) / 2.;
				const double height = grid_y_.dxe()(0)*0.50;
				const double volume = boost::math::constants::pi<double>()*(re*re - ri*ri)*height;

				const int point = list_points_south_(i);
				sum_std += boost::math::pow<2>(v(point) - mean)*volume;
			}

			// North (zero gradient)
			for (unsigned int i = 1; i < nx_ - 1; i++)
			{
				const double ri = (grid_x_.x()(i) + grid_x_.x()(i - 1)) / 2.;
				const double re = (grid_x_.x()(i) + grid_x_.x()(i + 1)) / 2.;
				const double height = grid_y_.dxw()(ny_ - 1)*0.50;
				const double volume = boost::math::constants::pi<double>()*(re*re - ri*ri)*height;

				const int point = list_points_north_(i);
				sum_std += boost::math::pow<2>(v(point) - mean)*volume;
			}

			// West (zero gradient)
			for (unsigned int i = 1; i < ny_ - 1; i++)
			{
				const double ri = (grid_x_.x()(nx_ - 1) + grid_x_.x()(nx_ - 2)) / 2.;
				const double re = grid_x_.x()(nx_ - 1);
				const double height = grid_y_.dxc()(i)*0.50;
				const double volume = boost::math::constants::pi<double>()*(re*re - ri*ri)*height;

				const int point = list_points_west_(i);
				sum_std += boost::math::pow<2>(v(point) - mean)*volume;
			}

			// East (gas side)
			for (unsigned int i = 1; i < ny_ - 1; i++)
			{
				const double ri = grid_x_.x()(0);
				const double re = (grid_x_.x()(0) + grid_x_.x()(1)) / 2.;
				const double height = grid_y_.dxc()(i)*0.50;
				const double volume = boost::math::constants::pi<double>()*(re*re - ri*ri)*height;

				const int point = list_points_east_(i);
				sum_std += boost::math::pow<2>(v(point) - mean)*volume;
			}

			// Corner: north/west
			{
				const double ri = (grid_x_.x()(nx_ - 1) + grid_x_.x()(nx_ - 2)) / 2.;
				const double re = grid_x_.x()(nx_ - 1);
				const double height = grid_y_.dxw()(ny_ - 1)*0.50;
				const double volume = boost::math::constants::pi<double>()*(re*re - ri*ri)*height;

				const int point = list_points_north_(nx_ - 1);
				sum_std += boost::math::pow<2>(v(point) - mean)*volume;
			}

			// Corner: north/east
			{
				const double ri = grid_x_.x()(0);
				const double re = (grid_x_.x()(0) + grid_x_.x()(1)) / 2.;
				const double height = grid_y_.dxw()(ny_ - 1)*0.50;
				const double volume = boost::math::constants::pi<double>()*(re*re - ri*ri)*height;

				const int point = list_points_north_(0);
				sum_std += boost::math::pow<2>(v(point) - mean)*volume;
			}

			// Corner: south/west
			{
				const double ri = (grid_x_.x()(nx_ - 1) + grid_x_.x()(nx_ - 2)) / 2.;
				const double re = grid_x_.x()(nx_ - 1);
				const double height = grid_y_.dxe()(0)*0.50;
				const double volume = boost::math::constants::pi<double>()*(re*re - ri*ri)*height;

				const int point = list_points_south_(0);
				sum_std += boost::math::pow<2>(v(point) - mean)*volume;
			}

			// Corner: south/east
			{
				const double ri = grid_x_.x()(0);
				const double re = (grid_x_.x()(0) + grid_x_.x()(1)) / 2.;
				const double height = grid_y_.dxe()(0)*0.50;
				const double volume = boost::math::constants::pi<double>()*(re*re - ri*ri)*height;

				const int point = list_points_south_(nx_ - 1);
				sum_std += boost::math::pow<2>(v(point) - mean)*volume;
			}

			const double Ri = grid_x_.x()(0);
			const double Re = grid_x_.x()(nx_ - 1);
			const double volume_total = boost::math::constants::pi<double>()*(Re*Re - Ri*Ri)*grid_y_.L();
			const double coefficient = double(np_ - 1) / double(np_);
			return std::sqrt(sum_std / volume_total / coefficient);
		}
	}

	// Print XML Files
	void Reactor2D::PrintXMLFile(const std::string file_name, const double t)
	{
		const unsigned int n_additional = 5;

		std::ofstream fXML;
		fXML.open(file_name.c_str(), std::ios::out);
		fXML.setf(std::ios::scientific);

		fXML << "<?xml version=\"1.0\" encoding=\"utf-8\"?>" << std::endl;
		fXML << "<opensmoke version=\"0.1a\">" << std::endl;

		fXML << "<Type> CVI::Reactor2D </Type>" << std::endl;

		fXML << "<additional>" << std::endl;
		fXML << n_additional << std::endl;
		fXML << "temperature [K] 2" << std::endl;
		fXML << "pressure [Pa] 3" << std::endl;
		fXML << "mol-weight [kg/kmol] 4" << std::endl;
		fXML << "density [kg/m3] 5" << std::endl;
		fXML << "epsilon [-] 6" << std::endl;
		fXML << "</additional>" << std::endl;

		fXML << "<t-p-mw>" << std::endl;
		fXML << "0 1 2" << std::endl;
		fXML << "</t-p-mw>" << std::endl;

		fXML << "<time>" << std::endl;
		fXML << t << std::endl;
		fXML << "</time>" << std::endl;

		fXML << "<grid-x>" << std::endl;
		fXML << grid_x_.x().size() << std::endl;
		for (int i = 0; i < grid_x_.x().size(); i++)
			fXML << grid_x_.x()(i) << std::endl;
		fXML << "</grid-x>" << std::endl;

		fXML << "<grid-y>" << std::endl;
		fXML << grid_y_.x().size() << std::endl;
		for (int i = 0; i < grid_y_.x().size(); i++)
			fXML << grid_y_.x()(i) << std::endl;
		fXML << "</grid-y>" << std::endl;

		fXML << "<mass-fractions>" << std::endl;
		fXML << thermodynamicsMap_.NumberOfSpecies() << std::endl;
		for (unsigned int i = 0; i < thermodynamicsMap_.NumberOfSpecies(); i++)
			fXML << thermodynamicsMap_.NamesOfSpecies()[i] << " " << thermodynamicsMap_.MW(i) << " " << n_additional + (i + 1) << std::endl;
		fXML << "</mass-fractions>" << std::endl;

		fXML << "<surface-fractions>" << std::endl;
		fXML << surf_nc_ << std::endl;
		for (unsigned int i = 0; i < surf_nc_; i++)
			fXML << thermodynamicsSurfaceMap_.NamesOfSpecies()[nc_+i] << " " << thermodynamicsSurfaceMap_.MW(nc_+i) << " " << n_additional + nc_ + (i + 1) << std::endl;
		fXML << "</surface-fractions>" << std::endl;

		fXML << "<surface-phases>" << std::endl;
		fXML << surf_np_ << std::endl;
		for (unsigned int i = 0; i < surf_np_; i++)
			fXML << "Gamma" << i << " " << n_additional + nc_ + surf_nc_ + (i + 1) << std::endl;
		fXML << "</surface-phases>" << std::endl;

		fXML << "<profiles>" << std::endl;
		{
			for (unsigned int k = 0; k < ny_; k++)
				for (unsigned int i = 0; i < nx_; i++)
				{
					const int point = k*nx_ + i;

					fXML << T_(point) << " ";
					fXML << P_(point) << " ";
					fXML << mw_(point) << " ";
					fXML << rho_gas_(point) << " ";
					fXML << epsilon_(point) << " ";

					for (unsigned int j = 0; j < nc_; j++)
						fXML << Y_[i](j) << " ";

					for (unsigned int j = 0; j < surf_nc_; j++)
						fXML << Z_[i](j) << " ";

					for (unsigned int j = 0; j < surf_np_; j++)
						fXML << GammaFromEqn_[i](j) << " ";

					fXML << std::endl;
				}
		}
		fXML << "</profiles>" << std::endl;

		fXML << "<profiles-size>" << std::endl;
		fXML << np_ << " " << nc_ + surf_nc_ + surf_np_ + n_additional << std::endl;
		fXML << "</profiles-size>" << std::endl;
		fXML << "</opensmoke>" << std::endl;
		fXML.close();
	}

	void Reactor2D::SetInitialConditionsFromBackupFile(const boost::filesystem::path path_to_backup_file)
	{
		std::cout << "--------------------------------------------------------------------------" << std::endl;
		std::cout << "Reading initial conditions from backup file..." << std::endl;
		std::cout << "--------------------------------------------------------------------------" << std::endl;

		ReadFromBackupFile(	path_to_backup_file, time_starting_point_, 
							thermodynamicsSurfaceMap_,
							grid_x_.x(), grid_y_.x(),
							T_, P_, epsilon_,
							Y_, Z_, GammaFromEqn_ );

		for (unsigned int i = 0; i < np_; i++)
			for (unsigned int j = 0; j < surf_np_; j++)
				if (site_non_conservation_[j] == true)
					Gamma_[i](j) = GammaFromEqn_[i](j);
	}

}
