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
	Reactor2D::Reactor2D(	OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>& thermodynamicsMap,
							OpenSMOKE::KineticsMap_CHEMKIN<double>& kineticsMap,
							OpenSMOKE::TransportPropertiesMap_CHEMKIN<double>& transportMap,
							CVI::PorousMedium& porousMedium,
							CVI::PorosityDefect& porosityDefect,
							OpenSMOKE::Grid1D& grid_x, OpenSMOKE::Grid1D& grid_y,
							CVI::PlugFlowReactor& plugFlowReactor) :

	thermodynamicsMap_(thermodynamicsMap),
	kineticsMap_(kineticsMap),
	transportMap_(transportMap),
	porousMedium_(porousMedium),
	porosityDefect_(porosityDefect),
	grid_x_(grid_x),
	grid_y_(grid_y),
	plugFlowReactor_(plugFlowReactor)
	{
		t_old_ = 0.;
		time_smoothing_ = 10.;

		n_steps_video_ = 10;
		count_video_ = n_steps_video_;
		n_steps_file_ = 3;
		count_file_ = n_steps_file_;
		count_tecplot_ = 0;

		time_total_ = 48.*3600.;
		dae_time_interval_ = 3600.;
		tecplot_time_interval_ = 3600.;

		ns_ = thermodynamicsMap_.NumberOfSpecies();
		block_ = ns_ + 1;
		nx_ = grid_x_.Np();
		ny_ = grid_y_.Np();
		np_ = nx_*ny_;
		ne_ = block_*np_;
		band_size_ = block_*(nx_+1)-1;

		planar_symmetry_ = true;

		output_folder_ = "Output";

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
		OpenSMOKE::ChangeDimensions(ns_, &aux_X, true);
		OpenSMOKE::ChangeDimensions(ns_, &aux_Y, true);
		OpenSMOKE::ChangeDimensions(ns_, &aux_C, true);
		OpenSMOKE::ChangeDimensions(ns_, &aux_R, true);
		aux_eigen.resize(ns_);
		
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
		for (int i = 0; i < np_; i++)
			X_[i].resize(ns_);

		// Mass fractions [-]
		Y_.resize(np_);
		for (int i = 0; i < np_; i++)
			Y_[i].resize(ns_);

		// Formation rates in gas pahse [kg/m3/s]
		omega_homogeneous_.resize(np_);
		for (int i = 0; i < np_; i++)
			omega_homogeneous_[i].resize(ns_);

		// Heterogeneous formation rates [kg/m3/s]
		omega_heterogeneous_.resize(np_);
		for (int i = 0; i < np_; i++)
			omega_heterogeneous_[i].resize(ns_);

		// Deposition rate [kg/m3/s]
		omega_deposition_per_unit_volume_.resize(np_);

		// Deposition rate [kg/m2/s]
		omega_deposition_per_unit_area_.resize(np_);

		// Effective diffusion coefficiennts [m2/s]
		gamma_star_.resize(np_);
		for (int i = 0; i < np_; i++)
			gamma_star_[i].resize(ns_);

		// Time derivatives: mass fractions [1/s]
		dY_over_dt_.resize(np_);
		for (int i = 0; i < np_; i++)
			dY_over_dt_[i].resize(ns_);

		// Time derivatives: porosity [1/s]
		depsilon_over_dt_.resize(np_);


		// Gas side Mass fractions [-]
		Y_gas_side_.resize(2*nx_+ny_);
		for (int i = 0; i < 2 * nx_ + ny_; i++)
			Y_gas_side_[i].resize(ns_);

		// Gas side temperature and pressure
		T_gas_side_.resize(2 * nx_ + ny_);
		P_gas_side_.resize(2 * nx_ + ny_);

		// Increment of bulk density due to the single reactions
		delta_rhobulk_.resize(np_);
		delta_rhobulk_due_to_single_reaction_.resize(porousMedium_.r().size());
		delta_rhobulk_due_to_single_reaction_over_rhobulk_.resize(porousMedium_.r().size());
		for (int i = 0; i < porousMedium_.r().size(); i++)
		{
			delta_rhobulk_due_to_single_reaction_[i].resize(np_);
			delta_rhobulk_due_to_single_reaction_over_rhobulk_[i].resize(np_);
			delta_rhobulk_due_to_single_reaction_[i].setZero();
			delta_rhobulk_due_to_single_reaction_over_rhobulk_[i].setZero();
		}

		// Set to zero relevant fields
		for (int i = 0; i < np_; i++)
			omega_homogeneous_[i].setZero();
		for (int i = 0; i < np_; i++)
			omega_heterogeneous_[i].setZero();
		omega_deposition_per_unit_volume_.setZero();
		omega_deposition_per_unit_area_.setZero();
	}

	void Reactor2D::SetPlanarSymmetry(const bool flag)
	{
		planar_symmetry_ = flag;
	}

	void Reactor2D::SetAlgebraicDifferentialEquations()
	{
		id_equations_.resize(ne_);

		unsigned int count = 0;

		// South side
		for (unsigned int i = 0; i < nx_; i++)
		{
			// Species
			for (unsigned int j = 0; j < ns_; j++)
				id_equations_[count++] = false;

			// Porosity
			id_equations_[count++] = true;
		}
			

		// Internal strips
		for (unsigned int j = 1; j < ny_ - 1; j++)
		{
			// West side
			{
				// Species
				for (unsigned int j = 0; j < ns_; j++)
					id_equations_[count++] = false;

				// Porosity
				id_equations_[count++] = true;
			}

			// Internal point
			for (unsigned int i = 1; i < nx_-1; i++)
			{
				// Species
				for (unsigned int j = 0; j < ns_; j++)
					id_equations_[count++] = true;

				// Porosity
				id_equations_[count++] = true;
			}

			// East side (gas side)
			{
				// Species
				for (unsigned int j = 0; j < ns_; j++)
					id_equations_[count++] = false;

				// Porosity
				id_equations_[count++] = true;
			}
		}

		// North side
		for (unsigned int i = 0; i < nx_; i++)
		{
			// Species
			for (unsigned int j = 0; j < ns_; j++)
				id_equations_[count++] = false;

			// Porosity
			id_equations_[count++] = true;
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
		// Set gas side variables
		for (unsigned int i = 0; i < 2 * nx_ + ny_; i++)
		{
			for (unsigned int j = 0; j < ns_; j++)
				Y_gas_side_[i](j) = omega_gas[i](j);

			P_gas_side_(i) = P_gas;
			T_gas_side_(i) = T_gas;
		}

		// Set initial fields consistent along the east side
		if (plugFlowReactor_.geometric_pattern() == CVI::PlugFlowReactor::GeometricPattern::ONE_SIDE)
		{
			for (unsigned int i = 0; i < ny_; i++)
			{
				for (unsigned int j = 0; j < ns_; j++)
					Y_[list_points_east_(i)](j) = omega_gas[nx_+i](j);

				P_(list_points_east_(i)) = P_gas;
				T_(list_points_east_(i)) = T_gas;
			}
		}
		// Set initial fields consistent along the south, east, and north sides
		else if (plugFlowReactor_.geometric_pattern() == CVI::PlugFlowReactor::GeometricPattern::THREE_SIDES)
		{
			for (unsigned int i = 0; i < nx_; i++)
			{
				for (unsigned int j = 0; j < ns_; j++)
					Y_[list_points_south_(i)](j) = omega_gas[i](j);

				P_(list_points_south_(i)) = P_gas;
				T_(list_points_south_(i)) = T_gas;
			}

			for (unsigned int i = 0; i < ny_; i++)
			{
				for (unsigned int j = 0; j < ns_; j++)
					Y_[list_points_east_(i)](j) = omega_gas[nx_ + i](j);

				P_(list_points_east_(i)) = P_gas;
				T_(list_points_east_(i)) = T_gas;
			}

			for (unsigned int i = 0; i < nx_; i++)
			{
				for (unsigned int j = 0; j < ns_; j++)
					Y_[list_points_north_(i)](j) = omega_gas[2*nx_+ny_-1-i](j);

				P_(list_points_north_(i)) = P_gas;
				T_(list_points_north_(i)) = T_gas;
			}
		}
	}

	void Reactor2D::SetInitialConditions(const double T_gas, const double P_gas, const Eigen::VectorXd& omega_gas)
	{
		for (int i = 0; i < np_; i++)
			for (unsigned int j = 0; j < ns_; j++)
				Y_[i](j) = omega_gas(j);

		P_.setConstant(P_gas);
		T_.setConstant(T_gas);
		epsilon_.setConstant(porousMedium_.porosity());

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
	}

	void Reactor2D::SetTimeTotal(const double time_total)
	{
		time_total_ = time_total;
	}

	void Reactor2D::SetDaeTimeInterval(const double time_interval)
	{
		dae_time_interval_ = time_interval;
	}

	void Reactor2D::SetTecplotTimeInterval(const double time_interval)
	{
		tecplot_time_interval_ = time_interval;
	}

	void Reactor2D::SetStepsVideo(const int steps_video)
	{
		n_steps_video_ = steps_video;
		count_video_ = n_steps_video_;

	}

	void Reactor2D::SetStepsFile(const int steps_file)
	{
		n_steps_file_ = steps_file;
		count_file_ = n_steps_file_;
	}

	void Reactor2D::Properties()
	{
		for (int i = 0; i < np_; i++)
		{
			// Thermodynamics
			{
				thermodynamicsMap_.SetPressure(P_(i));
				thermodynamicsMap_.SetTemperature(T_(i));

				aux_Y.CopyFrom(Y_[i].data());
				thermodynamicsMap_.MoleFractions_From_MassFractions(aux_X, mw_(i), aux_Y);
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
				rho_bulk_(i) = porousMedium_.density_bulk();
				eta_bulk_(i) = porousMedium_.eta_bulk();
				eta_knudsen_(i) = porousMedium_.eta_knudsen();
				eta_viscous_(i) = porousMedium_.eta_viscous();

				// Mixture diffusion coefficients
				{
					porousMedium_.EffectiveDiffusionCoefficients(X_[i]);
					for (unsigned int j = 0; j < ns_; j++)
						gamma_star_[i](j) = porousMedium_.gamma_effective()(j);
				}
			}

			// Kinetics
			{
				// Homogeneous phase
				if (porousMedium_.homogeneous_reactions() == true)
				{
					kineticsMap_.SetTemperature(T_(i));
					kineticsMap_.SetPressure(P_(i));
					kineticsMap_.ReactionRates(aux_C);
					kineticsMap_.FormationRates(&aux_R);
					ElementByElementProduct(aux_R, thermodynamicsMap_.MW(), &aux_R); // [kg/m3/s]
					aux_R.CopyTo(omega_homogeneous_[i].data());
				}

				// Heterogeneous phase
				if (porousMedium_.heterogeneous_reactions() == true)
				{
					aux_C.CopyTo(aux_eigen.data());
					porousMedium_.FormationRates(aux_eigen);
					for (unsigned int j = 0; j < ns_; j++)
						omega_heterogeneous_[i](j) = porousMedium_.Rgas()(j)*thermodynamicsMap_.MW()[j+1];								// [kg/m3/s]
					omega_deposition_per_unit_area_(i) = porousMedium_.r_deposition_per_unit_area()*porousMedium_.mw_carbon();	// [kg/m2/s]
					omega_deposition_per_unit_volume_(i) = porousMedium_.r_deposition_per_unit_volume()*porousMedium_.mw_carbon();		// [kg/m3/s]
			}
			}
		}
	}

	void Reactor2D::SubEquations_MassFractions_BoundaryConditions_WestSide(const double t)
	{
		for (unsigned int i = 0; i < ny_; i++)
		{
			const int point = list_points_west_(i);
			for (unsigned int j = 0; j < ns_; j++)
				dY_over_dt_[point](j) = Y_[point](j) - Y_[point + 1](j);
		}
	}

	void Reactor2D::SubEquations_MassFractions_BoundaryConditions_EastSide(const double t)
	{
		if (plugFlowReactor_.internal_boundary_layer_correction() == true && t > time_smoothing_)
		{
			for (unsigned int i = 0; i < ny_; i++)
			{
				const int point = list_points_east_(i);

				const double kc = plugFlowReactor_.mass_transfer_coefficient(T_(point), P_(point), rho_gas_(point), grid_y_.x()(i)+plugFlowReactor_.inert_length());
				const double gamma_over_dx = gamma_star_[point](porousMedium_.index_CH4()) / grid_x_.dxw()(nx_ - 1);

				for (unsigned int j = 0; j < ns_; j++)
					dY_over_dt_[point](j) = Y_[point](j) - (gamma_over_dx*Y_[point - 1](j) + kc*Y_gas_side_[nx_ + i](j)) / (gamma_over_dx + kc);
			}
		}
		else
		{
			for (unsigned int i = 0; i < ny_; i++)
			{
				const int point = list_points_east_(i);
				for (unsigned int j = 0; j < ns_; j++)
					dY_over_dt_[point](j) = Y_[point](j) - Y_gas_side_[nx_ + i](j);
			}
		}
	}

	void Reactor2D::SubEquations_MassFractions_BoundaryConditions_NorthSide(const double t)
	{
		if (plugFlowReactor_.geometric_pattern() == CVI::PlugFlowReactor::ONE_SIDE)
		{
			for (unsigned int i = 0; i < nx_; i++)
			{
				const int point = list_points_north_(i);
				for (unsigned int j = 0; j < ns_; j++)
					dY_over_dt_[point](j) = Y_[point](j) - Y_[point - nx_](j);
			}
		}
		else if (plugFlowReactor_.geometric_pattern() == CVI::PlugFlowReactor::THREE_SIDES)
		{
			for (unsigned int i = 0; i < nx_; i++)
			{
				const int point = list_points_north_(i);
				for (unsigned int j = 0; j < ns_; j++)
					dY_over_dt_[point](j) = Y_[point](j) - Y_gas_side_[2 * nx_ + ny_ - 1 - i](j);
			}
		}
	}

	void Reactor2D::SubEquations_MassFractions_BoundaryConditions_SouthSide(const double t)
	{
		if (plugFlowReactor_.geometric_pattern() == CVI::PlugFlowReactor::ONE_SIDE)
		{
			for (unsigned int i = 0; i < nx_; i++)
			{
				const int point = list_points_south_(i);
				for (unsigned int j = 0; j < ns_; j++)
					dY_over_dt_[point](j) = Y_[point](j) - Y_[point + nx_](j);
			}
		}
		else if (plugFlowReactor_.geometric_pattern() == CVI::PlugFlowReactor::THREE_SIDES)
		{
			for (unsigned int i = 0; i < nx_; i++)
			{
				const int point = list_points_south_(i);
				for (unsigned int j = 0; j < ns_; j++)
					dY_over_dt_[point](j) = Y_[point](j) - Y_gas_side_[i](j);
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

				for (unsigned int j = 0; j < ns_; j++)
				{
					const double c_east = 0.50* (gamma_star_[east](j)*rho_gas_[east] + gamma_star_[center](j)*rho_gas_[center]);
					const double c_west = 0.50* (gamma_star_[west](j)*rho_gas_[west] + gamma_star_[center](j)*rho_gas_[center]);
					const double c_north = 0.50* (gamma_star_[north](j)*rho_gas_[north] + gamma_star_[center](j)*rho_gas_[center]);
					const double c_south = 0.50* (gamma_star_[south](j)*rho_gas_[south] + gamma_star_[center](j)*rho_gas_[center]);

					      double diffusion_x = (c_east*(Y_[east](j) - Y_[center](j)) / grid_x_.dxe()(i)-c_west*(Y_[center](j) - Y_[west](j)) / grid_x_.dxw()(i)) / grid_x_.dxc_over_2()(i);
					const double diffusion_y = (c_north*(Y_[north](j) - Y_[center](j)) / grid_y_.dxe()(k)-c_south*(Y_[center](j) - Y_[south](j)) / grid_y_.dxw()(k)) / grid_y_.dxc_over_2()(k);
					const double diffusion = diffusion_x + diffusion_y;

					const double homogeneous_reactions = epsilon_(center)*omega_homogeneous_[center](j);
					const double heterogeneous_reactions = omega_heterogeneous_[center](j) + Y_[center](j)*omega_deposition_per_unit_volume_(center);

					if (planar_symmetry_ == false)
					{
						const double dY_over_dr = (Y_[east](j) - Y_[west](j)) / (grid_x_.x()[i + 1] - grid_x_.x()[i - 1]);
						const double diffusion_radial = gamma_star_[center](j)*rho_gas_[center] * dY_over_dr / grid_x_.x()[i];
						diffusion_x += diffusion_radial;
					}

					dY_over_dt_[center](j) = diffusion + homogeneous_reactions + heterogeneous_reactions;
					dY_over_dt_[center](j) /= (rho_gas_(center)*epsilon_(center));
				}
			}
	}

	void Reactor2D::SubEquations_Porosity()
	{
		for (int i = 0; i < np_; i++)
			depsilon_over_dt_(i) = -omega_deposition_per_unit_volume_(i) / porousMedium_.rho_graphite();
	}

	void Reactor2D::Recover_Unknowns(const double* y)
	{
		unsigned int count = 0;
		for (int i = 0; i < np_; i++)
		{
			// Species
			for (unsigned int j = 0; j < ns_; j++)
				Y_[i](j) = y[count++];

			// Porosity
			epsilon_(i) = y[count++];
		}
	}

	void Reactor2D::Recover_Residuals(double* dy)
	{
		unsigned int count = 0;
		for (int i = 0; i < np_; i++)
		{
			// Species
			for (unsigned int j = 0; j < ns_; j++)
				dy[count++] = dY_over_dt_[i](j);
			
			// Porosity
			dy[count++] = depsilon_over_dt_(i);
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
		unsigned int count = 0;
		for (int i = 0; i < np_; i++)
		{
			for (unsigned int j = 0; j < ns_; j++)
				v[count++] = Y_[i](j);

			v[count++] = epsilon_[i];
		}
	}

	void Reactor2D::CorrectedUnknownsVector(double* v)
	{
		unsigned int count = 0;
		for (int i = 0; i < np_; i++)
		{
			for (unsigned int j = 0; j < ns_; j++)
				Y_[i](j) = v[count++];

			epsilon_(i) = v[count++];

			// Normalize
			const double sum = Y_[i].sum();
			for (unsigned int j = 0; j < ns_; j++)
				Y_[i](j) /= sum;
		}

		Properties();
	}

	void Reactor2D::MinimumUnknownsVector(double* v)
	{
		const double zero = 0.;

		unsigned int count = 0;
		for (int i = 0; i < np_; i++)
		{
			for (unsigned int j = 0; j < ns_; j++)
				v[count++] = zero;

			v[count++] = zero;
		}
	}

	void Reactor2D::MaximumUnknownsVector(double* v)
	{
		const double one = 1.;

		unsigned int count = 0;
		for (int i = 0; i < np_; i++)
		{
			for (unsigned int j = 0; j < ns_; j++)
				v[count++] = one;

			v[count++] = one;
		}
	}

	void Reactor2D::Equations(const double t, const double* y, double* dy)
	{
		// Recover unknowns
		Recover_Unknowns(y);

		// Properties
		Properties();

		// Equations
		SubEquations_MassFractions(t);
		SubEquations_Porosity();

		// Recover residuals
		Recover_Residuals(dy);
	}

	int Reactor2D::SolveFromScratch(DaeSMOKE::DaeSolver_Parameters& dae_parameters)
	{
		// Print initial solution
		Properties();
		PrintTecplot(0., (output_folder_ / "Solution.tec.0").string().c_str());

		// Loop
		unsigned int number_intervals = time_total_/dae_time_interval_;
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
			PrintSolution(tf, (output_folder_ / solution_file).string().c_str());
			PrintDiffusionCoefficients(tf, (output_folder_ / diffusion_coefficients_file).string().c_str());
			PrintHeterogeneousRates(tf, (output_folder_ / heterogeneous_rates_file).string().c_str());
			PrintHomogeneousRates(tf, (output_folder_ / homogeneous_rates_file).string().c_str());
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
			flag = DaeSMOKE::Solve_Band_BzzDae<Reactor2D, OpenSMOKE_Reactor2D_BzzDaeSystem>(this, dae_object_, dae_parameters, t0, tf);
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
		OpenSMOKE::OpenSMOKEVectorDouble yy(ns_);
		OpenSMOKE::OpenSMOKEVectorDouble xx(ns_);

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
			for (unsigned int j = 0; j < ns_; j++)
				fOutput << "\"" << thermodynamicsMap_.NamesOfSpecies()[j] << "\", ";
			
			// Heterogeneous reaction rates
			{
				fOutput << "\"rDep[kg/m2/s]\"" << ", ";
				fOutput << "\"rDep[kg/m3/s]\"" << ", ";
				fOutput << "\"rDep[m/s]\"" << ", ";
				for (int j = 0; j < porousMedium_.r().size(); j++)
				{
					std::stringstream number; number << j + 1;
					fOutput << "\"" << ("rDep" + number.str() + "[kg/m2/s]") << "\", ";
					fOutput << "\"" << ("rDep" + number.str() + "[m/s]") << "\", ";
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
				for (int j = 0; j < porousMedium_.r().size(); j++)
				{
					std::stringstream number; number << j + 1;
					fOutput << "\"" << ("dRhoBulk" + number.str() + "[kg/m3]") << "\", ";
					fOutput << "\"" << ("dRhoBulk" + number.str()) << "\", ";
				}
			}
			
			// Finalize
			fOutput << std::endl;

			// Titles
			fOutput << "Zone I=" << nx_ << ", J=" << ny_ << ", F=POINT" << std::endl;
		}

		for (int k = 0; k < ny_; k++)
			for (int i = 0; i < nx_; i++)
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
				thermodynamicsMap_.MoleFractions_From_MassFractions(aux_X, mw_(point), aux_Y);
				aux_X.CopyTo(X_[point].data());
				const double cTot = P_(point) / PhysicalConstants::R_J_kmol / T_(point); // [kmol/m3]
				Product(cTot, aux_X, &aux_C);

				// Heterogeneous reactions
				for (unsigned int j = 0; j < ns_; j++)
					aux_eigen(j) = aux_C[j + 1];
				porousMedium_.FormationRates(aux_eigen);


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
				for (unsigned int j = 0; j < ns_; j++)
					fOutput << std::setprecision(9) << std::setw(20) << Y_[point](j);

				// Heterogeneous reaction rates
				{
					fOutput << std::setprecision(9) << std::setw(20) << porousMedium_.r_deposition_per_unit_area()*porousMedium_.mw_carbon();		// [kg/m2/s]
					fOutput << std::setprecision(9) << std::setw(20) << porousMedium_.r_deposition_per_unit_volume()*porousMedium_.mw_carbon();		// [kg/m3/s]
					fOutput << std::setprecision(9) << std::setw(20) << porousMedium_.r_deposition_per_unit_volume()*porousMedium_.mw_carbon() /
						(porousMedium_.Sv() / porousMedium_.rho_graphite());																		// [m/s]
					for (int j = 0; j < porousMedium_.r().size(); j++)
					{
						fOutput << std::setprecision(9) << std::setw(20) << porousMedium_.r()(j)*porousMedium_.mw_carbon();									// [kg/m2/s]
						fOutput << std::setprecision(9) << std::setw(20) << porousMedium_.r()(j)*porousMedium_.mw_carbon() / porousMedium_.rho_graphite();	// [m/s]
					}
				}

				// Heterogeneous formation rates
				{
					fOutput << std::setprecision(9) << std::setw(20) << -porousMedium_.Rgas()(porousMedium_.index_CH4());
					fOutput << std::setprecision(9) << std::setw(20) << -porousMedium_.Rgas()(porousMedium_.index_C2H4());
					fOutput << std::setprecision(9) << std::setw(20) << -porousMedium_.Rgas()(porousMedium_.index_C2H2());
					fOutput << std::setprecision(9) << std::setw(20) << -porousMedium_.Rgas()(porousMedium_.index_C6H6());
					fOutput << std::setprecision(9) << std::setw(20) << -porousMedium_.Rgas()(porousMedium_.index_H2());
				}

				// Hydrogen inhibition factors
				{
					fOutput << std::setprecision(9) << std::setw(20) << porousMedium_.I_CH4();
					fOutput << std::setprecision(9) << std::setw(20) << porousMedium_.I_C2H4();
					fOutput << std::setprecision(9) << std::setw(20) << porousMedium_.I_C2H2();
					fOutput << std::setprecision(9) << std::setw(20) << porousMedium_.I_C6H6();
				}

				// Heterogeneous reaction rates: contributions to the bulk density
				{
					fOutput << std::setprecision(9) << std::setw(20) << delta_rhobulk_(point);	// [kg/m3]
					for (int j = 0; j < porousMedium_.r().size(); j++)
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

	void Reactor2D::PrintSolution(const double t, const std::string name_file)
	{
		std::ofstream fOutput(name_file.c_str(), std::ios::out);

		{
			unsigned int count = 1;
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "time[s]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "x[mm]", count);
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
			for (unsigned int j = 0; j < ns_; j++)
				OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, thermodynamicsMap_.NamesOfSpecies()[j] + "_x", count);
			for (unsigned int j = 0; j < ns_; j++)
				OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, thermodynamicsMap_.NamesOfSpecies()[j] + "_w", count);

			fOutput << std::endl;
		}

		for (int i = 0; i < np_; i++)
		{
			OpenSMOKE::OpenSMOKEVectorDouble yy(ns_);
			OpenSMOKE::OpenSMOKEVectorDouble xx(ns_);

			for (unsigned int j = 0; j < ns_; j++)
				yy[j + 1] = Y_[i](j);

			thermodynamicsMap_.SetPressure(P_(i));
			thermodynamicsMap_.SetTemperature(T_(i));
			thermodynamicsMap_.MoleFractions_From_MassFractions(xx, mw_(i), yy);

			fOutput << std::setprecision(9) << std::setw(20) << t;
			fOutput << std::setprecision(9) << std::setw(20) << i;//grid_.x()[i] * 1000.;
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
					for (unsigned int k = 0; k < ns_; k++)
						sum += thermodynamicsMap_.atomic_composition()(k, j) * xx(k + 1);
					fOutput << std::setprecision(9) << std::setw(20) << sum;
				}
			}

			// Species (mole fractions and mass fractions)
			{
				for (unsigned int j = 0; j < ns_; j++)
					fOutput << std::setprecision(9) << std::setw(20) << xx(j + 1);
				for (unsigned int j = 0; j < ns_; j++)
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

			for (unsigned int j = 0; j < ns_; j++)
			{
				OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, thermodynamicsMap_.NamesOfSpecies()[j] + "_ToE", count);
				OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, thermodynamicsMap_.NamesOfSpecies()[j] + "_FiE", count);
				OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, thermodynamicsMap_.NamesOfSpecies()[j] + "_KnE", count);
				OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, thermodynamicsMap_.NamesOfSpecies()[j] + "_Fi", count);
				OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, thermodynamicsMap_.NamesOfSpecies()[j] + "_Kn", count);
			}

			fOutput << std::endl;
		}

		for (int i = 0; i < np_; i++)
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
				for (unsigned int j = 0; j < ns_; j++)
					aux_eigen(j) = aux_X[j + 1];
				porousMedium_.EffectiveDiffusionCoefficients(aux_eigen);
			}

			fOutput << std::setprecision(9) << std::setw(20) << t;
			fOutput << std::setprecision(9) << std::setw(20) << i;// grid_.x()[i] * 1000.;
			fOutput << std::setprecision(9) << std::setw(20) << T_(i);
			fOutput << std::setprecision(9) << std::setw(20) << P_(i);

			// Diffusion coefficients
			{
				for (unsigned int j = 0; j < ns_; j++)
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

			for (unsigned int j = 0; j < ns_; j++)
				OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, thermodynamicsMap_.NamesOfSpecies()[j] + "[kmol/m3/s]", count);

			fOutput << std::endl;
		}

		for (int i = 0; i < np_; i++)
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
			fOutput << std::setprecision(9) << std::setw(20) << i;// grid_.x()[i] * 1000.;
			fOutput << std::setprecision(9) << std::setw(20) << T_(i);
			fOutput << std::setprecision(9) << std::setw(20) << P_(i);

			for (unsigned int j = 0; j < ns_; j++)
					fOutput << std::setprecision(9) << std::setw(20) << aux_R[j+1];

			fOutput << std::endl;
		}

		fOutput.close();
	}

	void Reactor2D::PrintHeterogeneousRates(const double t, const std::string name_file)
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
			for (int j = 0; j < porousMedium_.r().size(); j++)
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

		for (int i = 0; i < np_; i++)
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

				for (unsigned int j = 0; j < ns_; j++)
					aux_eigen(j) = aux_C[j + 1];
				porousMedium_.FormationRates(aux_eigen);
			}

			fOutput << std::setprecision(9) << std::setw(width) << t;
			fOutput << std::setprecision(9) << std::setw(width) << i;// grid_.x()[i] * 1000.;
			fOutput << std::setprecision(9) << std::setw(width) << T_(i);
			fOutput << std::setprecision(9) << std::setw(width) << P_(i);

			// Reaction rates
			{
				fOutput << std::setprecision(9) << std::setw(width) << porousMedium_.r_deposition_per_unit_area();
				for (int j = 0; j < porousMedium_.r().size(); j++)
					fOutput << std::setprecision(9) << std::setw(width) << porousMedium_.r()(j);
			}

			// Hydrogen inhibition factors
			{
				fOutput << std::setprecision(9) << std::setw(width) << porousMedium_.I_CH4();
				fOutput << std::setprecision(9) << std::setw(width) << porousMedium_.I_C2H4();
				fOutput << std::setprecision(9) << std::setw(width) << porousMedium_.I_C2H2();
				fOutput << std::setprecision(9) << std::setw(width) << porousMedium_.I_C6H6();

			}

			// Heterogeneous formation rates
			{
				fOutput << std::setprecision(9) << std::setw(width) << porousMedium_.Rgas()(porousMedium_.index_CH4());
				fOutput << std::setprecision(9) << std::setw(width) << porousMedium_.Rgas()(porousMedium_.index_C2H4());
				fOutput << std::setprecision(9) << std::setw(width) << porousMedium_.Rgas()(porousMedium_.index_C2H2());
				fOutput << std::setprecision(9) << std::setw(width) << porousMedium_.Rgas()(porousMedium_.index_C6H6());
				fOutput << std::setprecision(9) << std::setw(width) << porousMedium_.Rgas()(porousMedium_.index_H2());
			}

			fOutput << std::endl;
		}

		fOutput.close();
	}

	void Reactor2D::PrintLabelMonitoringFile()
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

	void Reactor2D::Print(const double t, const double* y)
	{
		if (count_video_ == n_steps_video_)
		{
			const double min_epsilon = epsilon_.minCoeff();
			const double max_epsilon = epsilon_.maxCoeff();

			const double rho_bulk_mean = AreaAveraged(rho_bulk_);
			const double rho_bulk_std = AreaStandardDeviation(rho_bulk_mean, rho_bulk_);

			std::cout << std::left << std::setprecision(9) << std::setw(20) << t/3600.;
			std::cout << std::left << std::setprecision(9) << std::setw(20) << t;
			std::cout << std::left << std::setw(16) << std::fixed << std::setprecision(6) << min_epsilon;
			std::cout << std::left << std::setw(16) << std::fixed << std::setprecision(6) << max_epsilon;
			std::cout << std::left << std::setw(16) << std::fixed << std::setprecision(6) << rho_bulk_mean;
			std::cout << std::left << std::setw(16) << std::fixed << std::setprecision(6) << rho_bulk_std;
			std::cout << std::endl;

			count_video_ = 0;
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

			fMonitoring_ << std::left << std::setw(width) << std::scientific << std::setprecision(6) << r_deposition_per_unit_area_mean/porousMedium_.rho_graphite()*1000.;		// [mm/s]
			fMonitoring_ << std::left << std::setw(width) << std::scientific << std::setprecision(6) << r_deposition_per_unit_area_std / porousMedium_.rho_graphite()*1000.;	// [mm/s]

			fMonitoring_ << std::left << std::setw(width) << std::scientific << std::setprecision(6) << r_deposition_per_unit_area_mean *1000.*3600.;	// [g/m2/h]
			fMonitoring_ << std::left << std::setw(width) << std::scientific << std::setprecision(6) << r_deposition_per_unit_area_std *1000.*3600.;		// [g/m2/h]

			fMonitoring_ << std::left << std::setw(width) << std::scientific << std::setprecision(6) << r_deposition_per_unit_volume_mean *1000.*3600.;	// [g/m3/h]
			fMonitoring_ << std::left << std::setw(width) << std::scientific << std::setprecision(6) << r_deposition_per_unit_volume_std *1000.*3600.;	// [g/m3/h]

			fMonitoring_ << std::endl;

			count_file_ = 0;
		}

		// Post-processing
		{
			const double delta_time = t - t_old_;
			const double coefficient = porousMedium_.mw_carbon()*delta_time;

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
				thermodynamicsMap_.MoleFractions_From_MassFractions(aux_X, mw_(i), aux_Y);
				aux_X.CopyTo(X_[i].data());
				const double cTot = P_(i) / PhysicalConstants::R_J_kmol / T_(i); // [kmol/m3]
				Product(cTot, aux_X, &aux_C);

				// Heterogeneous reactions
				for (unsigned int j = 0; j < ns_; j++)
					aux_eigen(j) = aux_C[j + 1];
				porousMedium_.FormationRates(aux_eigen);

				// Single contributions
				delta_rhobulk_(i) = 0.;
				for (int k = 0; k < porousMedium_.r().size(); k++)
				{
					delta_rhobulk_due_to_single_reaction_[k](i) += porousMedium_.r_deposition_per_unit_volume_per_single_reaction()(k)*coefficient;
					delta_rhobulk_(i) += delta_rhobulk_due_to_single_reaction_[k](i);
				}

				// Single contributions (normalized)
				for (int k = 0; k < porousMedium_.r().size(); k++)
					delta_rhobulk_due_to_single_reaction_over_rhobulk_[k](i) = delta_rhobulk_due_to_single_reaction_[k](i) / (delta_rhobulk_(i)+1.e-12);
			}
		}

		if ( t>= (count_tecplot_+1)*tecplot_time_interval_)
		{
			std::stringstream current_index; current_index << count_tecplot_;
			std::string tecplot_file = "Solution.tec." + current_index.str();
			PrintTecplot(t, (output_folder_ / tecplot_file).string().c_str());

			count_tecplot_++;
		}

		t_old_ = t;
		count_video_++;
		count_file_++;
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
		for (int i = 0; i < NumberOfEquations(); i++)
			y_plus[i] = y[i];

		// Fixed values
		Equations(t, y, dy_original);

		// Derivatives with respect to y[kd]
		for (int kd = 0; kd < NumberOfEquations(); kd++)
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

	double Reactor2D::AreaAveraged(const Eigen::VectorXd& v)
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

			const double area_total = grid_x_.L()*grid_y_.L();
			return sum / area_total;
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

			const double Ri = grid_x_.x()(0);
			const double Re = grid_x_.x()(nx_-1);
			const double volume_total = boost::math::constants::pi<double>()*(Re*Re-Ri*Ri)*grid_y_.L();
			return sum / volume_total;
		}
	}

	double Reactor2D::AreaStandardDeviation(const double mean, const Eigen::VectorXd& v)
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
			const double area_total = grid_x_.L()*grid_y_.L();
			return std::sqrt(sum_std / area_total / coefficient);
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
}