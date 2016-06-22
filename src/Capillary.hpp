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
#include "Interface_Capillary_OpenSMOKEppDae.h"

namespace CVI
{
	Capillary::Capillary(	OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>& thermodynamicsMap,
				OpenSMOKE::KineticsMap_CHEMKIN<double>& kineticsMap,
				OpenSMOKE::TransportPropertiesMap_CHEMKIN<double>& transportMap,
				CVI::HeterogeneousMechanism& heterogeneousMechanism,
				OpenSMOKE::Grid1D& grid) :

	thermodynamicsMap_(thermodynamicsMap),
	kineticsMap_(kineticsMap),
	transportMap_(transportMap),
	heterogeneousMechanism_(heterogeneousMechanism),
	grid_(grid)
	{
		n_steps_video_ = 10;
		count_video_ = n_steps_video_;
		n_steps_file_ = 3;
		count_file_ = n_steps_file_;
		count_tecplot_ = 0;

		ns_ = thermodynamicsMap_.NumberOfSpecies();
		block_ = ns_ + 1;
		np_ = grid_.Np();
		ne_ = block_*np_;
		band_size_ = 2 * block_ - 1;

		diameter_initial_ = 0.;

		time_total_ = 48.*3600.;
		dae_time_interval_ = 3600.;
		tecplot_time_interval_ = 3600.;

		output_folder_ = "Output";

		output_matlab_folder_ = output_folder_ / "Matlab";
		OpenSMOKE::CreateDirectory(output_matlab_folder_);

		output_diffusion_folder_ = output_folder_ / "DiffusionCoefficients";
		OpenSMOKE::CreateDirectory(output_diffusion_folder_);

		output_heterogeneous_folder_ = output_folder_ / "HeterogeneousReactions";
		OpenSMOKE::CreateDirectory(output_heterogeneous_folder_);

		output_homogeneous_folder_ = output_folder_ / "HomogeneousReactions";
		OpenSMOKE::CreateDirectory(output_homogeneous_folder_);

		fMonitoring_.open((output_folder_ / "monitor.out").string().c_str(), std::ios::out);
		fMonitoring_.setf(std::ios::scientific);
		PrintLabelMonitoringFile();

		MemoryAllocation();
		SetAlgebraicDifferentialEquations();
	}

	void Capillary::MemoryAllocation()
	{
		OpenSMOKE::ChangeDimensions(ns_, &aux_X, true);
		OpenSMOKE::ChangeDimensions(ns_, &aux_Y, true);
		OpenSMOKE::ChangeDimensions(ns_, &aux_C, true);
		OpenSMOKE::ChangeDimensions(ns_, &aux_R, true);
		OpenSMOKE::ChangeDimensions(ns_, &aux_gamma, true);
		aux_eigen.resize(ns_);
		
		// Densities [kg/m3]
		rho_gas_.resize(np_);
		
		// Diameter [m]
		diameter_.resize(np_);

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

		// Diffusion fluxes [???]
		j_star_.resize(np_ - 1);
		for (int i = 0; i < np_ - 1; i++)
			j_star_[i].resize(ns_);

		// Correction diffusion flux [???]
		jc_star_.resize(np_ - 1);

		// Spatial derivatives: mole fractions [1/m]
		dX_over_dx_.resize(np_);
		for (int i = 0; i < np_; i++)
			dX_over_dx_[i].resize(ns_);

		// Spatial derivatives: mass fractions [1/m]
		dY_over_dx_.resize(np_);
		for (int i = 0; i < np_; i++)
			dY_over_dx_[i].resize(ns_);

		// Second order spatial derivatives: mass fractions [1/m2]
		d2Y_over_dx2_.resize(np_);
		for (int i = 0; i < np_; i++)
			d2Y_over_dx2_[i].resize(ns_);

		// Time derivatives: mass fractions [1/s]
		dY_over_dt_.resize(np_);
		for (int i = 0; i < np_; i++)
			dY_over_dt_[i].resize(ns_);

		// Spatial derivatives: mass diffusion coefficients [m2/s/m]
		dgamma_star_over_dx_.resize(np_);
		for (int i = 0; i < np_; i++)
			dgamma_star_over_dx_[i].resize(ns_);

		// Spatial derivatives: density of gaseous phase
		drho_gas_over_dx_.resize(np_);

		// Time derivatives: porosity [1/s]
		ddiameter_over_dt_.resize(np_);

		// Reset to zero
		for (int i = 0; i < np_-1; i++)
			j_star_[i].setZero();
		for (int i = 0; i < np_; i++)
			omega_homogeneous_[i].setZero();
		for (int i = 0; i < np_; i++)
			omega_heterogeneous_[i].setZero();
		omega_deposition_per_unit_volume_.setZero();
		omega_deposition_per_unit_area_.setZero();
	}


	void Capillary::SetAlgebraicDifferentialEquations()
	{
		id_equations_.resize(ne_);

		unsigned int count = 0;

		// Gas side (capillary mouth)
		{
			for (unsigned int j = 0; j < ns_; j++)
				id_equations_[count++] = false;
			id_equations_[count++] = true;
		}

		// Internal points
		for (int i = 1; i < grid_.Ni(); i++)
		{
			for (unsigned int j = 0; j < ns_; j++)
				id_equations_[count++] = true;
			id_equations_[count++] = true;
		}

		// Internal boundary (capillary depth)
		{
			for (unsigned int j = 0; j < ns_; j++)
				id_equations_[count++] = false;
			id_equations_[count++] = true;
		}
	}

	void Capillary::SetGasSide(const double T_gas, const double P_gas, const Eigen::VectorXd& omega_gas)
	{
		Y_gas_side_.resize(ns_);
		for (unsigned int j = 0; j < ns_; j++)
			Y_gas_side_(j) = omega_gas(j);

		P_gas_side_ = P_gas;
		T_gas_side_ = T_gas;

		// Capillary mouth (gas side)
		for (unsigned int j = 0; j < ns_; j++)
			Y_[0](j) = omega_gas(j);
		P_(0) = P_gas;
		T_(0) = T_gas;
	}

	void Capillary::SetInitialConditions(const double T_gas, const double P_gas, const double diameter, const Eigen::VectorXd& omega_gas)
	{
		for (int i = 0; i < np_; i++)
			for (unsigned int j = 0; j < ns_; j++)
				Y_[i](j) = omega_gas(j);

		P_.setConstant(P_gas);
		T_.setConstant(T_gas);
		diameter_.setConstant(diameter);
		diameter_initial_ = diameter;
	}

	void Capillary::SetTimeTotal(const double time_total)
	{
		time_total_ = time_total;
	}

	void Capillary::SetDaeTimeInterval(const double time_interval)
	{
		dae_time_interval_ = time_interval;
	}

	void Capillary::SetTecplotTimeInterval(const double time_interval)
	{
		tecplot_time_interval_ = time_interval;
	}

	void Capillary::SetStepsVideo(const int steps_video)
	{
		n_steps_video_ = steps_video;
		count_video_ = n_steps_video_;

	}

	void Capillary::SetStepsFile(const int steps_file)
	{
		n_steps_file_ = steps_file;
		count_file_ = n_steps_file_;
	}

	void Capillary::Properties()
	{
		for (int i = 0; i < np_; i++)
		{
			// Thermodynamics
			{
				thermodynamicsMap_.SetPressure(P_(i));
				thermodynamicsMap_.SetTemperature(T_(i));

				//aux_Y.CopyFrom(Y_[i].data());
				for (unsigned int j = 0; j < ns_; j++)
					aux_Y[j + 1] = Y_[i](j);
				thermodynamicsMap_.MoleFractions_From_MassFractions(aux_X, mw_(i), aux_Y);
				//aux_X.CopyTo(X_[i].data());
				for (unsigned int j = 0; j < ns_; j++)
					X_[i](j) = aux_X[j + 1];

				// Concentrations [kmol/m3]
				const double cTot = P_(i) / PhysicalConstants::R_J_kmol / T_(i); // [kmol/m3]
				Product(cTot, aux_X, &aux_C);

				// Mixture density
				rho_gas_(i) = cTot*mw_(i);	// [kg/m3]
			}

			// Transport properties (porous medium)
			{
				transportMap_.SetPressure(P_(i));
				transportMap_.SetTemperature(T_(i));

				// Mixture diffusion coefficients
				{
					transportMap_.MassDiffusionCoefficients(aux_gamma, aux_X, false);
					for (unsigned int j = 0; j < ns_; j++)
						gamma_star_[i](j) = aux_gamma[j + 1];
				}
			}

			// Kinetics
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
					aux_R.CopyTo(omega_homogeneous_[i].data());
				}

				// Heterogeneous phase
				if (heterogeneousMechanism_.heterogeneous_reactions() == true)
				{
					for (unsigned int j = 0; j < ns_; j++)
						aux_eigen(j) = aux_C[j + 1];
					heterogeneousMechanism_.FormationRates(4./diameter_(i), aux_eigen);
					for (unsigned int j = 0; j < ns_; j++)
						omega_heterogeneous_[i](j) = heterogeneousMechanism_.Rgas()(j)*thermodynamicsMap_.MW()[j+1];				// [kg/m3/s]
					omega_deposition_per_unit_area_(i) = heterogeneousMechanism_.r_deposition_per_unit_area()*heterogeneousMechanism_.mw_carbon();		// [kg/m2/s]
					omega_deposition_per_unit_volume_(i) = heterogeneousMechanism_.r_deposition_per_unit_volume()*heterogeneousMechanism_.mw_carbon();		// [kg/m3/s]
				}
			}
		}
	}

	void Capillary::DiffusionFluxes()
	{
		// Effective diffusion fluxes
		{
			Eigen::VectorXd dummy(np_); dummy.setZero();
			grid_.Derivative(OpenSMOKE::DERIVATIVE_1ST_BACKWARD, dummy, Y_, &dY_over_dx_);
			for (int i = 0; i < grid_.Ni(); i++)
			{
				j_star_[i].setZero();
				for (unsigned int j = 0; j < ns_; j++)
					j_star_[i](j) = -0.50 * (rho_gas_(i)*gamma_star_[i](j) + rho_gas_(i + 1)*gamma_star_[i + 1](j)) *dY_over_dx_[i + 1](j);
			}
		}

		// Correction diffusion velocity
		{
			jc_star_.setZero();
			for (int i = 0; i < grid_.Ni(); i++)
			{
				for (unsigned int j = 0; j < ns_; j++)
					jc_star_(i) -= j_star_[i](j);
			}
		}

		// Total diffusion velocity
		{
			for (int i = 0; i < grid_.Ni(); i++)
			{
				for (unsigned int j = 0; j < ns_; j++)
					j_star_[i](j) += 0.50*(Y_[i](j) + Y_[i + 1](j))*jc_star_(i);
			}
		}
	}

	void Capillary::SpatialDerivatives()
	{
		Eigen::VectorXd dummy(np_);

		// First-order derivatives
		grid_.Derivative(OpenSMOKE::DERIVATIVE_1ST_BACKWARD, dummy, Y_, &dY_over_dx_);
		grid_.Derivative(OpenSMOKE::DERIVATIVE_1ST_FORWARD, dummy, gamma_star_, &dgamma_star_over_dx_);
		grid_.Derivative(OpenSMOKE::DERIVATIVE_1ST_FORWARD, dummy, rho_gas_, &drho_gas_over_dx_);

		// Second order derivatives
		grid_.SecondDerivative(Y_, &d2Y_over_dx2_);
	}

	void Capillary::SubEquations_MassFractions()
	{
		// Gas side
		for (unsigned int j = 0; j < ns_; j++)
			dY_over_dt_[0](j) = Y_[0](j) - Y_gas_side_(j);

		for (int i = 1; i < grid_.Ni(); i++)
		{
			for (unsigned int j = 0; j < ns_; j++)
			{
				dY_over_dt_[i](j) = 	gamma_star_[i](j)*rho_gas_(i)*d2Y_over_dx2_[i](j) +
							(gamma_star_[i](j)*drho_gas_over_dx_(i) + dgamma_star_over_dx_[i](j)*rho_gas_(i))*dY_over_dx_[i](j) +
							omega_homogeneous_[i](j) + omega_heterogeneous_[i](j) + Y_[i](j)*omega_deposition_per_unit_volume_(i);

				dY_over_dt_[i](j) /= rho_gas_(i);
			}

		}

		// Internal side (symmetry)
		for (unsigned int j = 0; j < ns_; j++)
			dY_over_dt_[grid_.Ni()](j) = Y_[grid_.Ni()](j) - Y_[grid_.Ni()-1](j);
	}

	void Capillary::SubEquations_Diameter()
	{
		// Internal points
		for (int i = 0; i < np_; i++)
			ddiameter_over_dt_(i) = -2.*omega_deposition_per_unit_area_(i) / heterogeneousMechanism_.rho_graphite();
	}

	void Capillary::Recover_Unknowns(const double* y)
	{
		unsigned int count = 0;
		for (int i = 0; i < np_; i++)
		{
			// Species
			for (unsigned int j = 0; j < ns_; j++)
				Y_[i](j) = y[count++];

			// Porosity
			diameter_(i) = y[count++];
		}
	}

	void Capillary::Recover_Residuals(double* dy)
	{
		unsigned int count = 0;
		for (int i = 0; i < np_; i++)
		{
			// Species
			for (unsigned int j = 0; j < ns_; j++)
				dy[count++] = dY_over_dt_[i](j);
			
			// Porosity
			dy[count++] = ddiameter_over_dt_(i);
		}
	}

	void Capillary::AlgebraicDifferentialVector(double* v)
	{
		int count = 0;
		for (unsigned int i = 0; i < id_equations_.size(); i++)
		{
			if (id_equations_[i] == true)  v[count++] = 1.;
			if (id_equations_[i] == false) v[count++] = 0.;
		}
	}

	void Capillary::UnknownsVector(double* v)
	{
		unsigned int count = 0;
		for (int i = 0; i < np_; i++)
		{
			for (unsigned int j = 0; j < ns_; j++)
				v[count++] = Y_[i](j);

			v[count++] = diameter_[i];
		}
	}

	void Capillary::CorrectedUnknownsVector(double* v)
	{
		unsigned int count = 0;
		for (int i = 0; i < np_; i++)
		{
			for (unsigned int j = 0; j < ns_; j++)
				Y_[i](j) = v[count++];

			diameter_(i) = v[count++];

			// Normalize
			//const double sum = Y_[i].sum();
			//for (unsigned int j = 0; j < ns_; j++)
			//	Y_[i](j) /= sum;
		}
	}

	void Capillary::MinimumUnknownsVector(double* v)
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

	void Capillary::MaximumUnknownsVector(double* v)
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

	void Capillary::Equations(const double t, const double* y, double* dy)
	{
		// Recover unknowns
		Recover_Unknowns(y);

		// Properties
		Properties();

		// Spatial derivatives
		SpatialDerivatives();

		// Fluxes
		DiffusionFluxes();

		// Equations
		SubEquations_MassFractions();
		SubEquations_Diameter();

		// Recover residuals
		Recover_Residuals(dy);
	}

	int Capillary::SolveFromScratch(DaeSMOKE::DaeSolver_Parameters& dae_parameters)
	{
		std::ofstream fMonitoring;
		fMonitoring.open((output_folder_ / "log.out").string().c_str(), std::ios::out);
		fMonitoring.setf(std::ios::scientific);

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
			PrintDiffusionCoefficients(tf, (output_diffusion_folder_ / diffusion_coefficients_file).string().c_str());
			PrintHeterogeneousRates(tf, (output_heterogeneous_folder_ / heterogeneous_rates_file).string().c_str());
			PrintHomogeneousRates(tf, (output_homogeneous_folder_ / homogeneous_rates_file).string().c_str());
		}

		return true;
	}

	int Capillary::Solve(DaeSMOKE::DaeSolver_Parameters& dae_parameters, const double t0, const double tf)
	{
		int flag = DaeSMOKE::Solve_Band_OpenSMOKEppDae<Capillary, OpenSMOKE_Capillary_DaeSystem>(this, dae_parameters, t0, tf);
		return flag;
	}

	void Capillary::PrintSolution(const double t, const std::string name_file)
	{
		std::ofstream fOutput(name_file.c_str(), std::ios::out);

		{
			unsigned int count = 1;
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "time[s]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "x[mm]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "dummy", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "T[K]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "P[Pa]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "D[mm]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "s[micron]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "Sv[1/m]", count);
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
			fOutput << std::setprecision(9) << std::setw(20) << grid_.x()[i] * 1000.;
			fOutput << std::setprecision(9) << std::setw(20) << 0.;
			fOutput << std::setprecision(9) << std::setw(20) << T_(i);
			fOutput << std::setprecision(9) << std::setw(20) << P_(i);
			
			fOutput << std::setprecision(9) << std::setw(20) << diameter_(i)*1000.;
			fOutput << std::setprecision(9) << std::setw(20) << 0.50*(diameter_initial_-diameter_(i))*1.e6;
			fOutput << std::setprecision(9) << std::setw(20) << 4./diameter_(i);

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

	void Capillary::PrintDiffusionCoefficients(const double t, const std::string name_file)
	{
		std::ofstream fOutput(name_file.c_str(), std::ios::out);

		{
			unsigned int count = 1;
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "time[s]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "x[mm]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "T[K]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "P[Pa]", count);

			for (unsigned int j = 0; j < ns_; j++)
				OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, thermodynamicsMap_.NamesOfSpecies()[j], count);

			fOutput << std::endl;
		}

		for (int i = 0; i < np_; i++)
		{
			// Calculate the diffusion coefficients
			{
				thermodynamicsMap_.SetPressure(P_(i));
				thermodynamicsMap_.SetTemperature(T_(i));
				transportMap_.SetPressure(P_(i));
				transportMap_.SetTemperature(T_(i));

				aux_Y.CopyFrom(Y_[i].data());
				thermodynamicsMap_.MoleFractions_From_MassFractions(aux_X, mw_(i), aux_Y);
				aux_X.CopyTo(X_[i].data());
				transportMap_.MassDiffusionCoefficients(aux_gamma, aux_X, false);
			}

			fOutput << std::setprecision(9) << std::setw(20) << t;
			fOutput << std::setprecision(9) << std::setw(20) << grid_.x()[i] * 1000.;
			fOutput << std::setprecision(9) << std::setw(20) << T_(i);
			fOutput << std::setprecision(9) << std::setw(20) << P_(i);

			// Diffusion coefficients
			for (unsigned int j = 0; j < ns_; j++)
			{
				fOutput << std::setprecision(9) << std::setw(20) << aux_gamma[j+1];
			}

			fOutput << std::endl;
		}

		fOutput.close();
	}

	void Capillary::PrintHomogeneousRates(const double t, const std::string name_file)
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
			fOutput << std::setprecision(9) << std::setw(20) << grid_.x()[i] * 1000.;
			fOutput << std::setprecision(9) << std::setw(20) << T_(i);
			fOutput << std::setprecision(9) << std::setw(20) << P_(i);

			for (unsigned int j = 0; j < ns_; j++)
					fOutput << std::setprecision(9) << std::setw(20) << aux_R[j+1];

			fOutput << std::endl;
		}

		fOutput.close();
	}

	void Capillary::PrintHeterogeneousRates(const double t, const std::string name_file)
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

		for (int i = 0; i < np_; i++)
		{
			// Calculate the reaction and formation rates of heterogeneous reactions
			{
				heterogeneousMechanism_.SetTemperature(T_(i));
				heterogeneousMechanism_.SetPressure(P_(i));

				thermodynamicsMap_.SetPressure(P_(i));
				thermodynamicsMap_.SetTemperature(T_(i));
				aux_Y.CopyFrom(Y_[i].data());
				thermodynamicsMap_.MoleFractions_From_MassFractions(aux_X, mw_(i), aux_Y);
				aux_X.CopyTo(X_[i].data());
				const double cTot = P_(i) / PhysicalConstants::R_J_kmol / T_(i); // [kmol/m3]
				Product(cTot, aux_X, &aux_C);

				for (unsigned int j = 0; j < ns_; j++)
					aux_eigen(j) = aux_C[j + 1];
				heterogeneousMechanism_.FormationRates(4./diameter_(i), aux_eigen);
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

	void Capillary::Print(const double t, const double* y)
	{
		if (count_video_ == n_steps_video_)
		{
			const double min_diameter = diameter_.minCoeff();
			const double max_diameter = diameter_.maxCoeff();

			std::cout << std::left << std::setw(16) << std::scientific << t;
			std::cout << std::left << std::setw(14) << std::fixed << std::setprecision(3) << min_diameter;
			std::cout << std::left << std::setw(14) << std::fixed << std::setprecision(3) << max_diameter;
			std::cout << std::endl;

			count_video_ = 0;
		}

		if (count_file_ == n_steps_file_)
		{
			Eigen::VectorXd thickness(np_);
			Eigen::VectorXd Sv(np_);
			for (int i = 0; i < np_; i++)
			{
				thickness(i) = 0.50*(diameter_initial_-diameter_(i));
				Sv(i) = 4./diameter_(i);
			}

			const int width = 20;

			const double T_mean = AreaAveraged(T_);
			const double T_std = AreaStandardDeviation(T_mean, T_);

			const double P_mean = AreaAveraged(P_);
			const double P_std = AreaStandardDeviation(P_mean, P_);

			const double diameter_mean = AreaAveraged(diameter_);
			const double diameter_std = AreaStandardDeviation(diameter_mean, diameter_);

			const double thickness_mean = AreaAveraged(thickness);
			const double thickness_std = AreaStandardDeviation(thickness_mean, thickness);

			const double Sv_mean = AreaAveraged(Sv);
			const double Sv_std = AreaStandardDeviation(Sv_mean, Sv);

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

			fMonitoring_ << std::left << std::setw(width) << std::fixed << std::setprecision(4) << diameter_mean*1e3;	// [mm]
			fMonitoring_ << std::left << std::setw(width) << std::fixed << std::setprecision(4) << diameter_std*1e3;	// [mm]

			fMonitoring_ << std::left << std::setw(width) << std::fixed << std::setprecision(4) << thickness_mean*1e6;	// [micron]
			fMonitoring_ << std::left << std::setw(width) << std::fixed << std::setprecision(4) << thickness_std*1e6;	// [micron]

			fMonitoring_ << std::left << std::setw(width) << std::fixed << std::setprecision(4) << Sv_mean;	// [1/m]
			fMonitoring_ << std::left << std::setw(width) << std::fixed << std::setprecision(4) << Sv_std;	// [1/m]

			fMonitoring_ << std::left << std::setw(width) << std::scientific << std::setprecision(6) << r_deposition_per_unit_area_mean / heterogeneousMechanism_.rho_graphite()*1000.;		// [mm/s]
			fMonitoring_ << std::left << std::setw(width) << std::scientific << std::setprecision(6) << r_deposition_per_unit_area_std / heterogeneousMechanism_.rho_graphite()*1000.;	// [mm/s]

			fMonitoring_ << std::left << std::setw(width) << std::scientific << std::setprecision(6) << r_deposition_per_unit_area_mean *1000.*3600.;	// [g/m2/h]
			fMonitoring_ << std::left << std::setw(width) << std::scientific << std::setprecision(6) << r_deposition_per_unit_area_std *1000.*3600.;		// [g/m2/h]

			fMonitoring_ << std::left << std::setw(width) << std::scientific << std::setprecision(6) << r_deposition_per_unit_volume_mean *1000.*3600.;	// [g/m3/h]
			fMonitoring_ << std::left << std::setw(width) << std::scientific << std::setprecision(6) << r_deposition_per_unit_volume_std *1000.*3600.;	// [g/m3/h]

			fMonitoring_ << std::endl;

			count_file_ = 0;
		}

		count_video_++;
		count_file_++;
	}

	void Capillary::SparsityPattern(std::vector<unsigned int>& rows, std::vector<unsigned int>& cols)
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

	double Capillary::AreaAveraged(const Eigen::VectorXd& v)
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

	double Capillary::AreaStandardDeviation(const double mean, const Eigen::VectorXd& v)
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

	void Capillary::PrintLabelMonitoringFile()
	{
		unsigned int count = 1;
		const unsigned int width = 20;
		OpenSMOKE::PrintTagOnASCIILabel(width, fMonitoring_, "time[h]", count);
		OpenSMOKE::PrintTagOnASCIILabel(width, fMonitoring_, "time[s]", count);
		OpenSMOKE::PrintTagOnASCIILabel(width, fMonitoring_, "T[K]", count);
		OpenSMOKE::PrintTagOnASCIILabel(width, fMonitoring_, "Tstd[K]", count);
		OpenSMOKE::PrintTagOnASCIILabel(width, fMonitoring_, "Pa[Pa]", count);
		OpenSMOKE::PrintTagOnASCIILabel(width, fMonitoring_, "Pastd[Pa]", count);
		OpenSMOKE::PrintTagOnASCIILabel(width, fMonitoring_, "D[mm]", count);
		OpenSMOKE::PrintTagOnASCIILabel(width, fMonitoring_, "Dstd[mm]", count);
		OpenSMOKE::PrintTagOnASCIILabel(width, fMonitoring_, "s[micron]", count);
		OpenSMOKE::PrintTagOnASCIILabel(width, fMonitoring_, "sstd[micron]", count);
		OpenSMOKE::PrintTagOnASCIILabel(width, fMonitoring_, "Sv[1/m]", count);
		OpenSMOKE::PrintTagOnASCIILabel(width, fMonitoring_, "Svstd[1/m]", count);
		OpenSMOKE::PrintTagOnASCIILabel(width, fMonitoring_, "rDep[mm/s]", count);
		OpenSMOKE::PrintTagOnASCIILabel(width, fMonitoring_, "rDepstd[mm/s]", count);
		OpenSMOKE::PrintTagOnASCIILabel(width, fMonitoring_, "rDep[g/m2/h]", count);
		OpenSMOKE::PrintTagOnASCIILabel(width, fMonitoring_, "rDepstd[g/m2/h]", count);
		OpenSMOKE::PrintTagOnASCIILabel(width, fMonitoring_, "rDep[g/m3/h]", count);
		OpenSMOKE::PrintTagOnASCIILabel(width, fMonitoring_, "rDepstd[g/m3/h]", count);

		fMonitoring_ << std::endl;
	}
}
