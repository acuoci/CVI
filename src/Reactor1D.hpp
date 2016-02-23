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
							CVI::PorousMedium& porousMedium,
							OpenSMOKE::Grid1D& grid) :

	thermodynamicsMap_(thermodynamicsMap),
	kineticsMap_(kineticsMap),
	transportMap_(transportMap),
	porousMedium_(porousMedium),
	grid_(grid)
	{
		n_steps_video_ = 10;
		count_video_ = n_steps_video_;

		ns_ = thermodynamicsMap_.NumberOfSpecies();
		block_ = ns_ + 1;
		np_ = grid_.Np();
		ne_ = block_*np_;

		output_folder_ = "Output";

		MemoryAllocation();
		SetAlgebraicDifferentialEquations();
	}

	void Reactor1D::MemoryAllocation()
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
		omega_deposition_.resize(np_);

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
		depsilon_over_dt_.resize(np_);

		// Reset to zero
		for (int i = 0; i < np_-1; i++)
			j_star_[i].setZero();
		for (int i = 0; i < np_; i++)
			omega_homogeneous_[i].setZero();
		for (int i = 0; i < np_; i++)
			omega_heterogeneous_[i].setZero();
		omega_deposition_.setZero();
	}

	void Reactor1D::SetAlgebraicDifferentialEquations()
	{
		id_equations_.resize(ne_);

		unsigned int count = 0;

		// Interal boundary (symmetry plane)
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

		// Gas side
		{
			for (unsigned int j = 0; j < ns_; j++)
				id_equations_[count++] = false;
			id_equations_[count++] = true;
		}
	}

	void Reactor1D::SetGasSide(const double T_gas, const double P_gas, const OpenSMOKE::OpenSMOKEVectorDouble& omega_gas)
	{
		Y_gas_side_.resize(ns_);
		for (unsigned int j = 0; j < ns_; j++)
			Y_gas_side_(j) = omega_gas[j + 1];

		P_gas_side_ = P_gas;
		T_gas_side_ = T_gas;

		for (unsigned int j = 0; j < ns_; j++)
			Y_[grid_.Ni()](j) = omega_gas[j + 1];

		P_(grid_.Ni()) = P_gas;
		T_(grid_.Ni()) = T_gas;
	}

	void Reactor1D::SetInitialConditions(const double T_gas, const double P_gas, const OpenSMOKE::OpenSMOKEVectorDouble& omega_gas)
	{
		for (int i = 0; i < np_; i++)
			for (unsigned int j = 0; j < ns_; j++)
				Y_[i](j) = omega_gas[j + 1];

		P_.setConstant(P_gas);
		T_.setConstant(T_gas);
		epsilon_.setConstant(porousMedium_.porosity());
	}

	void Reactor1D::Properties()
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
					for (unsigned int j = 0; j < ns_; j++)
						aux_eigen(j) = aux_X[j + 1];
					porousMedium_.EffectiveDiffusionCoefficients(aux_eigen);
					for (unsigned int j = 0; j < ns_; j++)
						gamma_star_[i](j) = porousMedium_.gamma_effective()(j);
				}
			}

			// Kinetics
			{
				// Homogeneous phase
				{
					kineticsMap_.SetTemperature(T_(i));
					kineticsMap_.SetPressure(P_(i));
					kineticsMap_.ReactionRates(aux_C);
					kineticsMap_.FormationRates(&aux_R);
					ElementByElementProduct(aux_R, thermodynamicsMap_.MW(), &aux_R); // [kg/m3/s]
					aux_R.CopyTo(omega_homogeneous_[i].data());
				}

				// Heterogeneous phase
				{
					for (unsigned int j = 0; j < ns_; j++)
						aux_eigen(j) = aux_C[j + 1];
					porousMedium_.FormationRates(aux_eigen);
					for (unsigned int j = 0; j < ns_; j++)
						omega_heterogeneous_[i](j) = porousMedium_.Rgas()(j)*thermodynamicsMap_.MW()[j+1];			// [kg/m3/s]
					omega_deposition_(i) = porousMedium_.r_deposition_per_unit_volume()*porousMedium_.mw_carbon();	// [kg/m3/s]
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

	void Reactor1D::SpatialDerivatives()
	{
		Eigen::VectorXd dummy(np_);

		// First-order derivatives
		grid_.Derivative(OpenSMOKE::DERIVATIVE_1ST_BACKWARD, dummy, Y_, &dY_over_dx_);
		grid_.Derivative(OpenSMOKE::DERIVATIVE_1ST_FORWARD, dummy, gamma_star_, &dgamma_star_over_dx_);
		grid_.Derivative(OpenSMOKE::DERIVATIVE_1ST_FORWARD, dummy, rho_gas_, &drho_gas_over_dx_);

		// Second order derivatives
		grid_.SecondDerivative(Y_, &d2Y_over_dx2_);
	}

	void Reactor1D::SubEquations_MassFractions()
	{
		// Internal side (symmetry)
		for (unsigned int j = 0; j < ns_; j++)
			dY_over_dt_[0](j) = Y_[0](j) - Y_[1](j);

		for (int i = 1; i < grid_.Ni(); i++)
		{
			for (unsigned int j = 0; j < ns_; j++)
			{
				// Explicit derivatives
				dY_over_dt_[i](j) = gamma_star_[i](j)*rho_gas_(i)*d2Y_over_dx2_[i](j) + 
									(gamma_star_[i](j)*drho_gas_over_dx_(i) + dgamma_star_over_dx_[i](j)*rho_gas_(i))*dY_over_dx_[i](j) +
									epsilon_(i)*omega_homogeneous_[i](j) + omega_heterogeneous_[i](j) + Y_[i](j)*omega_deposition_(i);

				// Fluxes
				// dY_over_dt_[i](j) = -(j_star_[i](j) - j_star_[i - 1](j)) / grid_.dxc_over_2()(i) +
				//					epsilon_(i)*omega_homogeneous_[i](j) + omega_heterogeneous_[i](j) + Y_[i](j)*omega_deposition_(i);

				dY_over_dt_[i](j) /= (rho_gas_(i)*epsilon_(i));
			}

		}

		// Gas side
		for (unsigned int j = 0; j < ns_; j++)
			dY_over_dt_[grid_.Ni()](j) = Y_[grid_.Ni()](j) - Y_gas_side_(j);

		/*
		{
			const double dc = 2e-3;
			const double v = 1.;
			const double Re = rho_gas_(np_ - 1)*v*dc / 1.8e-5;
			const double Pr = 1.;
			const double Nu = 0.023*std::pow(Re, 0.8) * std::pow(Pr, 0.333);
			const double hc = Nu*1.e-5 / dc;

			std::cout << Re << " " << Nu << " " << hc << std::endl;
			for (unsigned int j = 0; j < ns_; j++)
			{
				const double Diff_over_dx = gamma_star_[np_ - 1](j) / (grid_.x()(np_ - 1) - grid_.x()(np_ - 2));
				dY_over_dt_[grid_.Ni()](j) = Y_[grid_.Ni()](j) - (Y_gas_side_(j)*hc + Y_[grid_.Ni() - 1](j)*Diff_over_dx) / (hc + Diff_over_dx);
			}
		}
		*/
	}

	void Reactor1D::SubEquations_Porosity()
	{
		// Internal points
		for (int i = 0; i < np_; i++)
			depsilon_over_dt_(i) = -omega_deposition_(i) / porousMedium_.rho_carbon();
	}

	void Reactor1D::Recover_Unknowns(const double* y)
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

	void Reactor1D::Recover_Residuals(double* dy)
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
		for (int i = 0; i < np_; i++)
		{
			for (unsigned int j = 0; j < ns_; j++)
				v[count++] = Y_[i](j);

			v[count++] = epsilon_[i];
		}
	}

	void Reactor1D::CorrectedUnknownsVector(double* v)
	{
		unsigned int count = 0;
		for (int i = 0; i < np_; i++)
		{
			for (unsigned int j = 0; j < ns_; j++)
				Y_[i](j) = v[count++];

			epsilon_(i) = v[count++];

			// Normalize
			//const double sum = Y_[i].sum();
			//for (unsigned int j = 0; j < ns_; j++)
			//	Y_[i](j) /= sum;
		}
	}

	void Reactor1D::MinimumUnknownsVector(double* v)
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

	void Reactor1D::MaximumUnknownsVector(double* v)
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

	void Reactor1D::Equations(const double t, const double* y, double* dy)
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
		SubEquations_Porosity();

		// Recover residuals
		Recover_Residuals(dy);
	}

	int Reactor1D::SolveFromScratch(DaeSMOKE::DaeSolver_Parameters& dae_parameters)
	{
		std::ofstream fMonitoring;
		fMonitoring.open((output_folder_ / "monitoring.out").string().c_str(), std::ios::out);
		fMonitoring.setf(std::ios::scientific);

		// Loop
		double time_interval = 3600.*10.;		// [s]
		unsigned int number_intervals = 12;
		for (unsigned int k = 1; k <= number_intervals; k++)
		{
			const double t0 = (k - 1)*time_interval;
			const double tf = t0 + time_interval;

			// Solve
			int flag = Solve(dae_parameters, t0, tf);
			if (flag < 0)
				return flag;

			// Write current solution
			std::stringstream number; number << k;
			std::string solution_file = "Solution." + number.str() + ".out";
			std::string diffusion_coefficients_file = "DiffusionCoefficients." + number.str() + ".out";
			std::string heterogeneous_rates_file = "HeterogeneousRates." + number.str() + ".out";
			std::string homogeneous_rates_file = "HomogeneousRates." + number.str() + ".out";
			PrintSolution(tf, (output_folder_ / solution_file).string().c_str());
			PrintDiffusionCoefficients(tf, (output_folder_ / diffusion_coefficients_file).string().c_str());
			PrintHeterogeneousRates(tf, (output_folder_ / heterogeneous_rates_file).string().c_str());
			PrintHomogeneousRates(tf, (output_folder_ / homogeneous_rates_file).string().c_str());
			//PrintXMLFile((output_folder_ / "Output.xml").string().c_str());
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
			fOutput << std::setprecision(9) << std::setw(20) << grid_.x()[i] * 1000.;
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

	void Reactor1D::PrintDiffusionCoefficients(const double t, const std::string name_file)
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
			fOutput << std::setprecision(9) << std::setw(20) << grid_.x()[i] * 1000.;
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

	void Reactor1D::PrintHomogeneousRates(const double t, const std::string name_file)
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

	void Reactor1D::PrintHeterogeneousRates(const double t, const std::string name_file)
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
			fOutput << std::setprecision(9) << std::setw(width) << grid_.x()[i] * 1000.;
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

	void Reactor1D::Print(const double t, const double* y)
	{
		if (count_video_ == n_steps_video_)
		{
			const double min_epsilon = epsilon_.minCoeff();
			const double max_epsilon = epsilon_.maxCoeff();

			std::cout << std::left << std::setw(16) << std::scientific << t;
			std::cout << std::left << std::setw(14) << std::fixed << std::setprecision(3) << min_epsilon;
			std::cout << std::left << std::setw(14) << std::fixed << std::setprecision(3) << max_epsilon;
			std::cout << std::endl;

			count_video_ = 0;
		}

		count_video_++;
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
}