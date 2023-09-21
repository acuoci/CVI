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
	PorousMedium::PorousMedium(	OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap,
								OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap,
								OpenSMOKE::TransportPropertiesMap_CHEMKIN& transportMap,
								PorousSubstrateType porous_substrate_type, const double rf, const double rho_fiber, const double epsilon0) :

	thermodynamicsMap_(thermodynamicsMap),
	kineticsMap_(kineticsMap),
	transportMap_(transportMap),
	rf_(rf),
	rho_fiber_(rho_fiber),
	epsilon0_(epsilon0),
	porous_substrate_type_(porous_substrate_type)
	{
		Initialize();
	}

	PorousMedium::PorousMedium(	OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap,
								OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap,
								OpenSMOKE::TransportPropertiesMap_CHEMKIN& transportMap,
								OpenSMOKE::OpenSMOKE_Dictionary& dictionary ) :
	thermodynamicsMap_(thermodynamicsMap),
	kineticsMap_(kineticsMap),
	transportMap_(transportMap)
	{
		// Read fiber radius
		{
			double value;
			std::string units;
			if (dictionary.CheckOption("@FiberRadius") == true)
			{
				dictionary.ReadMeasure("@FiberRadius", value, units);
				if (units == "m")				rf_ = value;
				else if (units == "cm")			rf_ = value / 1.e2;
				else if (units == "mm")			rf_ = value / 1.e3;
				else if (units == "micron")		rf_ = value / 1.e6;
				else OpenSMOKE::FatalErrorMessage("@FiberRadius: Unknown fiber radius units. Available units: m | cm | mm | micron");
			}
		}

		// Read fiber density
		{
			double value;
			std::string units;
			if (dictionary.CheckOption("@FiberDensity") == true)
			{
				dictionary.ReadMeasure("@FiberDensity", value, units);
				if (units == "kg/m3")			rho_fiber_ = value;
				else if (units == "g/cm3")		rho_fiber_ = value*1.e3;
				else OpenSMOKE::FatalErrorMessage("@FiberDensity: Unknown fiber density units. Available units: kg/m3 | g/cm3");
			}
		}

		// Read initial porosity
		if (dictionary.CheckOption("@InitialPorosity") == true)
			dictionary.ReadDouble("@InitialPorosity", epsilon0_);

		// Read initial porosity
		mass_diffusion_multiplier_ = 1.;
		if (dictionary.CheckOption("@MassDiffusionMultiplier") == true)
			dictionary.ReadDouble("@MassDiffusionMultiplier", mass_diffusion_multiplier_);

		// Threshold porosity
		epsilon_threshold_ = 1e-2;
			if (dictionary.CheckOption("@ThresholdPorosity") == true)
				dictionary.ReadDouble("@ThresholdPorosity", epsilon_threshold_);

		// Thrrshold porosity
		epsilon_smoothing_coefficient_ = 200.;
		if (dictionary.CheckOption("@SmoothingCoefficientPorosity") == true)
			dictionary.ReadDouble("@SmoothingCoefficientPorosity", epsilon_smoothing_coefficient_);

		// Read porous substrate type
		{
			std::string value;
			if (dictionary.CheckOption("@PorousSubstrate") == true)
			{
				dictionary.ReadString("@PorousSubstrate", value);
				if (value == "polynomial")						porous_substrate_type_ = CVI::POLYNOMIAL;
				else if (value == "random")						porous_substrate_type_ = CVI::RANDOM;
				else if (value == "random_hardcore")			porous_substrate_type_ = CVI::RANDOM_HARDCORE;
				else if (value == "polynomial_onehalf")			porous_substrate_type_ = CVI::POLINOMIAL_ONEHALF;
				else if (value == "from_spheres_to_cylinders")	porous_substrate_type_ = CVI::FROM_SPHERES_TO_CYLINDERS;
				else if (value == "deutschmann_correlation")	porous_substrate_type_ = CVI::DEUTSCHMANN_CORRELATION;
				else if (value == "tang_felt_correlation")	porous_substrate_type_ = CVI::TANG_FELT_CORRELATION;
				else if (value == "tang_cloth_correlation")	porous_substrate_type_ = CVI::TANG_CLOTH_CORRELATION;
				else OpenSMOKE::FatalErrorMessage("@PorousSubstrate: Substrates available: polynomial | random | random_hardcore | polynomial_onehalf | from_spheres_to_cylinders | deutschmann_correlation | tang_felt_correlation | tang_cloth_correlation");
			}
		}

		// Porous substrate correction coefficient
		porous_substrate_correction_coefficient_ = 1.;
		if (dictionary.CheckOption("@PorousSubstrateCorrectionCoefficient") == true)
			dictionary.ReadDouble("@PorousSubstrateCorrectionCoefficient", porous_substrate_correction_coefficient_);

		Initialize();
	}

	void PorousMedium::Initialize()
	{
		// Number of species
		ns_ = thermodynamicsMap_.NumberOfSpecies();

		// Memory allocation
		gamma_knudsen_.resize(ns_);
		gamma_knudsen_effective_.resize(ns_);
		gamma_fick_.resize(ns_);
		gamma_fick_effective_.resize(ns_);
		gamma_effective_.resize(ns_);

		// Set initial porosity
		SetPorosity(epsilon0_);

		// Summary on the screen
		Summary();
	}

	void PorousMedium::Summary()
	{
		std::cout << std::endl;
		std::cout << "-------------------------------------------------------" << std::endl;
		std::cout << "                      Porous Medium                    " << std::endl;
		std::cout << "-------------------------------------------------------" << std::endl;
		std::cout << " * Porosity [-]:                     " << epsilon0_ << std::endl;
		std::cout << " * Fiber density [kg/m3]:            " << rho_fiber_ << std::endl;
		std::cout << " * Fiber radius [micron]:            " << rf_*1e6 << std::endl;
		std::cout << " * Pore radius [micron]:             " << rp()*1e6 << std::endl;
		std::cout << " * Surface per unit of volume [1/m]: " << Sv() << std::endl;
		std::cout << " * Permeability [m2]:                " << permeability() << std::endl;
		std::cout << " * Tortuosity (ordinary) [-]:        " << eta_bulk() << std::endl;
		std::cout << " * Tortuosity (Knudsen) [-]:         " << eta_knudsen() << std::endl;
		std::cout << " * Tortuosity (viscous) [-]:         " << eta_viscous() << std::endl;
		std::cout << "--------------------------------------------------------" << std::endl;
		std::cout << std::endl;
	}

	void PorousMedium::SetTemperature(const double T)
	{
		T_ = T;
	}

	void PorousMedium::SetPressure(const double P_Pa)
	{
		P_Pa_ = P_Pa;
	}
	
	void PorousMedium::SetPorosity(const double epsilon)
	{
		epsilon_ = epsilon;
	}

	double PorousMedium::density_bulk(const double rho_graphite)
	{
		return rho_fiber_*(1. - epsilon0_) + rho_graphite*(epsilon0_ - epsilon_);
	}

	double PorousMedium::Sv()
	{
		const double small = 1.e-16;
		const double epsilon = epsilon_ + small;

		double sv = 0.;
		if (porous_substrate_type_ == POLYNOMIAL)
		{
			const double ratio = epsilon / epsilon0_;
			sv = 2. / rf_*((2. - epsilon0_)*ratio - std::pow(ratio, 2.));
		}
		else if (porous_substrate_type_ == RANDOM)
		{
			sv = -2. / rf_*epsilon*std::log(epsilon);
		}
		else if (porous_substrate_type_ == RANDOM_HARDCORE)
		{
			const double ratio = epsilon / epsilon0_;
			sv = 2. / rf_*ratio*(1. - epsilon)*(1. / epsilon0_*std::log(1. / ratio) + 1.);
		}
		else if (porous_substrate_type_ == POLINOMIAL_ONEHALF)
		{
			const double ratio = epsilon / epsilon0_;
			sv = 2. / rf_*((2.-3./2.*epsilon)*ratio-(1.-epsilon0_/2.)*ratio*ratio);
		}
		else if (porous_substrate_type_ == FROM_SPHERES_TO_CYLINDERS)
		{
			const double c = 3. / 4. / PhysicalConstants::pi;
			sv = 2. / rf_*(1. - epsilon)*(std::pow(-c*std::log(1 - epsilon), 2. / 3.) + epsilon) / (std::pow(-c*std::log(1 - epsilon0_), 2. / 3.) + epsilon0_);
		}
		else if (porous_substrate_type_ == DEUTSCHMANN_CORRELATION)
		{
			const double epsilon2 = epsilon*epsilon;
			const double epsilon3 = epsilon2*epsilon;
			const double epsilon4 = epsilon2*epsilon2;
			const double epsilon5 = epsilon3*epsilon2;
			sv = -3.855762523198E+05*epsilon5 + 8.558541857891E+05*epsilon4 - 6.109196594973E+05*epsilon3 - 4.351758023548E+04*epsilon2 + 2.196529832093E+05*epsilon;
		}
		else if (porous_substrate_type_ == TANG_CLOTH_CORRELATION)
		{
			const double epsilon2 = epsilon*epsilon;
			const double epsilon3 = epsilon2*epsilon;
			const double epsilon4 = epsilon2*epsilon2;
			
			sv = 1647969.227*epsilon4 -802950.891*epsilon3 -29868.7377*epsilon2 +55980.45249*epsilon;
		}
		else if (porous_substrate_type_ == TANG_FELT_CORRELATION)
		{
			const double epsilon2 = epsilon*epsilon;
			const double epsilon3 = epsilon2*epsilon;
			const double epsilon4 = epsilon2*epsilon2;
			
			sv = -1.335107101E+05*epsilon2 +1.479173135E+05*epsilon;
			
		}

		return sv*porous_substrate_correction_coefficient_;
	}

	double PorousMedium::rp()
	{
		if (porous_substrate_type_ == DEUTSCHMANN_CORRELATION)
		{
			const double epsilon = epsilon_;
			const double epsilon2 = epsilon*epsilon;
			const double epsilon3 = epsilon2*epsilon;
			const double epsilon4 = epsilon2*epsilon2;
			const double epsilon5 = epsilon3*epsilon2;

			return 1.765902994055E-06*epsilon5 - 6.507266939515E-06*epsilon4 + 1.394874313698E-05*epsilon3 - 2.219345484444E-05*epsilon2 + 2.741053445908E-05*epsilon;
		}
		else
		{
			return 2. / Sv()*epsilon_;
		}
		
	}

	double PorousMedium::eta_bulk()
	{
		return 1./std::pow(epsilon_, 2./3.);
	}

	double PorousMedium::eta_knudsen()
	{
		return 1.444 / epsilon_;
	}

	double PorousMedium::eta_viscous()
	{
		return 2.76 / std::pow(epsilon_, 2. / 3.) * std::pow(std::log(epsilon_), 2.);
	}

	double PorousMedium::permeability()
	{
		return epsilon_*std::pow(rp(), 2.) / 8. / eta_viscous();
	}

	void PorousMedium::KnudsenDiffusionCoefficients()
	{
		const double coefficient = 2. / 3.*rp()*std::sqrt(8.*PhysicalConstants::R_J_kmol*T_ / PhysicalConstants::pi);
		
		for (unsigned int i = 0; i < ns_; i++)
			gamma_knudsen_(i) = coefficient/std::sqrt(thermodynamicsMap_.MW(i));
	}

	void PorousMedium::KnudsenEffectiveDiffusionCoefficients()
	{
		KnudsenDiffusionCoefficients();

		const double correction = epsilon_ / eta_knudsen();
		gamma_knudsen_effective_ = gamma_knudsen_*correction;
	}

	void PorousMedium::FickDiffusionCoefficients(const Eigen::VectorXd& mole_fractions)
	{
		// Set temperature and pressure of transport map
		transportMap_.SetTemperature(T_);
		transportMap_.SetPressure(P_Pa_);

		// Fick diffusion coefficients
		{
			OpenSMOKE::OpenSMOKEVectorDouble gamma_fick_os(ns_);
			OpenSMOKE::OpenSMOKEVectorDouble mole_fractions_os(ns_);

			for (unsigned int i = 0; i < ns_; i++)
				mole_fractions_os[i + 1] = mole_fractions(i);
			transportMap_.MassDiffusionCoefficients(gamma_fick_os.GetHandle(), mole_fractions_os.GetHandle(), false);
			for (unsigned int i = 0; i < ns_; i++)
				gamma_fick_(i) = gamma_fick_os[i + 1];
		}
	}

	void PorousMedium::FickEffectiveDiffusionCoefficients(const Eigen::VectorXd& mole_fractions)
	{
		FickDiffusionCoefficients(mole_fractions);

		const double correction = epsilon_ / eta_bulk();
		gamma_fick_effective_ = gamma_fick_*correction;
	}

	void PorousMedium::EffectiveDiffusionCoefficients(const Eigen::VectorXd& mole_fractions)
	{
		// Fick effective diffusion coefficients
		FickEffectiveDiffusionCoefficients(mole_fractions);

		// Knudsen effective diffusion coefficients
		KnudsenEffectiveDiffusionCoefficients();

		// Combination of diffusion coefficients
		for (unsigned int i = 0; i < ns_; i++)
			gamma_effective_(i) = 1. / (1. / gamma_fick_effective_(i) + 1. / gamma_knudsen_effective_(i));

		for (unsigned int i = 0; i < ns_; i++)
			gamma_effective_(i) *= mass_diffusion_multiplier_; 
	}

	double PorousMedium::lambda_gas()
	{
		return 1.5207e-11*T_*T_*T_-4.8574e-08*T_*T_+1.0184e-04*T_-3.9333e-04;	// [W/m/K]
	}

	double PorousMedium::lambda_fiber()
	{
		return 1.1176 - 6.217e-4*T_ + 8.3e-7*T_*T_;	// [W/m/K]
	}

	double PorousMedium::lambda_graphite()
	{
		return -3.466 + 0.0271*T_ - 2.05e-5*T_*T_ + 5.3e-9*T_*T_*T_;	// [W/m/K]
	}

	double PorousMedium::lambda()
	{
		return	epsilon_*lambda_gas()			+ 
				(1.-epsilon0_)*lambda_fiber()	+ 
				(epsilon0_-epsilon_)*lambda_graphite();	// [W/m/K]
	}

	double PorousMedium::cp_gas()
	{
		return 1.9327e-10*T_*T_*T_*T_-7.9999e-07*T_*T_*T_
				+1.1407e-03*T_*T_-4.4890e-01*T_ + 1.0575e+03;	// [J/kg/K]
	}

	double PorousMedium::cp_fiber()
	{
		return cp_graphite();	// [W/m/K]
	}

	double PorousMedium::cp_graphite()
	{
		return -42.4678 + 2.8516*T_ - 0.001*T_*T_;	// [W/m/K]
	}

	double PorousMedium::cp_times_rho(const double rho_gas, const double rho_graphite)
	{
		return	epsilon_*rho_gas*cp_gas() +
				(1. - epsilon0_)*rho_fiber_*cp_fiber() +
				(epsilon0_ - epsilon_)*rho_graphite*cp_graphite();	// [J/m3/K]
	}
}
