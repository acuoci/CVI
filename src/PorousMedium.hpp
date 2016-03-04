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
	PorousMedium::PorousMedium(	OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>& thermodynamicsMap,
								OpenSMOKE::KineticsMap_CHEMKIN<double>& kineticsMap,
								OpenSMOKE::TransportPropertiesMap_CHEMKIN<double>& transportMap,
								PorousSubstrateType porous_substrate_type, const double rf, const double rho_fiber, const double epsilon0,
								const HeterogeneousMechanism heterogeneous_mechanism_type,
								const HydrogenInhibitionType hydrogen_inhibition_type) :

	thermodynamicsMap_(thermodynamicsMap),
	kineticsMap_(kineticsMap),
	transportMap_(transportMap),
	rf_(rf),
	rho_fiber_(rho_fiber),
	epsilon0_(epsilon0),
	porous_substrate_type_(porous_substrate_type),
	heterogeneous_mechanism_type_(heterogeneous_mechanism_type),
	hydrogen_inhibition_type_(hydrogen_inhibition_type)
	{
		Initialize();
	}

	PorousMedium::PorousMedium(	OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>& thermodynamicsMap,
								OpenSMOKE::KineticsMap_CHEMKIN<double>& kineticsMap,
								OpenSMOKE::TransportPropertiesMap_CHEMKIN<double>& transportMap,
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
				else OpenSMOKE::FatalErrorMessage("@FiberRadius: Unknown fiber radius units");
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
				else OpenSMOKE::FatalErrorMessage("@FiberDensity: Unknown fiber density units");
			}
		}

		// Read graphite density
		{
			double value;
			std::string units;
			if (dictionary.CheckOption("@GraphiteDensity") == true)
			{
				dictionary.ReadMeasure("@GraphiteDensity", value, units);
				if (units == "kg/m3")			rho_graphite_ = value;
				else if (units == "g/cm3")		rho_graphite_ = value*1.e3;
				else OpenSMOKE::FatalErrorMessage("@GraphiteDensity: Unknown graphite density units");
			}
		}

		// Read initial porosity
		if (dictionary.CheckOption("@InitialPorosity") == true)
			dictionary.ReadDouble("@InitialPorosity", epsilon0_);

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
				else OpenSMOKE::FatalErrorMessage("@PorousSubstrate: Substrates available: polynomial | random | random_hardcore | polynomial_onehalf | from_spheres_to_cylinders");
			}
		}

		// Read heterogeneous mechanism type
		{
			std::string value;
			if (dictionary.CheckOption("@HeterogeneousMechanism") == true)
			{
				dictionary.ReadString("@HeterogeneousMechanism", value);
				if (value == "Ibrahim-Paolucci")		heterogeneous_mechanism_type_ = CVI::IBRAHIM_PAOLUCCI;
				else OpenSMOKE::FatalErrorMessage("Heterogeneous mechanisms available: Ibrahim-Paolucci");
			}
		}

		// Read hydrogen inhibition type
		{
			std::string value;
			if (dictionary.CheckOption("@HydrogenInhibition") == true)
			{
				dictionary.ReadString("@HydrogenInhibition", value);
				if (value == "none")			hydrogen_inhibition_type_ = CVI::NONE;
				else if (value == "Becker")		hydrogen_inhibition_type_ = CVI::BECKER;
				else OpenSMOKE::FatalErrorMessage("Hydrogen inhibitions available: none | Becker");
			}
		}

		Initialize();
	}

	void PorousMedium::Initialize()
	{
		// Constants
		mw_carbon_ = 12.010999679565430;	// [kg/kmol]

		// Number of species
		ns_ = thermodynamicsMap_.NumberOfSpecies();

		// Number of heterogeneous reactions
		if (heterogeneous_mechanism_type_ == IBRAHIM_PAOLUCCI)
			nr_ = 4;

		// Memory allocation
		gamma_knudsen_.resize(ns_);
		gamma_knudsen_effective_.resize(ns_);
		gamma_fick_.resize(ns_);
		gamma_fick_effective_.resize(ns_);
		gamma_effective_.resize(ns_);
		Rgas_.resize(ns_);
		r_.resize(nr_);
		r_deposition_per_unit_area_per_single_reaction_.resize(nr_);
		r_deposition_per_unit_volume_per_single_reaction_.resize(nr_);

		// Set initial porosity
		SetPorosity(epsilon0_);

		// Indices of relevant species
		index_CH4_ = thermodynamicsMap_.IndexOfSpecies("CH4") - 1;
		index_C2H4_ = thermodynamicsMap_.IndexOfSpecies("C2H4") - 1;
		index_C2H2_ = thermodynamicsMap_.IndexOfSpecies("C2H2") - 1;
		index_C6H6_ = thermodynamicsMap_.IndexOfSpecies("C6H6") - 1;
		index_H2_ = thermodynamicsMap_.IndexOfSpecies("H2") - 1;

		// Set default hydrogen inhibition coefficients
		I_CH4_ = 1.;
		I_C2H4_ = 1.;
		I_C2H2_ = 1.;
		I_C6H6_ = 1.;

		// Summary on the screen
		Summary();
	}

	void PorousMedium::Summary()
	{
		std::cout << std::endl;
		std::cout << "-------------------------------------------------------" << std::endl;
		std::cout << "                      Porous Medium                    " << std::endl;
		std::cout << "-------------------------------------------------------" << std::endl;
		std::cout << " * Graphite density [kg/m3]:         " << rho_graphite_ << std::endl;
		std::cout << " * Bulk density density [kg/m3]:     " << density_bulk() << std::endl;
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

	void PorousMedium::SetGraphiteDensity(const double rho_graphite)
	{
		rho_graphite_ = rho_graphite;
	}

	double PorousMedium::density_bulk()
	{
		return rho_fiber_*(1. - epsilon0_) + rho_graphite_*(epsilon0_ - epsilon_);
	}

	double PorousMedium::Sv()
	{
		if (epsilon_ < 1e-3)
			return 1.;

		if (porous_substrate_type_ == POLYNOMIAL)
		{
			const double ratio = epsilon_ / epsilon0_;
			return 2. / rf_*((2. - epsilon0_)*ratio - std::pow(ratio, 2.));
		}
		else if (porous_substrate_type_ == RANDOM)
		{
			return -2. / rf_*epsilon_*std::log(epsilon_);
		}
		else if (porous_substrate_type_ == RANDOM_HARDCORE)
		{
			const double ratio = epsilon_ / epsilon0_;
			return 2. / rf_*ratio*(1. - epsilon_)*(1. / epsilon0_*std::log(1. / ratio) + 1.);
		}
		else if (porous_substrate_type_ == POLINOMIAL_ONEHALF)
		{
			const double ratio = epsilon_ / epsilon0_;
			return 2. / rf_*((2.-3./2.*epsilon_)*ratio-(1.-epsilon0_/2.)*ratio*ratio);
		}
		else if (porous_substrate_type_ == FROM_SPHERES_TO_CYLINDERS)
		{
			const double c = 3. / 4. / PhysicalConstants::pi;
			return 2. / rf_*(1. - epsilon_)*(std::pow(-c*std::log(1 - epsilon_), 2. / 3.) + epsilon_) / (std::pow(-c*std::log(1 - epsilon0_), 2. / 3.) + epsilon0_);
		}
		else
			return 0.;
	}

	double PorousMedium::rp()
	{
		return 2. / Sv()*epsilon_;
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
			gamma_knudsen_(i) = coefficient/std::sqrt(thermodynamicsMap_.MW()[i+1]);
	}

	void PorousMedium::KnudsenEffectiveDiffusionCoefficients()
	{
		KnudsenDiffusionCoefficients();

		const double correction = epsilon_ / eta_knudsen();
		gamma_knudsen_effective_ = gamma_knudsen_*correction;
	}

	void PorousMedium::FickDiffusionCoefficients(const Eigen::VectorXd& mole_fractions)
	{
		// Set temperature and pressure of transort map
		transportMap_.SetTemperature(T_);
		transportMap_.SetPressure(P_Pa_);

		// Fick diffusion coefficients
		{
			OpenSMOKE::OpenSMOKEVectorDouble gamma_fick_os(ns_);
			OpenSMOKE::OpenSMOKEVectorDouble mole_fractions_os(ns_);

			for (unsigned int i = 0; i < ns_; i++)
				mole_fractions_os[i + 1] = mole_fractions(i);
			transportMap_.MassDiffusionCoefficients(gamma_fick_os, mole_fractions_os, false);
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
	}

	void PorousMedium::HydrogenInhibitionCoefficients(const Eigen::VectorXd& C)
	{
		const double CH4  = C(index_CH4_);
		const double C2H4 = C(index_C2H4_);
		const double C2H2 = C(index_C2H2_);
		const double C6H6 = C(index_C6H6_);
		const double H2   = C(index_H2_);

		const double eps = 1.e-16;
		
		HeterogeneousMechanism ;
		if (hydrogen_inhibition_type_ == BECKER)
		{
			I_CH4_ = CH4  > eps ? 1.112 / (1.112 + H2 / CH4) : 0.;
			I_C2H4_ = C2H4 > eps ? 1.104 / (1.104 + H2 / C2H4) : 0.;
			I_C2H2_ = C2H2 > eps ? 3.269 / (3.269 + H2 / C2H2) : 0.;
			I_C6H6_ = C6H6 > eps ? 0.589 / (0.589 + H2 / C6H6) : 0.;
		}
	}

	void PorousMedium::FormationRates(const Eigen::VectorXd& C)
	{
		// Area per unit of volume [1/m]
		const double Sv_ = Sv();

		// Hydrogen inhibition factors
		if (hydrogen_inhibition_type_ != NONE)
			HydrogenInhibitionCoefficients(C);

		// Heterogeneous reactions
		if (heterogeneous_mechanism_type_ == IBRAHIM_PAOLUCCI)
		{
			const double CH4  = C(index_CH4_);
			const double C2H4 = C(index_C2H4_);
			const double C2H2 = C(index_C2H2_);
			const double C6H6 = C(index_C6H6_);
			const double H2   = C(index_H2_);

			// Frequency factors are in [m/s]
			const double kCH4  = 0.;
			const double kC2H4 = 72.4*exp(-155000 / PhysicalConstants::R_J_mol / T_);
			const double kC2H2 = 13.5*exp(-126000 / PhysicalConstants::R_J_mol / T_);
			const double kC6H6 = 4.75e5*exp(-217000 / PhysicalConstants::R_J_mol / T_);

			// Reaction rates (without corrections) [kmol/m2/s]
			const double rCH4  = kCH4*CH4;
			const double rC2H4 = kC2H4*C2H4;
			const double rC2H2 = kC2H2*C2H2;
			const double rC6H6 = kC6H6*C6H6;

			// Reaction rates [kmol/m2/s]
			r_(0) = rCH4*I_CH4_;
			r_(1) = rC2H4*I_C2H4_;
			r_(2) = rC2H2*I_C2H2_;
			r_(3) = rC6H6*I_C6H6_;

			// Consumption rate of homogeneous species [kmol/m3/s]
			Rgas_.setZero();
			Rgas_(index_CH4_)  = -Sv_*r_(0);
			Rgas_(index_C2H4_) = -Sv_*r_(1);
			Rgas_(index_C2H2_) = -Sv_*r_(2);
			Rgas_(index_C6H6_) = -Sv_*r_(3);
			Rgas_(index_H2_)   =  Sv_*(2.*r_(0) + 2.*r_(1) + 1.*r_(2) + 3.*r_(3));
		}

		// Heterogeneous deposition rate [kmol/m2/s]
		r_deposition_per_unit_area_per_single_reaction_(0) = r_(0);
		r_deposition_per_unit_area_per_single_reaction_(1) = 2.*r_(1);
		r_deposition_per_unit_area_per_single_reaction_(2) = 2.*r_(2);
		r_deposition_per_unit_area_per_single_reaction_(3) = 6.*r_(3);

		// Heterogeneous deposition rate [kmol/m3/s]
		r_deposition_per_unit_volume_per_single_reaction_ = Sv_*r_deposition_per_unit_area_per_single_reaction_;
	
		// Heterogeneous deposition rate [kmol/m2/s]
		r_deposition_per_unit_area_ = r_(0) + 2.*r_(1) + 2.*r_(2) + 6.*r_(3);

		// Heterogeneous deposition rate [kmol/m3/s]
		r_deposition_per_unit_volume_ = Sv_*r_deposition_per_unit_area_;
	}
}