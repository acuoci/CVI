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
	HeterogeneousMechanism::HeterogeneousMechanism(	OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>& thermodynamicsMap,
								OpenSMOKE::KineticsMap_CHEMKIN<double>& kineticsMap,
								OpenSMOKE::TransportPropertiesMap_CHEMKIN<double>& transportMap,
								const bool homogeneous_reactions, const bool heterogeneous_reactions, 
								const HeterogeneousMechanismType heterogeneous_mechanism_type,
								const HydrogenInhibitionType hydrogen_inhibition_type) :

	thermodynamicsMap_(thermodynamicsMap),
	kineticsMap_(kineticsMap),
	transportMap_(transportMap),
	homogeneous_reactions_(homogeneous_reactions),
	heterogeneous_reactions_(heterogeneous_reactions),
	heterogeneous_mechanism_type_(heterogeneous_mechanism_type),
	hydrogen_inhibition_type_(hydrogen_inhibition_type)
	{
		Initialize();
	}

	HeterogeneousMechanism::HeterogeneousMechanism(	OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>& thermodynamicsMap,
								OpenSMOKE::KineticsMap_CHEMKIN<double>& kineticsMap,
								OpenSMOKE::TransportPropertiesMap_CHEMKIN<double>& transportMap,
								OpenSMOKE::OpenSMOKE_Dictionary& dictionary ) :
	thermodynamicsMap_(thermodynamicsMap),
	kineticsMap_(kineticsMap),
	transportMap_(transportMap)
	{
		// Read homogeneous reactions true/false
		{
			if (dictionary.CheckOption("@HomogeneousReactions") == true)
				dictionary.ReadBool("@HomogeneousReactions", homogeneous_reactions_);
		}

		// Read heterogeneous reactions true/false
		{
			if (dictionary.CheckOption("@HeterogeneousReactions") == true)
				dictionary.ReadBool("@HeterogeneousReactions", heterogeneous_reactions_);
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
				else OpenSMOKE::FatalErrorMessage("@GraphiteDensity: Unknown graphite density units. Available units: kg/m3 | g/cm3");
			}
		}

		// Read initial porosity
		heterogeneous_reaction_rates_multiplier_ = 1.;
		if (dictionary.CheckOption("@HeterogeneousReactionRatesMultiplier") == true)
			dictionary.ReadDouble("@HeterogeneousReactionRatesMultiplier", heterogeneous_reaction_rates_multiplier_);

		// Read heterogeneous mechanism type
		{
			std::string value;
			if (dictionary.CheckOption("@HeterogeneousMechanism") == true)
			{
				dictionary.ReadString("@HeterogeneousMechanism", value);
				if (value == "Huttinger")					heterogeneous_mechanism_type_ = CVI::HUTTINGER;
				else if (value == "Ziegler")				heterogeneous_mechanism_type_ = CVI::ZIEGLER;
				else if (value == "Vignoles")				heterogeneous_mechanism_type_ = CVI::VIGNOLES;
				else if (value == "Huttinger-extended")		heterogeneous_mechanism_type_ = CVI::HUTTINGER_EXTENDED;
				else if (value == "Ziegler-extended")		heterogeneous_mechanism_type_ = CVI::ZIEGLER_EXTENDED;
				else if (value == "Vignoles-extended")		heterogeneous_mechanism_type_ = CVI::VIGNOLES_EXTENDED;
				else OpenSMOKE::FatalErrorMessage("Heterogeneous mechanisms available: Huttinger | Ziegler | Vignoles |  Huttinger-extended | Ziegler-extended | Vignoles-extended");
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

	void HeterogeneousMechanism::Initialize()
	{
		// Constants
		mw_carbon_ = 12.010999679565430;	// [kg/kmol]

		// Number of species
		ns_ = thermodynamicsMap_.NumberOfSpecies();

		// Number of heterogeneous reactions
		if (heterogeneous_mechanism_type_ == CVI::HUTTINGER)
			nr_ = 4;
		else if (heterogeneous_mechanism_type_ == CVI::ZIEGLER)
			nr_ = 4;
		else if (heterogeneous_mechanism_type_ == CVI::VIGNOLES)
			nr_ = 4;
		else if (heterogeneous_mechanism_type_ == CVI::HUTTINGER_EXTENDED)
			nr_ = 5;
		else if (heterogeneous_mechanism_type_ == CVI::ZIEGLER_EXTENDED)
			nr_ = 5;
		else if (heterogeneous_mechanism_type_ == CVI::VIGNOLES_EXTENDED)
			nr_ = 5;

		// Memory allocation
		Rgas_.resize(ns_);
		r_.resize(nr_);
		r_deposition_per_unit_area_per_single_reaction_.resize(nr_);
		r_deposition_per_unit_volume_per_single_reaction_.resize(nr_);
		tags_.resize(nr_);

		// Indices of relevant species
		if (heterogeneous_mechanism_type_ == CVI::HUTTINGER)
		{
			index_CH4_ = thermodynamicsMap_.IndexOfSpecies("CH4") - 1;
			index_C2H4_ = thermodynamicsMap_.IndexOfSpecies("C2H4") - 1;
			index_C2H2_ = thermodynamicsMap_.IndexOfSpecies("C2H2") - 1;
			index_C6H6_ = thermodynamicsMap_.IndexOfSpecies("C6H6") - 1;
			index_C10H8_ = thermodynamicsMap_.IndexOfSpecies("C10H8") - 1;
			index_H2_ = thermodynamicsMap_.IndexOfSpecies("H2") - 1;

			tags_[0] = "CH4";	tags_[1] = "C2H4";	tags_[2] = "C2H2";	tags_[3] = "C6H6";
		}
		else if (heterogeneous_mechanism_type_ == CVI::ZIEGLER)
		{
			index_CH4_ = thermodynamicsMap_.IndexOfSpecies("CH4") - 1;
			index_C2H4_ = thermodynamicsMap_.IndexOfSpecies("C2H4") - 1;
			index_C2H2_ = thermodynamicsMap_.IndexOfSpecies("C2H2") - 1;
			index_C14H10_ = thermodynamicsMap_.IndexOfSpecies("C14H10") - 1;
			index_C10H8_ = thermodynamicsMap_.IndexOfSpecies("C10H8") - 1;
			index_H2_ = thermodynamicsMap_.IndexOfSpecies("H2") - 1;

			tags_[0] = "CH4";	tags_[1] = "C2H4";	tags_[2] = "C2H2";	tags_[3] = "C14H10";
		}
		else if (heterogeneous_mechanism_type_ == CVI::VIGNOLES)
		{
			index_CH4_ = thermodynamicsMap_.IndexOfSpecies("CH4") - 1;
			index_C6H6_ = thermodynamicsMap_.IndexOfSpecies("C6H6") - 1;
			index_C2H2_ = thermodynamicsMap_.IndexOfSpecies("C2H2") - 1;
			index_C14H10_ = thermodynamicsMap_.IndexOfSpecies("C14H10") - 1;
			index_C10H8_ = thermodynamicsMap_.IndexOfSpecies("C10H8") - 1;
			index_H2_ = thermodynamicsMap_.IndexOfSpecies("H2") - 1;

			tags_[0] = "CH4";	tags_[1] = "C6H6";	tags_[2] = "C2H2";	tags_[3] = "C14H10";
		}
		else if (heterogeneous_mechanism_type_ == CVI::HUTTINGER_EXTENDED)
		{
			index_CH4_ = thermodynamicsMap_.IndexOfSpecies("CH4") - 1;
			index_C2H4_ = thermodynamicsMap_.IndexOfSpecies("C2H4") - 1;
			index_C2H2_ = thermodynamicsMap_.IndexOfSpecies("C2H2") - 1;
			index_C6H6_ = thermodynamicsMap_.IndexOfSpecies("C6H6") - 1;
			index_C10H8_ = thermodynamicsMap_.IndexOfSpecies("C10H8") - 1;
			index_H2_ = thermodynamicsMap_.IndexOfSpecies("H2") - 1;

			tags_[0] = "CH4";	tags_[1] = "C2H4";	tags_[2] = "C2H2";	tags_[3] = "C6H6";	tags_[4] = "C10H8";
		}
		else if (heterogeneous_mechanism_type_ == CVI::ZIEGLER_EXTENDED)
		{
			index_CH4_ = thermodynamicsMap_.IndexOfSpecies("CH4") - 1;
			index_C2H4_ = thermodynamicsMap_.IndexOfSpecies("C2H4") - 1;
			index_C2H2_ = thermodynamicsMap_.IndexOfSpecies("C2H2") - 1;
			index_C14H10_ = thermodynamicsMap_.IndexOfSpecies("C14H10") - 1;
			index_C10H8_ = thermodynamicsMap_.IndexOfSpecies("C10H8") - 1;
			index_H2_ = thermodynamicsMap_.IndexOfSpecies("H2") - 1;

			tags_[0] = "CH4";	tags_[1] = "C2H4";	tags_[2] = "C2H2";	tags_[3] = "C14H10";	tags_[4] = "C10H8";
		}
		else if (heterogeneous_mechanism_type_ == CVI::VIGNOLES_EXTENDED)
		{
			index_CH4_ = thermodynamicsMap_.IndexOfSpecies("CH4") - 1;
			index_C6H6_ = thermodynamicsMap_.IndexOfSpecies("C6H6") - 1;
			index_C2H2_ = thermodynamicsMap_.IndexOfSpecies("C2H2") - 1;
			index_C14H10_ = thermodynamicsMap_.IndexOfSpecies("C14H10") - 1;
			index_C10H8_ = thermodynamicsMap_.IndexOfSpecies("C10H8") - 1;
			index_H2_ = thermodynamicsMap_.IndexOfSpecies("H2") - 1;

			tags_[0] = "CH4";	tags_[1] = "C6H6";	tags_[2] = "C2H2";	tags_[3] = "C14H10";	tags_[4] = "C10H8";
		}

		// Set default hydrogen inhibition coefficients
		I_CH4_ = 1.;
		I_C2H4_ = 1.;
		I_C2H2_ = 1.;
		I_C6H6_ = 1.;
		I_C10H8_ = 1.;
		I_C14H10_ = 1.;

		// Summary on the screen
		Summary();
	}

	void HeterogeneousMechanism::Summary()
	{
		std::cout << std::endl;
		std::cout << "-------------------------------------------------------" << std::endl;
		std::cout << "                      Porous Medium                    " << std::endl;
		std::cout << "-------------------------------------------------------" << std::endl;
		std::cout << " * Graphite density [kg/m3]:         " << rho_graphite_ << std::endl;
		std::cout << "--------------------------------------------------------" << std::endl;
		std::cout << std::endl;
	}

	void HeterogeneousMechanism::SetTemperature(const double T)
	{
		T_ = T;
	}

	void HeterogeneousMechanism::SetPressure(const double P_Pa)
	{
		P_Pa_ = P_Pa;
	}
	
	void HeterogeneousMechanism::SetGraphiteDensity(const double rho_graphite)
	{
		rho_graphite_ = rho_graphite;
	}

	void HeterogeneousMechanism::HydrogenInhibitionCoefficients(const Eigen::VectorXd& C)
	{
		const double CH4  = C(index_CH4_);
		const double C2H4 = C(index_C2H4_);
		const double C2H2 = C(index_C2H2_);
		const double C6H6 = C(index_C6H6_);
		const double C10H8 = C(index_C10H8_);
		const double H2   = C(index_H2_);

		const double eps = 1.e-16;
		
		if (hydrogen_inhibition_type_ == BECKER)
		{
			I_CH4_ = CH4  > eps ? 1.112 / (1.112 + H2 / CH4) : 0.;
			I_C2H4_ = C2H4 > eps ? 1.104 / (1.104 + H2 / C2H4) : 0.;
			I_C2H2_ = C2H2 > eps ? 3.269 / (3.269 + H2 / C2H2) : 0.;
			I_C6H6_ = C6H6 > eps ? 0.589 / (0.589 + H2 / C6H6) : 0.;
			I_C10H8_ = C10H8 > eps ? 0.589 / (0.589 + H2 / C10H8) : 0.;
			I_C14H10_ = 1.;
		}
	}

	void HeterogeneousMechanism::FormationRates(const double Sv, const Eigen::VectorXd& C)
	{
		// Hydrogen inhibition factors
		if (hydrogen_inhibition_type_ != NONE)
			HydrogenInhibitionCoefficients(C);

		// Heterogeneous reactions
		if (heterogeneous_mechanism_type_ == HUTTINGER)
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

			// Correction (if any) 
			r_ *= heterogeneous_reaction_rates_multiplier_;

			// Consumption rate of homogeneous species [kmol/m3/s]
			Rgas_.setZero();
			Rgas_(index_CH4_)  = -Sv*r_(0);
			Rgas_(index_C2H4_) = -Sv*r_(1);
			Rgas_(index_C2H2_) = -Sv*r_(2);
			Rgas_(index_C6H6_) = -Sv*r_(3);
			Rgas_(index_H2_)   =  Sv*(2.*r_(0) + 2.*r_(1) + 1.*r_(2) + 3.*r_(3));

			// Heterogeneous deposition rate [kmol/m2/s]
			r_deposition_per_unit_area_per_single_reaction_(0) = r_(0);
			r_deposition_per_unit_area_per_single_reaction_(1) = 2.*r_(1);
			r_deposition_per_unit_area_per_single_reaction_(2) = 2.*r_(2);
			r_deposition_per_unit_area_per_single_reaction_(3) = 6.*r_(3);

			// Heterogeneous deposition rate [kmol/m2/s]
			r_deposition_per_unit_area_ = r_(0) + 2.*r_(1) + 2.*r_(2) + 6.*r_(3);
		}
		else if (heterogeneous_mechanism_type_ == HUTTINGER_EXTENDED)
		{
			const double CH4 = C(index_CH4_);
			const double C2H4 = C(index_C2H4_);
			const double C2H2 = C(index_C2H2_);
			const double C6H6 = C(index_C6H6_);
			const double C10H8 = C(index_C10H8_);
			const double H2 = C(index_H2_);

			// Frequency factors are in [m/s]
			const double kCH4 = 0.;
			const double kC2H4 = 72.4*exp(-155000 / PhysicalConstants::R_J_mol / T_);
			const double kC2H2 = 13.5*exp(-126000 / PhysicalConstants::R_J_mol / T_);
			const double kC6H6 = 4.75e5*exp(-217000 / PhysicalConstants::R_J_mol / T_);
			const double kC10H8 = 4.75e5*exp(-217000 / PhysicalConstants::R_J_mol / T_);

			// Reaction rates (without corrections) [kmol/m2/s]
			const double rCH4 = kCH4*CH4;
			const double rC2H4 = kC2H4*C2H4;
			const double rC2H2 = kC2H2*C2H2;
			const double rC6H6 = kC6H6*C6H6;
			const double rC10H8 = kC10H8*C10H8;

			// Reaction rates [kmol/m2/s]
			r_(0) = rCH4*I_CH4_;
			r_(1) = rC2H4*I_C2H4_;
			r_(2) = rC2H2*I_C2H2_;
			r_(3) = rC6H6*I_C6H6_;
			r_(4) = rC10H8*I_C10H8_;

			// Correction (if any) 
			r_ *= heterogeneous_reaction_rates_multiplier_;

			// Consumption rate of homogeneous species [kmol/m3/s]
			Rgas_.setZero();
			Rgas_(index_CH4_) = -Sv*r_(0);
			Rgas_(index_C2H4_) = -Sv*r_(1);
			Rgas_(index_C2H2_) = -Sv*r_(2);
			Rgas_(index_C6H6_) = -Sv*r_(3);
			Rgas_(index_C10H8_) = -Sv*r_(4);
			Rgas_(index_H2_) = Sv*(2.*r_(0) + 2.*r_(1) + 1.*r_(2) + 3.*r_(3) + 4.*r_(4));

			// Heterogeneous deposition rate [kmol/m2/s]
			r_deposition_per_unit_area_per_single_reaction_(0) = r_(0);
			r_deposition_per_unit_area_per_single_reaction_(1) = 2.*r_(1);
			r_deposition_per_unit_area_per_single_reaction_(2) = 2.*r_(2);
			r_deposition_per_unit_area_per_single_reaction_(3) = 6.*r_(3);
			r_deposition_per_unit_area_per_single_reaction_(4) = 10.*r_(4);

			// Heterogeneous deposition rate [kmol/m2/s]
			r_deposition_per_unit_area_ = r_(0) + 2.*r_(1) + 2.*r_(2) + 6.*r_(3) + 10.*r_(4);
		}

		else if (heterogeneous_mechanism_type_ == ZIEGLER)
		{
			const double CH4 = C(index_CH4_);
			const double C2H4 = C(index_C2H4_);
			const double C2H2 = C(index_C2H2_);
			const double C14H10 = C(index_C14H10_);
			const double H2 = C(index_H2_);

			// Frequency factors are in [1/s]
			const double kCH4 = 0.;
			const double kC2H4 = 7.7e10*exp(-281000 / PhysicalConstants::R_J_mol / T_);
			const double kC2H2 = 4.1e10*exp(-281000 / PhysicalConstants::R_J_mol / T_);
			const double kC14H10 = 3.5e10*exp(-192000 / PhysicalConstants::R_J_mol / T_);

			// Reaction rates (without corrections) [kmol/m2/s]
			const double rCH4 = kCH4 / Sv * CH4;
			const double rC2H4 = kC2H4 / Sv * C2H4;
			const double rC2H2 = kC2H2 / Sv * C2H2;
			const double rC14H10 = kC14H10 / Sv * C14H10;

			// Reaction rates [kmol/m2/s]
			r_(0) = rCH4*I_CH4_;
			r_(1) = rC2H4*I_C2H4_;
			r_(2) = rC2H2*I_C2H2_;
			r_(3) = rC14H10*I_C14H10_;

			// Correction (if any) 
			r_ *= heterogeneous_reaction_rates_multiplier_;

			// Consumption rate of homogeneous species [kmol/m3/s]
			Rgas_.setZero();
			Rgas_(index_CH4_) = -Sv*r_(0);
			Rgas_(index_C2H4_) = -Sv*r_(1);
			Rgas_(index_C2H2_) = -Sv*r_(2);
			Rgas_(index_C14H10_) = -Sv*r_(3);
			Rgas_(index_H2_) = Sv*(2.*r_(0) + 2.*r_(1) + 1.*r_(2) + 5.*r_(3));

			// Heterogeneous deposition rate [kmol/m2/s]
			r_deposition_per_unit_area_per_single_reaction_(0) = r_(0);
			r_deposition_per_unit_area_per_single_reaction_(1) = 2.*r_(1);
			r_deposition_per_unit_area_per_single_reaction_(2) = 2.*r_(2);
			r_deposition_per_unit_area_per_single_reaction_(3) = 14.*r_(3);

			// Heterogeneous deposition rate [kmol/m2/s]
			r_deposition_per_unit_area_ = r_(0) + 2.*r_(1) + 2.*r_(2) + 14.*r_(3);
		}

		else if (heterogeneous_mechanism_type_ == ZIEGLER_EXTENDED)
		{
			const double CH4 = C(index_CH4_);
			const double C2H4 = C(index_C2H4_);
			const double C2H2 = C(index_C2H2_);
			const double C14H10 = C(index_C14H10_);
			const double C10H8 = C(index_C10H8_);
			const double H2 = C(index_H2_);

			// Frequency factors are in [1/s]
			const double kCH4 = 0.;
			const double kC2H4 = 7.7e10*exp(-281000 / PhysicalConstants::R_J_mol / T_);
			const double kC2H2 = 4.1e10*exp(-281000 / PhysicalConstants::R_J_mol / T_);
			const double kC14H10 = 3.5e10*exp(-192000 / PhysicalConstants::R_J_mol / T_);
			const double kC10H8 = 3.5e10*exp(-192000 / PhysicalConstants::R_J_mol / T_);

			// Reaction rates (without corrections) [kmol/m2/s]
			const double rCH4 = kCH4 / Sv * CH4;
			const double rC2H4 = kC2H4 / Sv * C2H4;
			const double rC2H2 = kC2H2 / Sv * C2H2;
			const double rC14H10 = kC14H10 / Sv * C14H10;
			const double rC10H8 = kC10H8 / Sv * C10H8;

			// Reaction rates [kmol/m2/s]
			r_(0) = rCH4*I_CH4_;
			r_(1) = rC2H4*I_C2H4_;
			r_(2) = rC2H2*I_C2H2_;
			r_(3) = rC14H10*I_C14H10_;
			r_(4) = rC10H8*I_C10H8_;

			// Correction (if any) 
			r_ *= heterogeneous_reaction_rates_multiplier_;

			// Consumption rate of homogeneous species [kmol/m3/s]
			Rgas_.setZero();
			Rgas_(index_CH4_) = -Sv*r_(0);
			Rgas_(index_C2H4_) = -Sv*r_(1);
			Rgas_(index_C2H2_) = -Sv*r_(2);
			Rgas_(index_C14H10_) = -Sv*r_(3);
			Rgas_(index_C10H8_) = -Sv*r_(4);
			Rgas_(index_H2_) = Sv*(2.*r_(0) + 2.*r_(1) + 1.*r_(2) + 5.*r_(3) + 4.*r_(4));

			// Heterogeneous deposition rate [kmol/m2/s]
			r_deposition_per_unit_area_per_single_reaction_(0) = r_(0);
			r_deposition_per_unit_area_per_single_reaction_(1) = 2.*r_(1);
			r_deposition_per_unit_area_per_single_reaction_(2) = 2.*r_(2);
			r_deposition_per_unit_area_per_single_reaction_(3) = 14.*r_(3);
			r_deposition_per_unit_area_per_single_reaction_(4) = 10.*r_(4);

			// Heterogeneous deposition rate [kmol/m2/s]
			r_deposition_per_unit_area_ = r_(0) + 2.*r_(1) + 2.*r_(2) + 14.*r_(3) + 10.*r_(4);
		}

		else if (heterogeneous_mechanism_type_ == VIGNOLES)
		{
			const double CH4 = C(index_CH4_);
			const double C2H2 = C(index_C2H2_);
			const double C6H6 = C(index_C6H6_);
			const double C14H10 = C(index_C14H10_);
			const double H2 = C(index_H2_);

			// Frequency factors are in [m/s]
			const double kCH4 = 0.;
			const double kC2H2 = 31.4*exp(-120000 / PhysicalConstants::R_J_mol / T_);
			const double kC6H6 = 1.41e4*exp(-120000 / PhysicalConstants::R_J_mol / T_);
			const double kC14H10 = 7e7*exp(-230000 / PhysicalConstants::R_J_mol / T_);

			// Reaction rates (without corrections) [kmol/m2/s]
			const double rCH4 = kCH4*CH4;
			const double rC2H2 = kC2H2*C2H2;
			const double rC6H6 = kC6H6*C6H6;
			const double rC14H10 = kC14H10*C14H10;

			// Reaction rates [kmol/m2/s]
			r_(0) = rCH4*I_CH4_;
			r_(1) = rC2H2*I_C2H2_;
			r_(2) = rC6H6*I_C6H6_;
			r_(3) = rC14H10*I_C14H10_;

			// Correction (if any) 
			r_ *= heterogeneous_reaction_rates_multiplier_;

			// Consumption rate of homogeneous species [kmol/m3/s]
			Rgas_.setZero();
			Rgas_(index_CH4_) = -Sv*r_(0);
			Rgas_(index_C2H2_) = -Sv*r_(1);
			Rgas_(index_C6H6_) = -Sv*r_(2);
			Rgas_(index_C14H10_) = -Sv*r_(3);
			Rgas_(index_H2_) = Sv*(2.*r_(0) + 1.*r_(1) + 3.*r_(2) + 5.*r_(3));

			// Heterogeneous deposition rate [kmol/m2/s]
			r_deposition_per_unit_area_per_single_reaction_(0) = r_(0);
			r_deposition_per_unit_area_per_single_reaction_(1) = 2.*r_(1);
			r_deposition_per_unit_area_per_single_reaction_(2) = 6.*r_(2);
			r_deposition_per_unit_area_per_single_reaction_(3) = 14.*r_(3);

			// Heterogeneous deposition rate [kmol/m2/s]
			r_deposition_per_unit_area_ = r_(0) + 2.*r_(1) + 6.*r_(2) + 14.*r_(3);
		}

		else if (heterogeneous_mechanism_type_ == VIGNOLES_EXTENDED)
		{
			const double CH4 = C(index_CH4_);
			const double C2H2 = C(index_C2H2_);
			const double C6H6 = C(index_C6H6_);
			const double C14H10 = C(index_C14H10_);
			const double C10H8 = C(index_C10H8_);
			const double H2 = C(index_H2_);

			// Frequency factors are in [m/s]
			const double kCH4 = 0.;
			const double kC2H2 = 31.4*exp(-120000 / PhysicalConstants::R_J_mol / T_);
			const double kC6H6 = 1.41e4*exp(-120000 / PhysicalConstants::R_J_mol / T_);
			const double kC14H10 = 7e7*exp(-230000 / PhysicalConstants::R_J_mol / T_);
			const double kC10H8 = 7e7*exp(-230000 / PhysicalConstants::R_J_mol / T_);

			// Reaction rates (without corrections) [kmol/m2/s]
			const double rCH4 = kCH4*CH4;
			const double rC2H2 = kC2H2*C2H2;
			const double rC6H6 = kC6H6*C6H6;
			const double rC14H10 = kC14H10*C14H10;
			const double rC10H8 = kC10H8*C10H8;

			// Reaction rates [kmol/m2/s]
			r_(0) = rCH4*I_CH4_;
			r_(1) = rC2H2*I_C2H2_;
			r_(2) = rC6H6*I_C6H6_;
			r_(3) = rC14H10*I_C14H10_;
			r_(4) = rC10H8*I_C10H8_;

			// Correction (if any) 
			r_ *= heterogeneous_reaction_rates_multiplier_;

			// Consumption rate of homogeneous species [kmol/m3/s]
			Rgas_.setZero();
			Rgas_(index_CH4_) = -Sv*r_(0);
			Rgas_(index_C2H2_) = -Sv*r_(1);
			Rgas_(index_C6H6_) = -Sv*r_(2);
			Rgas_(index_C14H10_) = -Sv*r_(3);
			Rgas_(index_C10H8_) = -Sv*r_(4);
			Rgas_(index_H2_) = Sv*(2.*r_(0) + 1.*r_(1) + 3.*r_(2) + 5.*r_(3) + 4.*r_(4));

			// Heterogeneous deposition rate [kmol/m2/s]
			r_deposition_per_unit_area_per_single_reaction_(0) = r_(0);
			r_deposition_per_unit_area_per_single_reaction_(1) = 2.*r_(1);
			r_deposition_per_unit_area_per_single_reaction_(2) = 6.*r_(2);
			r_deposition_per_unit_area_per_single_reaction_(3) = 14.*r_(3);
			r_deposition_per_unit_area_per_single_reaction_(4) = 10.*r_(4);

			// Heterogeneous deposition rate [kmol/m2/s]
			r_deposition_per_unit_area_ = r_(0) + 2.*r_(1) + 6.*r_(2) + 14.*r_(3) + 10.*r_(4);
		}
		
		// Heterogeneous deposition rate [kmol/m3/s]
		r_deposition_per_unit_volume_per_single_reaction_ = Sv*r_deposition_per_unit_area_per_single_reaction_;
	
		// Heterogeneous deposition rate [kmol/m3/s]
		r_deposition_per_unit_volume_ = Sv*r_deposition_per_unit_area_;
	}
}
