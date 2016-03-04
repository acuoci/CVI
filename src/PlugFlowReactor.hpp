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
#include "Interface_PlugFlowReactor_ODE.h"

namespace CVI
{
	PlugFlowReactor::PlugFlowReactor(	OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>& thermodynamicsMap,
										OpenSMOKE::KineticsMap_CHEMKIN<double>& kineticsMap,
										const double v, const double Dh) :

	thermodynamicsMap_(thermodynamicsMap),
	kineticsMap_(kineticsMap),
	v_(v),
	Dh_(Dh)
	{
		DefaultValues();
		Initialize();
	}

	PlugFlowReactor::PlugFlowReactor(	OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>& thermodynamicsMap,
										OpenSMOKE::KineticsMap_CHEMKIN<double>& kineticsMap,
										OpenSMOKE::OpenSMOKE_Dictionary& dictionary ) :

		thermodynamicsMap_(thermodynamicsMap),
		kineticsMap_(kineticsMap)
	{
		// Sets the default values
		DefaultValues();

		// Read the velocity (mandatory)
		{
			std::string units;
			if (dictionary.CheckOption("@Velocity") == true)
			{
				dictionary.ReadMeasure("@Velocity", v_, units);
				if (units == "m/s")				v_ = v_;
				else if (units == "cm/s")		v_ = v_ / 1.e2;
				else if (units == "mm/s")		v_ = v_ / 1.e3;
				else OpenSMOKE::FatalErrorMessage("Unknown @Velocity units");
			}
		}

		// Read the hydraulic diameter (mandatory)
		{
			std::string units;
			if (dictionary.CheckOption("@HydraulicDiameter") == true)
			{
				dictionary.ReadMeasure("@HydraulicDiameter", Dh_, units);
				if (units == "m")			Dh_ = Dh_;
				else if (units == "cm")		Dh_ = Dh_ / 1.e2;
				else if (units == "mm")		Dh_ = Dh_ / 1.e3;
				else OpenSMOKE::FatalErrorMessage("Unknown@HydraulicDiameter units");
			}
		}

		// Read the inlet length (optional)
		{
			std::string units;
			if (dictionary.CheckOption("@InletLength") == true)
			{
				dictionary.ReadMeasure("@InletLength", inert_length_, units);
				if (units == "m")			inert_length_ = inert_length_;
				else if (units == "cm")		inert_length_ = inert_length_ / 1.e2;
				else if (units == "mm")		inert_length_ = inert_length_ / 1.e3;
				else OpenSMOKE::FatalErrorMessage("Unknown @InletLength units");
			}
		}

		// Read the asymptotic Nusselt number (optional)
		{
			if (dictionary.CheckOption("@AsymptoticNusselt") == true)
				dictionary.ReadDouble("@AsymptoticNusselt", NuInf_);
		}

		// Read the internal boundary layer limitations (true/false) (optional)
		{
			if (dictionary.CheckOption("@InternalBoundaryLayer") == true)
				dictionary.ReadBool("@InternalBoundaryLayer", internal_boundary_layer_correction_);
		}

		// Read the geometric pattern (optional)
		{
			std::string value;
			if (dictionary.CheckOption("@GeometricPattern") == true)
			{
				dictionary.ReadString("@GeometricPattern", value);
				if (value == "OneSide")				geometric_pattern_ = CVI::PlugFlowReactor::ONE_SIDE;
				else if (value == "ThreeSides")		geometric_pattern_ = CVI::PlugFlowReactor::THREE_SIDES;
				else OpenSMOKE::FatalErrorMessage("@GeometricPattern available: OneSide | ThreeSides");
			}
		}

		// Initialize and allocate memory
		Initialize();
	}

	void PlugFlowReactor::DefaultValues()
	{
		inert_length_ = 0.;
		NuInf_ = 3.66;
		internal_boundary_layer_correction_ = false;
		geometric_pattern_ = ONE_SIDE;
	}
	
	void PlugFlowReactor::Initialize()
	{
		n_steps_video_ = 10;
		count_video_ = n_steps_video_;

		ns_ = thermodynamicsMap_.NumberOfSpecies();
		ne_ = ns_ + 1;

		output_folder_ = "OutputPFR";

		MemoryAllocation();
	}

	void PlugFlowReactor::MemoryAllocation()
	{
		OpenSMOKE::ChangeDimensions(ns_, &aux_X, true);
		OpenSMOKE::ChangeDimensions(ns_, &aux_Y, true);
		OpenSMOKE::ChangeDimensions(ns_, &aux_C, true);
		OpenSMOKE::ChangeDimensions(ns_, &aux_R, true);

		Y_.resize(ns_);
		dY_over_dt_.resize(ns_);
	}

	void PlugFlowReactor::SetInitialConditions(const double T_gas, const double P_gas, const Eigen::VectorXd& omega_gas)
	{
		for (unsigned int j = 0; j < ns_; j++)
			Y_(j) = omega_gas(j);

		P_ = P_gas;
		T_ = T_gas;

		csi_ = 0.;
	}

	void PlugFlowReactor::SetInertLength(const double inert_length)
	{
		inert_length_ = inert_length;
	}

	void PlugFlowReactor::SetAsymptoticNusseltNumber(const double NuInf)
	{
		NuInf_ = NuInf;
	}

	void PlugFlowReactor::SetInternalBoundaryLayerCorrection(const bool flag)
	{
		internal_boundary_layer_correction_ = flag;
	}

	void PlugFlowReactor::SetGeometricPattern(const GeometricPattern pattern)
	{
		geometric_pattern_ = pattern;
	}

	double PlugFlowReactor::mass_transfer_coefficient(const double T, const double P_Pa, const double rho, const double x)
	{
		const double logT = std::log(T);
		const double mu = std::exp(-2.179357e+001 + 3.339017e+000*logT - 3.532041e-001*logT*logT + 1.543465e-002*std::pow(logT, 3.));
		const double gamma = std::exp(-2.485390e+001 + 3.532615e+000*logT - 2.444656e-001*std::pow(logT, 2.) + 1.061160e-002*std::pow(logT, 3.)) / (P_Pa / 1e5);

		const double nuGas = mu / rho;
		const double ReGas = rho*v_*Dh_ / mu;
		const double ScGas = nuGas / gamma;

		const double zStar = 1. / (ReGas*ScGas)*x / Dh_;
		const double Sh = NuInf_ + 6.874*std::pow(1000.*zStar, -0.488)*std::exp(-57.2*zStar);
		const double kc = Sh*gamma / Dh_;

		return kc;
	}

	void PlugFlowReactor::Properties()
	{
		// Termodynamics
		thermodynamicsMap_.SetPressure(P_);
		thermodynamicsMap_.SetTemperature(T_);

		// Mole fractions
		double mw_;
		aux_Y.CopyFrom(Y_.data());
		thermodynamicsMap_.MoleFractions_From_MassFractions(aux_X, mw_, aux_Y);

		// Concentrations [kmol/m3]
		const double cTot = P_ / PhysicalConstants::R_J_kmol / T_; // [kmol/m3]
		Product(cTot, aux_X, &aux_C);

		// Mixture density
		rho_ = cTot*mw_;	// [kg/m3]
		
		// Homogeneous phase
		kineticsMap_.SetTemperature(T_);
		kineticsMap_.SetPressure(P_);
		kineticsMap_.ReactionRates(aux_C);
		kineticsMap_.FormationRates(&aux_R);
		ElementByElementProduct(aux_R, thermodynamicsMap_.MW(), &aux_R); // [kg/m3/s]
	}

	void PlugFlowReactor::SubEquations_MassFractions()
	{
		for (unsigned int j = 0; j < ns_; j++)
			dY_over_dt_(j) = aux_R[j+1]/rho_;
	}

	void PlugFlowReactor::SubEquations_Space()
	{
		dcsi_over_dt_ = v_;
	}

	void PlugFlowReactor::Recover_Unknowns(const double* y)
	{
		unsigned int count = 0;

		// Species
		for (unsigned int j = 0; j < ns_; j++)
			Y_(j) = y[count++];

		// Space
		csi_ = y[count++];
	}

	void PlugFlowReactor::Recover_Residuals(double* dy)
	{
		unsigned int count = 0;

		// Species
		for (unsigned int j = 0; j < ns_; j++)
			dy[count++] = dY_over_dt_(j);	

		// Space
		dy[count++] = dcsi_over_dt_;
	}

	void PlugFlowReactor::UnknownsVector(double* v)
	{
		unsigned int count = 0;

		// Species
		for (unsigned int j = 0; j < ns_; j++)
			v[count++] = Y_(j);

		// Space
		v[count++] = csi_;
	}

	void PlugFlowReactor::CorrectedUnknownsVector(double* v)
	{
		unsigned int count = 0;

		// Species
		for (unsigned int j = 0; j < ns_; j++)
			Y_(j) = v[count++];

		// Space
		csi_ = v[count++];

		// Normalize
		const double sum = Y_.sum();
		for (unsigned int j = 0; j < ns_; j++)
			Y_(j) /= sum;
	}

	void PlugFlowReactor::Equations(const double t, const double* y, double* dy)
	{
		// Recover unknowns
		Recover_Unknowns(y);

		// Properties
		Properties();

		// Equations
		SubEquations_MassFractions();
		SubEquations_Space();

		// Recover residuals
		Recover_Residuals(dy);
	}

	bool PlugFlowReactor::Solve(const double tau)
	{
		bool flag = PlugFlowReactor_Ode(this, 0., tau);
		return flag;
	}

	void PlugFlowReactor::PrintSolution(const std::string name_file)
	{	
	}

	void PlugFlowReactor::Print(const double t, const double* y)
	{
		// History
		{
			Eigen::VectorXd y_eigen(ne_-1);
			for (unsigned int i = 0; i < ne_-1; i++)
				y_eigen(i) = y[i];

			history_tau_.push_back(t);
			history_csi_.push_back(y[ne_-1]);
			history_Y_.push_back(y_eigen);
		}

		if (count_video_ == n_steps_video_)
		{
			std::cout << std::left << std::setw(16) << std::scientific << t;
			std::cout << std::left << std::setw(16) << std::scientific << csi_*1000.;
			std::cout << std::endl;

			count_video_ = 0;
		}

		count_video_++;
	}
}