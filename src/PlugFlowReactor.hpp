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
										OpenSMOKE::KineticsMap_CHEMKIN<double>& kineticsMap) :

	thermodynamicsMap_(thermodynamicsMap),
	kineticsMap_(kineticsMap)
	{
		n_steps_video_ = 10;
		count_video_ = n_steps_video_;

		ns_ = thermodynamicsMap_.NumberOfSpecies();
		ne_ = ns_;

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

	void PlugFlowReactor::SetInitialConditions(const double T_gas, const double P_gas, const OpenSMOKE::OpenSMOKEVectorDouble& omega_gas)
	{
		for (unsigned int j = 0; j < ns_; j++)
			Y_(j) = omega_gas[j + 1];

		P_ = P_gas;
		T_ = T_gas;
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
		// Species
		for (unsigned int j = 0; j < ns_; j++)
			dY_over_dt_(j) = aux_R[j+1]/rho_;
	}

	void PlugFlowReactor::Recover_Unknowns(const double* y)
	{
		unsigned int count = 0;

		// Species
		for (unsigned int j = 0; j < ns_; j++)
			Y_(j) = y[count++];
	}

	void PlugFlowReactor::Recover_Residuals(double* dy)
	{
		unsigned int count = 0;

		// Species
		for (unsigned int j = 0; j < ns_; j++)
			dy[count++] = dY_over_dt_(j);	
	}

	void PlugFlowReactor::UnknownsVector(double* v)
	{
		unsigned int count = 0;

		// Species
		for (unsigned int j = 0; j < ns_; j++)
			v[count++] = Y_(j);
	}

	void PlugFlowReactor::CorrectedUnknownsVector(double* v)
	{
		unsigned int count = 0;

		// Species
		for (unsigned int j = 0; j < ns_; j++)
			Y_(j) = v[count++];

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
		if (count_video_ == n_steps_video_)
		{
			std::cout << std::left << std::setw(16) << std::scientific << t;
			std::cout << std::endl;

			count_video_ = 0;
		}

		count_video_++;
	}
}