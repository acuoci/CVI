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

// OpenSMOKE++ Definitions
#include "OpenSMOKEpp"

// CHEMKIN maps
#include "maps/Maps_CHEMKIN"

// Reactor utilities
#include "reactors/utilities/Utilities"

// 1D grid
#include "grids/adaptive/Grid1D.h"

// Grammar CVI_Solver1D
#include "Grammar_CVI_Solver1D.h"

// Porous medium
#include "PorousMedium.h"

// Plug Flow Reactor
#include "PlugFlowReactor.h"

// Reactor 1D
#include "Reactor1D.h"

// Numerical parameters
#include "premixedlaminarflame1d/parameters/DAE_Parameters.h"

int main(int argc, char** argv)
{
	boost::filesystem::path executable_file = OpenSMOKE::GetExecutableFileName(argv);
	boost::filesystem::path executable_folder = executable_file.parent_path();

	OpenSMOKE::OpenSMOKE_logo("CVI_Solver1D", "Alberto Cuoci (alberto.cuoci@polimi.it)");

	std::string input_file_name_ = "input.dic";
	std::string main_dictionary_name_ = "CVI_Solver1D";

	// Program options from command line
	{
		namespace po = boost::program_options;
		po::options_description description("Options for the CVI_Solver1D");
		description.add_options()
			("help", "print help messages")
			("input", po::value<std::string>(), "name of the file containing the main dictionary (default \"input.dic\")")
			("dictionary", po::value<std::string>(), "name of the main dictionary to be used (default \"PlugFlowReactor\")");

		po::variables_map vm;
		try
		{
			po::store(po::parse_command_line(argc, argv, description), vm); // can throw 

			if (vm.count("help"))
			{
				std::cout << "Basic Command Line Parameters" << std::endl;
				std::cout << description << std::endl;
				return OPENSMOKE_SUCCESSFULL_EXIT;
			}

			if (vm.count("input"))
				input_file_name_ = vm["input"].as<std::string>();

			if (vm.count("dictionary"))
				main_dictionary_name_ = vm["dictionary"].as<std::string>();

			po::notify(vm); // throws on error, so do after help in case  there are any problems 
		}
		catch (po::error& e)
		{
			std::cerr << "Fatal error: " << e.what() << std::endl << std::endl;
			std::cerr << description << std::endl;
			return OPENSMOKE_FATAL_ERROR_EXIT;
		}
	}

	// Defines the grammar rules
	CVI::Grammar_CVI_Solver1D grammar_cvi_solver1d;

	// Define the dictionaries
	OpenSMOKE::OpenSMOKE_DictionaryManager dictionaries;
	dictionaries.ReadDictionariesFromFile(input_file_name_);
	dictionaries(main_dictionary_name_).SetGrammar(grammar_cvi_solver1d);

	// Kinetic scheme
	boost::filesystem::path path_kinetics_output;
	if (dictionaries(main_dictionary_name_).CheckOption("@KineticsFolder") == true)
	{
		dictionaries(main_dictionary_name_).ReadPath("@KineticsFolder", path_kinetics_output);
		OpenSMOKE::CheckKineticsFolder(path_kinetics_output);
	}
	else
	{
		std::string name_of_rapid_kinetics_subdictionary;
		if (dictionaries(main_dictionary_name_).CheckOption("@KineticsPreProcessor") == true)
			dictionaries(main_dictionary_name_).ReadDictionary("@KineticsPreProcessor", name_of_rapid_kinetics_subdictionary);

		OpenSMOKE::Grammar_RapidKineticMechanism grammar_rapid_kinetics;
		dictionaries(name_of_rapid_kinetics_subdictionary).SetGrammar(grammar_rapid_kinetics);

		boost::filesystem::path path_input_thermodynamics;
		if (dictionaries(name_of_rapid_kinetics_subdictionary).CheckOption("@Thermodynamics") == true)
			dictionaries(name_of_rapid_kinetics_subdictionary).ReadPath("@Thermodynamics", path_input_thermodynamics);

		boost::filesystem::path path_input_kinetics;
		if (dictionaries(name_of_rapid_kinetics_subdictionary).CheckOption("@Kinetics") == true)
			dictionaries(name_of_rapid_kinetics_subdictionary).ReadPath("@Kinetics", path_input_kinetics);

		boost::filesystem::path path_input_transport;
		if (dictionaries(name_of_rapid_kinetics_subdictionary).CheckOption("@Transport") == true)
			dictionaries(name_of_rapid_kinetics_subdictionary).ReadPath("@Transport", path_input_transport);

		if (dictionaries(name_of_rapid_kinetics_subdictionary).CheckOption("@Output") == true)
			dictionaries(name_of_rapid_kinetics_subdictionary).ReadPath("@Output", path_kinetics_output);

		OpenSMOKE::RapidKineticMechanismWithTransport(path_kinetics_output, path_input_transport.c_str(), path_input_thermodynamics.c_str(), path_input_kinetics.c_str());
	}

	// Read thermodynamics and kinetics maps
	OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>*		thermodynamicsMapXML;
	OpenSMOKE::KineticsMap_CHEMKIN<double>*				kineticsMapXML;
	OpenSMOKE::TransportPropertiesMap_CHEMKIN<double>*	transportMapXML;

	{
		rapidxml::xml_document<> doc;
		std::vector<char> xml_string;
		OpenSMOKE::OpenInputFileXML(doc, xml_string, path_kinetics_output / "kinetics.xml");

		double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
		thermodynamicsMapXML = new OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>(doc);
		kineticsMapXML = new OpenSMOKE::KineticsMap_CHEMKIN<double>(*thermodynamicsMapXML, doc);
		transportMapXML = new OpenSMOKE::TransportPropertiesMap_CHEMKIN<double>(doc);

		// Transport properties
		transportMapXML = new OpenSMOKE::TransportPropertiesMap_CHEMKIN<double>(doc);

		double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
		std::cout << "Time to read XML file: " << tEnd - tStart << std::endl;
	}

	// Monodimensional grid
	OpenSMOKE::Grid1D* grid;
	{
		int number_of_points;
		if (dictionaries(main_dictionary_name_).CheckOption("@Points") == true)
			dictionaries(main_dictionary_name_).ReadInt("@Points", number_of_points);

		double stretching_factor = 1.;
		if (dictionaries(main_dictionary_name_).CheckOption("@StretchingFactor") == true)
			dictionaries(main_dictionary_name_).ReadDouble("@StretchingFactor", stretching_factor);

		double length;
		std::string units;
		if (dictionaries(main_dictionary_name_).CheckOption("@Length") == true)
		{
			dictionaries(main_dictionary_name_).ReadMeasure("@Length", length, units);
			if (units == "m")			length = length;
			else if (units == "cm")		length = length / 1.e2;
			else if (units == "mm")		length = length / 1.e3;
			else if (units == "micron")	length = length / 1.e6;
			else OpenSMOKE::FatalErrorMessage("Unknown length units");
		}

		Eigen::VectorXd x(number_of_points);
		x(0) = 0.;
		for (unsigned int i = 1; i < number_of_points; i++)
			x(i) = x(i - 1) + length / double(number_of_points - 1);
		
//		grid = new OpenSMOKE::Grid1D(number_of_points, 0., length, stretching_factor);
		grid = new OpenSMOKE::Grid1D(x);
	}

	// Dae Options
	OpenSMOKE::DAE_Parameters* dae_parameters;
	dae_parameters = new OpenSMOKE::DAE_Parameters();
	if (dictionaries(main_dictionary_name_).CheckOption("@DaeParameters") == true)
	{
		std::string name_of_subdictionary;
		dictionaries(main_dictionary_name_).ReadDictionary("@DaeParameters", name_of_subdictionary);
		dae_parameters->SetupFromDictionary(dictionaries(name_of_subdictionary));
	}

	// Read inlet conditions
	double inlet_T;
	double inlet_P;
	OpenSMOKE::OpenSMOKEVectorDouble inlet_omega;
	{
		std::string dict_name;
		if (dictionaries(main_dictionary_name_).CheckOption("@InletStream") == true)
			dictionaries(main_dictionary_name_).ReadDictionary("@InletStream", dict_name);
			GetGasStatusFromDictionary(dictionaries(dict_name), *thermodynamicsMapXML, inlet_T, inlet_P, inlet_omega);
	}

	// Read initial conditions (uniform)
	double initial_T;
	double initial_P;
	OpenSMOKE::OpenSMOKEVectorDouble initial_omega;
	{
		std::string dict_name;
		if (dictionaries(main_dictionary_name_).CheckOption("@InitialConditions") == true)
			dictionaries(main_dictionary_name_).ReadDictionary("@InitialConditions", dict_name);
		GetGasStatusFromDictionary(dictionaries(dict_name), *thermodynamicsMapXML, initial_T, initial_P, initial_omega);
	}

	// Read fiber radius
	double rf = 0.;
	{
		double value;
		std::string units;
		if (dictionaries(main_dictionary_name_).CheckOption("@FiberRadius") == true)
		{
			dictionaries(main_dictionary_name_).ReadMeasure("@FiberRadius", value, units);
			if (units == "m")				rf = value;
			else if (units == "cm")			rf = value / 1.e2;
			else if (units == "mm")			rf = value / 1.e3;
			else if (units == "micron")		rf = value / 1.e6;
			else OpenSMOKE::FatalErrorMessage("Unknown fiber units");
		}
	}

	// Read initial porosity
	double epsilon0 = 0.;
	if (dictionaries(main_dictionary_name_).CheckOption("@InitialPorosity") == true)
		dictionaries(main_dictionary_name_).ReadDouble("@InitialPorosity", epsilon0);

	// Read porous substrate type
	CVI::PorousSubstrateType porous_substrate_type;
	{
		std::string value;
		if (dictionaries(main_dictionary_name_).CheckOption("@PorousSubstrate") == true)
		{
			dictionaries(main_dictionary_name_).ReadString("@PorousSubstrate", value);
			if (value == "polynomial")			porous_substrate_type = CVI::POLYNOMIAL;
			else if (value == "random")		porous_substrate_type = CVI::RANDOM;
			else if (value == "random_hardcore")		porous_substrate_type = CVI::RANDOM_HARDCORE;
			else if (value == "polynomial_onehalf")		porous_substrate_type = CVI::POLINOMIAL_ONEHALF;
			else if (value == "from_spheres_to_cylinders")		porous_substrate_type = CVI::FROM_SPHERES_TO_CYLINDERS;
			else OpenSMOKE::FatalErrorMessage("Substrates available: polynomial | random | random_hardcore | polynomial_onehalf | from_spheres_to_cylinders");
		}
	}

	// Read heterogeneous mechanism type
	CVI::HeterogeneousMechanism heterogeneous_mechanism_type;
	{
		std::string value;
		if (dictionaries(main_dictionary_name_).CheckOption("@HeterogeneousMechanism") == true)
		{
			dictionaries(main_dictionary_name_).ReadString("@HeterogeneousMechanism", value);
			if (value == "Ibrahim-Paolucci")		heterogeneous_mechanism_type = CVI::IBRAHIM_PAOLUCCI;
			else OpenSMOKE::FatalErrorMessage("Heterogeneous mechanisms available: Ibrahim-Paolucci");
		}
	}

	// Read hydrogen inhibition type
	CVI::HydrogenInhibitionType hydrogen_inhibition_type;
	{
		std::string value;
		if (dictionaries(main_dictionary_name_).CheckOption("@HydrogenInhibition") == true)
		{
			dictionaries(main_dictionary_name_).ReadString("@HydrogenInhibition", value);
			if (value == "none")			hydrogen_inhibition_type = CVI::NONE;
			else if (value == "Becker")		hydrogen_inhibition_type = CVI::BECKER;
			else OpenSMOKE::FatalErrorMessage("Hydrogen inhibitions available: none | Becker");
		}
	}

	// Read the plug flow residence time
	double tau = 0.;
	{
		std::string units;
		if (dictionaries(main_dictionary_name_).CheckOption("@PlugFlowResidenceTime") == true)
		{
			dictionaries(main_dictionary_name_).ReadMeasure("@PlugFlowResidenceTime", tau, units);
			if (units == "s")				tau = tau;
			else if (units == "ms")			tau = tau / 1.e3;
			else if (units == "min")		tau = tau * 60.;
			else OpenSMOKE::FatalErrorMessage("Unknown @PlugFlowResidenceTime units");
		}
	}

	// Plug flow ractor simulation
	OpenSMOKE::OpenSMOKEVectorDouble pfr_omega(thermodynamicsMapXML->NumberOfSpecies());
	{
		CVI::PlugFlowReactor* plug_flow_reactor = new CVI::PlugFlowReactor(*thermodynamicsMapXML, *kineticsMapXML);
		plug_flow_reactor->SetInitialConditions(inlet_T, inlet_P, inlet_omega);

		plug_flow_reactor->Solve(tau);
		
		for (unsigned int i = 0; i < thermodynamicsMapXML->NumberOfSpecies(); i++)
			pfr_omega[i + 1] = plug_flow_reactor->Y()(i);
	}

	for (unsigned int i = 0; i < thermodynamicsMapXML->NumberOfSpecies(); i++)
		std::cout << thermodynamicsMapXML->NamesOfSpecies()[i] << " " << inlet_omega[i + 1] << " " << pfr_omega[i + 1] << std::endl;;

	getchar();

	CVI::PorousMedium* porous_medium = new CVI::PorousMedium(	*thermodynamicsMapXML, *kineticsMapXML, *transportMapXML, 
																porous_substrate_type, rf, epsilon0, 
																heterogeneous_mechanism_type, hydrogen_inhibition_type);

	CVI::Reactor1D* reactor1d = new CVI::Reactor1D(*thermodynamicsMapXML, *kineticsMapXML, *transportMapXML, *porous_medium, *grid);

	reactor1d->SetInitialConditions(initial_T, initial_P, initial_omega);
	reactor1d->SetGasSide(inlet_T, inlet_P, pfr_omega);

	// Solve
	{
		time_t timerStart;
		time_t timerEnd;

		time(&timerStart);
		reactor1d->SolveFromScratch(*dae_parameters);
		time(&timerEnd);

		std::cout << "Total time: " << difftime(timerEnd, timerStart) << " s" << std::endl;
	}
	
	/*
	OpenSMOKE::OpenSMOKEVectorDouble  y(thermodynamicsMapXML->NumberOfSpecies() + 1);
	OpenSMOKE::OpenSMOKEVectorDouble dy(thermodynamicsMapXML->NumberOfSpecies() + 1);
	for (unsigned int i = 0; i < thermodynamicsMapXML->NumberOfSpecies(); i++)
		y[i + 1] = initial_omega[i + 1];
	y[thermodynamicsMapXML->NumberOfSpecies() + 1] = epsilon0;
	reactor1d->Equations(0., y.GetHandle(), dy.GetHandle());
	*/


	/*


	porous_medium->SetTemperature(inlet_T);
	porous_medium->SetPressure(inlet_P);
	
	OpenSMOKE::OpenSMOKEVectorDouble inlet_x;
	double inlet_MW;
	thermodynamicsMapXML->MoleFractions_From_MassFractions(inlet_x, inlet_MW, inlet_omega);

	Eigen::VectorXd mole_fractions(thermodynamicsMapXML->NumberOfSpecies());
	for (unsigned int i = 0; i < mole_fractions.size(); i++)
		mole_fractions(i) = inlet_x[i + 1];

	porous_medium->EffectiveDiffusionCoefficients(mole_fractions);

	for (unsigned int i = 0; i < mole_fractions.size(); i++)
		std::cout << i << " " << thermodynamicsMapXML->NamesOfSpecies()[i] << " " <<
		porous_medium->gamma_effective()(i) << " " <<
		porous_medium->gamma_fick_effective()(i) << " " <<
		porous_medium->gamma_knudsen_effective()(i) << " " <<
		porous_medium->gamma_fick()(i) << " " <<
		porous_medium->gamma_knudsen()(i) << " " << std::endl;

	const double cTot = inlet_P / PhysicalConstants::R_J_kmol / inlet_T;
	Eigen::VectorXd c = mole_fractions; c *= cTot;

	porous_medium->FormationRates(c);

	for (unsigned int i = 0; i < 4; i++)
		std::cout << porous_medium->r()(i) << std::endl;

	std::cout << porous_medium->I_CH4() << std::endl;
	std::cout << porous_medium->I_C2H4() << std::endl;
	std::cout << porous_medium->I_C2H2() << std::endl;
	std::cout << porous_medium->I_C6H6() << std::endl;

	for (unsigned int i = 0; i < mole_fractions.size(); i++)
		std::cout << i << " " << thermodynamicsMapXML->NamesOfSpecies()[i] << " " << c(i) << " " << porous_medium->R()(i) << std::endl;
		*/
	return 0;
}