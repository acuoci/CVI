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

// OpenSMOKE++ Definitions
#include "OpenSMOKEpp"

// CHEMKIN maps
#include "maps/Maps_CHEMKIN"

// Reactor utilities
#include "idealreactors/utilities/Utilities"
#include "surfacereactors/utilities/Utilities"
#include "utilities/ropa/OnTheFlyROPA.h"

// 1D grid
#include "utilities/grids/adaptive/Grid1D.h"

// Grammars
#include "Grammar_CVI_Solver.h"
#include "Grammar_CVI_PlugFlowReactorCoupled.h"
#include "Grammar_CVI_PorousMedium.h"
#include "Grammar_CVI_HeterogeneousMechanism.h"
#include "Grammar_CVI_PorosityDefect.h"

// Heterogeneous mechanism
#include "HeterogeneousMechanism.h"
#include "HeterogeneousDetailedMechanism.h"

// Porous medium
#include "PorousMedium.h"
#include "PorosityDefect.h"

// Plug Flow Reactor
#include "PlugFlowReactorCoupled.h"
#include "PlugFlowReactorCoupledProfiles.h"

// Disk from CFD simulation
#include "DiskFromCFD.h"

// Capillary 1D
#include "Capillary.h"

// Reactor 1D
#include "Reactor1D.h"

// Reactor 2D
#include "Reactor2D.h"

// Numerical parameters
#include "math/native-dae-solvers/parameters/DaeSolver_Parameters.h"

#include "math/native-nls-solvers/NonLinearSolver.h"
#include "math/native-nls-solvers/KernelSparse.h"

int main(int argc, char** argv)
{
	boost::filesystem::path executable_file = OpenSMOKE::GetExecutableFileName(argv);
	boost::filesystem::path executable_folder = executable_file.parent_path();

	OpenSMOKE::OpenSMOKE_logo("CVI_Solver", "Alberto Cuoci (alberto.cuoci@polimi.it)");

	std::string input_file_name_ = "input.dic";
	std::string main_dictionary_name_ = "CVI_Solver";

	// Program options from command line
	{
		namespace po = boost::program_options;
		po::options_description description("Options for the CVI_Solver");
		description.add_options()
			("help", "print help messages")
			("input", po::value<std::string>(), "name of the file containing the main dictionary (default \"input.dic\")")
			("dictionary", po::value<std::string>(), "name of the main dictionary to be used (default \"CVI_Solver1D\")");

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

	// Gaseous phase
	CVI::GaseousPhase gaseous_phase = CVI::GASEOUS_PHASE_FROM_PLUG_FLOW;

	// Defines the grammar rules
	CVI::Grammar_CVI_Solver						grammar_cvi_solver;
	CVI::Grammar_CVI_PlugFlowReactorCoupled		grammar_cvi_plug_flow_reactor_coupled;
	CVI::Grammar_CVI_PorousMedium				grammar_cvi_porous_medium;
	CVI::Grammar_CVI_PorosityDefect				grammar_defect_porous_medium;
	CVI::Grammar_CVI_HeterogeneousMechanism		grammar_cvi_heterogeneous_mechanism;

	// Define the dictionaries
	OpenSMOKE::OpenSMOKE_DictionaryManager dictionaries;
	dictionaries.ReadDictionariesFromFile(input_file_name_);
	dictionaries(main_dictionary_name_).SetGrammar(grammar_cvi_solver);

	// Type of problem
	enum CVIProblemType { CVI_CAPILLARY, CVI_REACTOR1D, CVI_REACTOR2D } problem_type;
	{
		std::string value;
		if (dictionaries(main_dictionary_name_).CheckOption("@Type") == true)
			dictionaries(main_dictionary_name_).ReadString("@Type", value);
		if (value == "1D")					problem_type = CVI_REACTOR1D;
		else if (value == "2D")				problem_type = CVI_REACTOR2D;
		else if (value == "Capillary")		problem_type = CVI_CAPILLARY;
		else OpenSMOKE::FatalErrorMessage("Wrong @Type: Capillary | 1D | 2D");
	}

	// Type of problem
	bool symmetry_planar = true;
	{
		std::string value;
		if (dictionaries(main_dictionary_name_).CheckOption("@Symmetry") == true)
			dictionaries(main_dictionary_name_).ReadString("@Symmetry", value);
		if (value == "Planar")				symmetry_planar = true;
		else if (value == "Cylindrical")	symmetry_planar = false;
		else OpenSMOKE::FatalErrorMessage("Wrong @Symmetry: Planar | Cylindrical");
	}

	// Type of problem
	CVI::PorosityTreatment porosity_treatment = CVI::PorosityTreatment::POROSITY_COUPLED;
	{
		if (dictionaries(main_dictionary_name_).CheckOption("@PorosityTreatment") == true)
		{
			std::string value;
			dictionaries(main_dictionary_name_).ReadString("@PorosityTreatment", value);
			if (value == "coupled")						porosity_treatment = CVI::PorosityTreatment::POROSITY_COUPLED;
			else if (value == "decoupled-cumulative")	porosity_treatment = CVI::PorosityTreatment::POROSITY_DECOUPLED_CUMULATIVE;
			else if (value == "decoupled-finalvalue")	porosity_treatment = CVI::PorosityTreatment::POROSITY_DECOUPLED_FINALVALUE;
			else OpenSMOKE::FatalErrorMessage("Wrong @PorosityTreatment: coupled | decoupled-cumulative | decoupled-finalvalue");
		}
	}

	// Plug flow reactor
	std::string dict_name_plug_flow;
	if (dictionaries(main_dictionary_name_).CheckOption("@PlugFlowReactor") == true)
	{
		dictionaries(main_dictionary_name_).ReadDictionary("@PlugFlowReactor", dict_name_plug_flow);
		dictionaries(dict_name_plug_flow).SetGrammar(grammar_cvi_plug_flow_reactor_coupled);
		gaseous_phase = CVI::GASEOUS_PHASE_FROM_PLUG_FLOW;
	}
	else
	{
		if (problem_type == CVI_REACTOR1D || problem_type == CVI_CAPILLARY)
			OpenSMOKE::FatalErrorMessage("@PlugFlowReactor is compulsory for 1D and Capillary problems");
	}

	// Read from Disk file
	boost::filesystem::path disk_file_name;
	if (dictionaries(main_dictionary_name_).CheckOption("@DiskFromCFD") == true)
	{
		dictionaries(main_dictionary_name_).ReadPath("@DiskFromCFD", disk_file_name);
		gaseous_phase = CVI::GASEOUS_PHASE_FROM_CFD;
	}

	// Profile
	bool is_temperature_profile = false;
	OpenSMOKE::FixedProfile* temperature_profile;
	{
		std::string name_of_gas_status_subdictionary;
		if (dictionaries(main_dictionary_name_).CheckOption("@TemperatureProfile") == true)
		{
			if (problem_type == CVI_REACTOR1D && gaseous_phase != CVI::GASEOUS_PHASE_FROM_PLUG_FLOW)
				OpenSMOKE::FatalErrorMessage("The @TemperatureProfile cannot be used without the @PlugFlowReactor option");

			if (problem_type == CVI_REACTOR2D && gaseous_phase != CVI::GASEOUS_PHASE_FROM_CFD)
				OpenSMOKE::FatalErrorMessage("The @TemperatureProfile cannot be used without the @DiskFromCFD option");

			if (problem_type == CVI_CAPILLARY)
				OpenSMOKE::FatalErrorMessage("The @TemperatureProfile cannot be used with the Capillary geometry");

			dictionaries(main_dictionary_name_).ReadDictionary("@TemperatureProfile", name_of_gas_status_subdictionary);

			OpenSMOKE::OpenSMOKEVectorDouble x, y;
			std::string x_variable, y_variable;
			GetXYProfileFromDictionary(dictionaries(name_of_gas_status_subdictionary), x, y, x_variable, y_variable);

			if (x_variable != "time")
				OpenSMOKE::FatalErrorMessage("The @TemperatureProfile must be defined versus the time");
			
			is_temperature_profile = true;
			temperature_profile = new OpenSMOKE::FixedProfile(x.Size(), x.GetHandle(), y.GetHandle());
		}
	}

	// Profile
	std::vector<std::string> impermeable_walls;
	if (dictionaries(main_dictionary_name_).CheckOption("@ImpermeableWalls") == true)
		dictionaries(main_dictionary_name_).ReadOption("@ImpermeableWalls", impermeable_walls);

	// Porous medium
	std::string dict_name_porous_medium;
	if (dictionaries(main_dictionary_name_).CheckOption("@PorousMedium") == true)
		dictionaries(main_dictionary_name_).ReadDictionary("@PorousMedium", dict_name_porous_medium);

	// Porous medium
	std::string dict_name_heterogeneous_mechanism;
	if (dictionaries(main_dictionary_name_).CheckOption("@HeterogeneousMechanism") == true)
		dictionaries(main_dictionary_name_).ReadDictionary("@HeterogeneousMechanism", dict_name_heterogeneous_mechanism);

	// Porosity defect
	CVI::PorosityDefect* porosity_defect = new CVI::PorosityDefect();
	if (dictionaries(main_dictionary_name_).CheckOption("@PorosityDefect") == true)
	{
		std::string dict_name_porosity_defect;
		dictionaries(main_dictionary_name_).ReadDictionary("@PorosityDefect", dict_name_porosity_defect);
		dictionaries(dict_name_porosity_defect).SetGrammar(grammar_defect_porous_medium);
		porosity_defect->ReadFromDictionary(dictionaries(dict_name_porosity_defect));
	}

	// Sets the grammars
	dictionaries(dict_name_porous_medium).SetGrammar(grammar_cvi_porous_medium);
	dictionaries(dict_name_heterogeneous_mechanism).SetGrammar(grammar_cvi_heterogeneous_mechanism);

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

		boost::filesystem::path path_input_surface_kinetics;
		if (dictionaries(name_of_rapid_kinetics_subdictionary).CheckOption("@Surface") == true)
			dictionaries(name_of_rapid_kinetics_subdictionary).ReadPath("@Surface", path_input_surface_kinetics);
		else OpenSMOKE::FatalErrorMessage("The @Surface option is strictly required by the @KineticsPreProcessor");

		if (dictionaries(name_of_rapid_kinetics_subdictionary).CheckOption("@Output") == true)
			dictionaries(name_of_rapid_kinetics_subdictionary).ReadPath("@Output", path_kinetics_output);

		OpenSMOKE::RapidKineticMechanismWithTransport(	path_kinetics_output, 
														path_input_transport.c_str(),
														path_input_thermodynamics.c_str(),
														path_input_kinetics.c_str(),
														path_input_surface_kinetics.c_str()	);
	}

	// Read thermodynamics and kinetics maps
	OpenSMOKE::ThermodynamicsMap_CHEMKIN*			thermodynamicsMapXML;
	OpenSMOKE::KineticsMap_CHEMKIN*					kineticsMapXML;
	OpenSMOKE::TransportPropertiesMap_CHEMKIN*		transportMapXML;
	OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN*	thermodynamicsSurfaceMapXML;
	OpenSMOKE::KineticsMap_Surface_CHEMKIN*			kineticsSurfaceMapXML;

	// Read the homogeneous kinetic scheme in XML format
	{
		// Read names of species
		boost::property_tree::ptree ptree;
		boost::property_tree::read_xml((path_kinetics_output / "kinetics.xml").string(), ptree);

		double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
		thermodynamicsMapXML = new OpenSMOKE::ThermodynamicsMap_CHEMKIN(ptree);
		kineticsMapXML = new OpenSMOKE::KineticsMap_CHEMKIN(*thermodynamicsMapXML, ptree);
		transportMapXML = new OpenSMOKE::TransportPropertiesMap_CHEMKIN(ptree);

		// Transport properties
		transportMapXML = new OpenSMOKE::TransportPropertiesMap_CHEMKIN(ptree);

		double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
		std::cout << "Time to read XML file: " << tEnd - tStart << std::endl;
	}

	// Read the surface kinetic scheme in XML format
	{
		// Read names of species
		boost::property_tree::ptree ptree;
		boost::property_tree::read_xml((path_kinetics_output / "kinetics.surface.xml").string(), ptree);

		double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
		thermodynamicsSurfaceMapXML = new OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN(ptree);
		kineticsSurfaceMapXML = new OpenSMOKE::KineticsMap_Surface_CHEMKIN(*thermodynamicsSurfaceMapXML, ptree);
		double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
		std::cout << "Time to read XML file: " << tEnd - tStart << std::endl;
	}

	bool detailed_heterogeneous_kinetics = true;;
	if (dictionaries(main_dictionary_name_).CheckOption("@DetailedSurfaceChemistry") == true)
	{
		dictionaries(main_dictionary_name_).ReadBool("@DetailedSurfaceChemistry", detailed_heterogeneous_kinetics);
	}

	std::string gas_dae_species = "none";
	if (dictionaries(main_dictionary_name_).CheckOption("@GasDaeSpecies") == true)
	{
		dictionaries(main_dictionary_name_).ReadString("@GasDaeSpecies", gas_dae_species);
	}

	std::string surface_dae_species = "none";
	if (dictionaries(main_dictionary_name_).CheckOption("@SurfaceDaeSpecies") == true)
	{
		dictionaries(main_dictionary_name_).ReadString("@SurfaceDaeSpecies", surface_dae_species);
	}

	boost::filesystem::path output_path;
	if (dictionaries(main_dictionary_name_).CheckOption("@Output") == true)
	{
		dictionaries(main_dictionary_name_).ReadPath("@Output", output_path);
		OpenSMOKE::CreateDirectory(output_path);
	}

	// On the fly ROPA
	OpenSMOKE::SurfaceOnTheFlyROPA* onTheFlyROPA;
	bool on_the_fly_ropa = false;
	if (detailed_heterogeneous_kinetics == true)
	{
		if (dictionaries(main_dictionary_name_).CheckOption("@OnTheFlyROPA") == true)
		{
			std::string name_of_options_subdictionary;
			dictionaries(main_dictionary_name_).ReadDictionary("@OnTheFlyROPA", name_of_options_subdictionary);

			on_the_fly_ropa = true;
			onTheFlyROPA = new OpenSMOKE::SurfaceOnTheFlyROPA(*thermodynamicsSurfaceMapXML, *kineticsSurfaceMapXML);
			onTheFlyROPA->SetupFromDictionary(dictionaries(name_of_options_subdictionary), path_kinetics_output);
		}
	}

	// Monodimensional grid along the x axis
	OpenSMOKE::Grid1D* grid_x;
	{
		int number_of_points;
		if (dictionaries(main_dictionary_name_).CheckOption("@XPoints") == true)
			dictionaries(main_dictionary_name_).ReadInt("@XPoints", number_of_points);

		double stretching_factor = 1.;
		if (dictionaries(main_dictionary_name_).CheckOption("@XStretchingFactor") == true)
			dictionaries(main_dictionary_name_).ReadDouble("@XStretchingFactor", stretching_factor);

		if (symmetry_planar == true)
		{
			double length;
			std::string units;
			dictionaries(main_dictionary_name_).ReadMeasure("@XLength", length, units);
			if (units == "m")			length = length;
			else if (units == "cm")		length = length / 1.e2;
			else if (units == "mm")		length = length / 1.e3;
			else if (units == "micron")	length = length / 1.e6;
			else OpenSMOKE::FatalErrorMessage("Unknown length units");

			const double beta = std::pow(stretching_factor, 1. / static_cast<double>(number_of_points - 2));
			grid_x = new OpenSMOKE::Grid1D(number_of_points, 0., length, 1./beta);
		}
		else
		{
			std::string units;

			double radius_internal;
			dictionaries(main_dictionary_name_).ReadMeasure("@RadiusInternal", radius_internal, units);
			if (units == "m")			radius_internal = radius_internal;
			else if (units == "cm")		radius_internal = radius_internal / 1.e2;
			else if (units == "mm")		radius_internal = radius_internal / 1.e3;
			else if (units == "micron")	radius_internal = radius_internal / 1.e6;
			else OpenSMOKE::FatalErrorMessage("Unknown length units");

			double radius_external;
			dictionaries(main_dictionary_name_).ReadMeasure("@RadiusExternal", radius_external, units);
			if (units == "m")			radius_external = radius_external;
			else if (units == "cm")		radius_external = radius_external / 1.e2;
			else if (units == "mm")		radius_external = radius_external / 1.e3;
			else if (units == "micron")	radius_external = radius_external / 1.e6;
			else OpenSMOKE::FatalErrorMessage("Unknown length units");

			if (radius_internal == 0.)
			{
				const double beta = std::pow(stretching_factor, 1. / static_cast<double>(number_of_points - 2));
				grid_x = new OpenSMOKE::Grid1D(number_of_points, radius_internal, radius_external, 1./beta);
			}
			else
			{
				if (stretching_factor != 1.)
				{
					if (number_of_points % 2 != 0)
						OpenSMOKE::FatalErrorMessage("In case of stretcheg grid along the x axis, the number of points must be odd.");

					const double beta = std::pow(stretching_factor, 1. / static_cast<double>((number_of_points - 1)/2-1));
					OpenSMOKE::Grid1D grid_left = OpenSMOKE::Grid1D((number_of_points+1)/2, radius_internal, (radius_external- radius_internal)/2., beta);
					OpenSMOKE::Grid1D grid_right = OpenSMOKE::Grid1D((number_of_points + 1) / 2, (radius_external - radius_internal) / 2., radius_external, 1./beta);
					
					std::vector<double> x;
					for (int i = 0; i < (number_of_points + 1) / 2; i++)
						x.push_back(grid_left.x()[i]);
					for (int i = 1; i < (number_of_points + 1) / 2; i++)
						x.push_back(grid_right.x()[i]);

					grid_x = new OpenSMOKE::Grid1D(x);
				}
				else
				{
					grid_x = new OpenSMOKE::Grid1D(number_of_points, radius_internal, radius_external, 1.);
				}
			}
		}
	}

	// Monodimensional grid along the y axis
	OpenSMOKE::Grid1D* grid_y;
	{
		int number_of_points;
		if (dictionaries(main_dictionary_name_).CheckOption("@YPoints") == true)
			dictionaries(main_dictionary_name_).ReadInt("@YPoints", number_of_points);

		double stretching_factor = 1.;
		if (dictionaries(main_dictionary_name_).CheckOption("@YStretchingFactor") == true)
			dictionaries(main_dictionary_name_).ReadDouble("@YStretchingFactor", stretching_factor);

		double length;
		std::string units;
		if (dictionaries(main_dictionary_name_).CheckOption("@YLength") == true)
		{
			dictionaries(main_dictionary_name_).ReadMeasure("@YLength", length, units);
			if (units == "m")			length = length;
			else if (units == "cm")		length = length / 1.e2;
			else if (units == "mm")		length = length / 1.e3;
			else if (units == "micron")	length = length / 1.e6;
			else OpenSMOKE::FatalErrorMessage("Unknown length units");
		}

		if (stretching_factor != 1.)
		{
			if (number_of_points % 2 == 0)
				OpenSMOKE::FatalErrorMessage("In case of stretcheg grid along the y axis, the number of points must be odd.");

			const double beta = std::pow(stretching_factor, 1. / static_cast<double>((number_of_points - 1) / 2 - 1));
			OpenSMOKE::Grid1D grid_left = OpenSMOKE::Grid1D((number_of_points + 1) / 2, 0, length / 2., beta);
			OpenSMOKE::Grid1D grid_right = OpenSMOKE::Grid1D((number_of_points + 1) / 2, length / 2., length, 1. / beta);

			std::vector<double> y;
			for (int i = 0; i < (number_of_points + 1) / 2; i++)
				y.push_back(grid_left.x()[i]);
			for (int i = 1; i < (number_of_points + 1) / 2; i++)
				y.push_back(grid_right.x()[i]);

			grid_y = new OpenSMOKE::Grid1D(y);
		}
		else
		{
			grid_y = new OpenSMOKE::Grid1D(number_of_points, 0., length, 1.);
		}
	}

	// Dae Options
	DaeSMOKE::DaeSolver_Parameters* dae_parameters;
	dae_parameters = new DaeSMOKE::DaeSolver_Parameters();
	if (dictionaries(main_dictionary_name_).CheckOption("@DaeParameters") == true)
	{
		std::string name_of_subdictionary;
		dictionaries(main_dictionary_name_).ReadDictionary("@DaeParameters", name_of_subdictionary);
		dae_parameters->SetupFromDictionary(dictionaries(name_of_subdictionary));
	}
	dae_parameters->SetMinimumMeanThreshold(0.);

	// Ode Options
	OdeSMOKE::OdeSolver_Parameters* ode_parameters;
	ode_parameters = new OdeSMOKE::OdeSolver_Parameters();
	if (dictionaries(main_dictionary_name_).CheckOption("@OdeParameters") == true)
	{
		std::string name_of_subdictionary;
		dictionaries(main_dictionary_name_).ReadDictionary("@OdeParameters", name_of_subdictionary);
		ode_parameters->SetupFromDictionary(dictionaries(name_of_subdictionary));
	}

	// Read the plug flow residence time 
	double residence_time = 2.;
	{
		std::string units;

		if (dictionaries(main_dictionary_name_).CheckOption("@ResidenceTime") == true)
		{
			dictionaries(main_dictionary_name_).ReadMeasure("@ResidenceTime", residence_time, units);
			if (units == "s")				residence_time = residence_time;
			else if (units == "ms")			residence_time = residence_time / 1.e3;
			else if (units == "min")		residence_time = residence_time * 60.;
			else OpenSMOKE::FatalErrorMessage("Unknown @ResidenceTime units");
		}
		else
		{
			if (problem_type == CVI_REACTOR1D || problem_type == CVI_CAPILLARY)
				OpenSMOKE::FatalErrorMessage("@ResidenceTime must be specified");
		}
	}

	// Uniform velocity inside the preform
	double vx = 0.;
	double vy = 0.;
	{
		std::string units;

		if (dictionaries(main_dictionary_name_).CheckOption("@XVelocity") == true)
		{
			dictionaries(main_dictionary_name_).ReadMeasure("@XVelocity", vx, units);
			if (units == "m/s")				vx = vx;
			else if (units == "cm/s")		vx /= 1.e2;
			else if (units == "mm/s")		vx /= 1.e3;
			else OpenSMOKE::FatalErrorMessage("Unknown @XVelocity units");
		}

		if (dictionaries(main_dictionary_name_).CheckOption("@YVelocity") == true)
		{
			dictionaries(main_dictionary_name_).ReadMeasure("@YVelocity", vy, units);
			if (units == "m/s")				vy = vy;
			else if (units == "cm/s")		vy /= 1.e2;
			else if (units == "mm/s")		vy /= 1.e3;
			else OpenSMOKE::FatalErrorMessage("Unknown @YVelocity units");
		}
	}

	// Read the capillary diameter
	double capillary_diameter = 1.e-3;
	{
		std::string units;

		if (dictionaries(main_dictionary_name_).CheckOption("@CapillaryDiameter") == true)
		{
			dictionaries(main_dictionary_name_).ReadMeasure("@CapillaryDiameter", capillary_diameter, units);
			if (units == "m")			capillary_diameter = capillary_diameter;
			else if (units == "cm")			capillary_diameter = capillary_diameter / 1.e2;
			else if (units == "mm")			capillary_diameter = capillary_diameter / 1.e3;
			else if (units == "micron")		capillary_diameter = capillary_diameter / 1.e6;
			else OpenSMOKE::FatalErrorMessage("Unknown @CapillaryDiameter units");
		}
		else
		{
			if (problem_type == CVI_CAPILLARY)
				OpenSMOKE::FatalErrorMessage("@CapillaryDiameter must be specified");
		}
	}

	// Read inlet conditions
	double inlet_T;
	double inlet_P;
	Eigen::VectorXd inlet_omega(thermodynamicsMapXML->NumberOfSpecies());
	{
		std::string dict_name;
		if (dictionaries(main_dictionary_name_).CheckOption("@InletStream") == true)
			dictionaries(main_dictionary_name_).ReadDictionary("@InletStream", dict_name);

		OpenSMOKE::OpenSMOKEVectorDouble aux(thermodynamicsMapXML->NumberOfSpecies());
		GetGasStatusFromDictionary(dictionaries(dict_name), *thermodynamicsMapXML, inlet_T, inlet_P, aux);
		for (int i = 1; i <= aux.Size(); i++)
			inlet_omega(i - 1) = aux[i];
	}

	// Read initial conditions (uniform)
	double initial_T;
	double initial_P;
	Eigen::VectorXd initial_omega(thermodynamicsMapXML->NumberOfSpecies());
	{
		std::string dict_name;
		if (dictionaries(main_dictionary_name_).CheckOption("@InitialConditions") == true)
			dictionaries(main_dictionary_name_).ReadDictionary("@InitialConditions", dict_name);

		OpenSMOKE::OpenSMOKEVectorDouble aux(thermodynamicsMapXML->NumberOfSpecies());
		GetGasStatusFromDictionary(dictionaries(dict_name), *thermodynamicsMapXML, initial_T, initial_P, aux);
		for (int i = 1; i <= aux.Size(); i++)
			initial_omega(i - 1) = aux[i];
	}

	// Read initial conditions
	OpenSMOKE::OpenSMOKEVectorDouble Z0(thermodynamicsSurfaceMapXML->number_of_site_species());
	Eigen::VectorXd Gamma0(thermodynamicsSurfaceMapXML->number_of_site_phases(0));
	std::vector<bool> SiteNonConservation(thermodynamicsSurfaceMapXML->number_of_site_phases(0));
	for (unsigned int k = 0; k<thermodynamicsSurfaceMapXML->number_of_site_phases(0); k++)
	{
		bool site_non_conservation;
		const std::string name = "Surface-" + thermodynamicsSurfaceMapXML->matrix_names_site_phases()[0][k];
		GetSurfaceCompositionFromDictionary
		(dictionaries(name), *thermodynamicsSurfaceMapXML, Z0, site_non_conservation);

		SiteNonConservation[k] = site_non_conservation;
		Gamma0(k) = thermodynamicsSurfaceMapXML->matrix_densities_site_phases()[0][k];
	}

	// Interval time
	double dae_time_interval = 0.;
	{
		std::string units;
		if (dictionaries(main_dictionary_name_).CheckOption("@DaeTimeInterval") == true)
		{
			dictionaries(main_dictionary_name_).ReadMeasure("@DaeTimeInterval", dae_time_interval, units);
			if (units == "s")				dae_time_interval = dae_time_interval;
			else if (units == "ms")			dae_time_interval = dae_time_interval / 1.e3;
			else if (units == "min")		dae_time_interval = dae_time_interval * 60.;
			else if (units == "h")			dae_time_interval = dae_time_interval * 3600.;
			else OpenSMOKE::FatalErrorMessage("Unknown @DaeTimeInterval units");
		}
	}

	// Interval time
	double ode_end_time = 1.;
	{
		std::string units;
		if (dictionaries(main_dictionary_name_).CheckOption("@OdeEndTime") == true)
		{
			dictionaries(main_dictionary_name_).ReadMeasure("@OdeEndTime", ode_end_time, units);
			if (units == "s")				ode_end_time = ode_end_time;
			else if (units == "ms")			ode_end_time = ode_end_time / 1.e3;
			else if (units == "min")		ode_end_time = ode_end_time * 60.;
			else if (units == "h")			ode_end_time = ode_end_time * 3600.;
			else OpenSMOKE::FatalErrorMessage("Unknown @OdeEndTime units");
		}
	}

	// Interval time
	double tecplot_time_interval = 0.;
	{
		std::string units;
		if (dictionaries(main_dictionary_name_).CheckOption("@TecplotTimeInterval") == true)
		{
			dictionaries(main_dictionary_name_).ReadMeasure("@TecplotTimeInterval", tecplot_time_interval, units);
			if (units == "s")				tecplot_time_interval = tecplot_time_interval;
			else if (units == "ms")			tecplot_time_interval = tecplot_time_interval / 1.e3;
			else if (units == "min")		tecplot_time_interval = tecplot_time_interval * 60.;
			else if (units == "h")			tecplot_time_interval = tecplot_time_interval * 3600.;
			else OpenSMOKE::FatalErrorMessage("Unknown @TecplotTimeInterval units");
		}
	}

	// Total time of simulation
	double time_total = 0.;
	{
		std::string units;
		if (dictionaries(main_dictionary_name_).CheckOption("@TimeTotal") == true)
		{
			dictionaries(main_dictionary_name_).ReadMeasure("@TimeTotal", time_total, units);
			if (units == "s")				time_total = time_total;
			else if (units == "ms")			time_total = time_total / 1.e3;
			else if (units == "min")		time_total = time_total * 60.;
			else if (units == "h")			time_total = time_total * 3600.;
			else OpenSMOKE::FatalErrorMessage("Unknown @TimeTotal units");
		}
	}

	// Steps video
	int steps_video = 0;
	{
		if (dictionaries(main_dictionary_name_).CheckOption("@StepsVideo") == true)
			dictionaries(main_dictionary_name_).ReadInt("@StepsVideo", steps_video);
	}

	// Steps file
	int steps_file = 0;
	{
		if (dictionaries(main_dictionary_name_).CheckOption("@StepsFile") == true)
			dictionaries(main_dictionary_name_).ReadInt("@StepsFile", steps_file);
	}

	// Steps update plug flow
	int steps_update_plug_flow = 0;
	{
		if (dictionaries(main_dictionary_name_).CheckOption("@StepsPlugFlow") == true)
			dictionaries(main_dictionary_name_).ReadInt("@StepsPlugFlow", steps_update_plug_flow);
	}


	// Derivatives
	OpenSMOKE::derivative_type derivative_mass_fractions = OpenSMOKE::DERIVATIVE_1ST_CENTERED;
	OpenSMOKE::derivative_type derivative_effective_diffusivity = OpenSMOKE::DERIVATIVE_1ST_CENTERED;
	OpenSMOKE::derivative_type derivative_bulk_density = OpenSMOKE::DERIVATIVE_1ST_CENTERED;
	{
		std::string dummy;
		
		if (dictionaries(main_dictionary_name_).CheckOption("@DerivativeMassFractions") == true)
		{
			dictionaries(main_dictionary_name_).ReadString("@DerivativeMassFractions", dummy);
			if (dummy == "Backward")
				derivative_mass_fractions = OpenSMOKE::DERIVATIVE_1ST_BACKWARD;
			else if (dummy == "Forward")
				derivative_mass_fractions = OpenSMOKE::DERIVATIVE_1ST_FORWARD;
			else if (dummy == "Centered")
				derivative_mass_fractions = OpenSMOKE::DERIVATIVE_1ST_CENTERED;
		}

		if (dictionaries(main_dictionary_name_).CheckOption("@DerivativeBulkDensity") == true)
		{
			dictionaries(main_dictionary_name_).ReadString("@DerivativeBulkDensity", dummy);
			if (dummy == "Backward")
				derivative_bulk_density = OpenSMOKE::DERIVATIVE_1ST_BACKWARD;
			else if (dummy == "Forward")
				derivative_bulk_density = OpenSMOKE::DERIVATIVE_1ST_FORWARD;
			else if (dummy == "Centered")
				derivative_bulk_density = OpenSMOKE::DERIVATIVE_1ST_CENTERED;
		}

		if (dictionaries(main_dictionary_name_).CheckOption("@DerivativeEffectiveDiffusivity") == true)
		{
			dictionaries(main_dictionary_name_).ReadString("@DerivativeEffectiveDiffusivity", dummy);
			if (dummy == "Backward")
				derivative_effective_diffusivity = OpenSMOKE::DERIVATIVE_1ST_BACKWARD;
			else if (dummy == "Forward")
				derivative_effective_diffusivity = OpenSMOKE::DERIVATIVE_1ST_FORWARD;
			else if (dummy == "Centered")
				derivative_effective_diffusivity = OpenSMOKE::DERIVATIVE_1ST_CENTERED;
		}
	}

	// Read from backupfile
	boost::filesystem::path path_backup;
	bool readjust_backup = false;
	if (dictionaries(main_dictionary_name_).CheckOption("@Backup") == true)
	{
		dictionaries(main_dictionary_name_).ReadPath("@Backup", path_backup);
		if (dictionaries(main_dictionary_name_).CheckOption("@ReadjustBackup") == true)
			dictionaries(main_dictionary_name_).ReadBool("@ReadjustBackup", readjust_backup);
	}

	// Solve the 2D problem
	if (problem_type == CVI_REACTOR2D && gaseous_phase == CVI::GASEOUS_PHASE_FROM_PLUG_FLOW)
	{
		// Plug flow ractor simulation
		std::vector<Eigen::VectorXd> Y_gas_side;
		CVI::PlugFlowReactorCoupled* plug_flow_reactor = new CVI::PlugFlowReactorCoupled(*thermodynamicsMapXML, *kineticsMapXML, dictionaries(dict_name_plug_flow));
		
		// Initial plug flow reactor
		{
			// Set initial conditions
			plug_flow_reactor->SetInitialConditions(inlet_T, inlet_P, inlet_omega);

			if (plug_flow_reactor->coupling() == true)
			{
				if (plug_flow_reactor->geometric_pattern() != CVI::PlugFlowReactorCoupled::ONE_SIDE)
					OpenSMOKE::FatalErrorMessage("The coupling is currently available only for one-side geometric configuration");

				// The initial composition in the felt is assumed to be uniform
				unsigned int np = 3;
				Eigen::VectorXd csi_external(np);
				csi_external(0)=0.;
				csi_external(1)=1.e3;
				csi_external(2)=2.e3;

				Eigen::MatrixXd omega_external(np, thermodynamicsMapXML->NumberOfSpecies());
				omega_external.setConstant(0.);
				for(unsigned int i=0;i<np;i++)
					for(unsigned int j=0;j<thermodynamicsMapXML->NumberOfSpecies();j++)
						omega_external(i,j) = initial_omega(j);
				
				plug_flow_reactor->SetExternalMassFractionsProfile(csi_external, omega_external);
			}

			// Solve the plug flow reactor
			plug_flow_reactor->Solve(residence_time);

			// Extract the plug flow history
			Eigen::VectorXd csi(plug_flow_reactor->history_csi().size());
			for (unsigned int i = 0; i < plug_flow_reactor->history_csi().size(); i++)
				csi(i) = plug_flow_reactor->history_csi()[i];

			// Assign the boundary conditions
			PlugFlowReactorCoupledProfiles* profiles = new PlugFlowReactorCoupledProfiles(csi);
			Y_gas_side.resize(2 * grid_x->Np() + grid_y->Np());
			for (unsigned int i = 0; i < Y_gas_side.size(); i++)
			{
				Y_gas_side[i].resize(thermodynamicsMapXML->NumberOfSpecies());
				Y_gas_side[i].setZero();
			}

			// Case: only east side
			if (plug_flow_reactor->geometric_pattern() == CVI::PlugFlowReactorCoupled::ONE_SIDE)
			{
				for (int i = 0; i < grid_y->Np(); i++)
				{
					const int point = i + grid_x->Np();
					const double coordinate = grid_y->x()(i)+plug_flow_reactor->inert_length();
					profiles->Interpolate(coordinate, plug_flow_reactor->history_Y(), Y_gas_side[point]);
				}
			}
			// Case: south, east, and north sides
			else if (plug_flow_reactor->geometric_pattern() == CVI::PlugFlowReactorCoupled::THREE_SIDES)
			{
				for (int i = 0; i < grid_x->Np(); i++)
				{
					const double coordinate = plug_flow_reactor->inert_length() + grid_x->x()(i);
					profiles->Interpolate(coordinate, plug_flow_reactor->history_Y(), Y_gas_side[i]);
				}

				for (int i = 0; i < grid_y->Np(); i++)
				{
					const int point = i + grid_x->Np();
					const double coordinate = plug_flow_reactor->inert_length() + grid_x->L() + grid_y->x()(i);
					profiles->Interpolate(coordinate, plug_flow_reactor->history_Y(), Y_gas_side[point]);
				}

				for (int i = 0; i < grid_x->Np(); i++)
				{
					const int point = 2 * grid_x->Np() - 1 + grid_y->Np() - i;
					const double coordinate = plug_flow_reactor->inert_length() + grid_x->L() + grid_y->L() + grid_x->x()(grid_x->Np() - 1 - i);
					profiles->Interpolate(coordinate, plug_flow_reactor->history_Y(), Y_gas_side[point]);
				}
			}

			// Write plug-flow reactor profile
			{
				const boost::filesystem::path plug_flow_file_path = output_path / "plugflow.out";
				plug_flow_reactor->Print(0.,plug_flow_file_path);
			}
		}

		// Set heterogeneous mechanism
		CVI::HeterogeneousMechanism* heterogeneous_mechanism = new CVI::HeterogeneousMechanism(*thermodynamicsMapXML, *kineticsMapXML, *transportMapXML, dictionaries(dict_name_heterogeneous_mechanism));

		// Set detailed heterogeneous mechanism
		CVI::HeterogeneousDetailedMechanism* heterogeneous_detailed_mechanism = new CVI::HeterogeneousDetailedMechanism(*thermodynamicsMapXML, *kineticsMapXML, *transportMapXML, *thermodynamicsSurfaceMapXML, *kineticsSurfaceMapXML, true, true);

		// Set porous medium
		CVI::PorousMedium* porous_medium = new CVI::PorousMedium(*thermodynamicsMapXML, *kineticsMapXML, *transportMapXML, dictionaries(dict_name_porous_medium));

		// Creates the reactor
		reactor2d = new CVI::Reactor2D(	*thermodynamicsMapXML, *kineticsMapXML, *transportMapXML,
										*thermodynamicsSurfaceMapXML, *kineticsSurfaceMapXML,
										*porous_medium, *porosity_defect,
										*heterogeneous_mechanism, *heterogeneous_detailed_mechanism,
										*grid_x, *grid_y, *plug_flow_reactor,
										detailed_heterogeneous_kinetics, SiteNonConservation, gas_dae_species, surface_dae_species, output_path, porosity_treatment);

		// Initial surface fractions
		Eigen::VectorXd initial_Z(thermodynamicsSurfaceMapXML->number_of_site_species());
		for (unsigned int i = 0; i<thermodynamicsSurfaceMapXML->number_of_site_species(); i++)
			initial_Z(i) = Z0[i + 1];

		// Set options
		reactor2d->SetPlanarSymmetry(symmetry_planar);
		reactor2d->SetSiteNonConservation(SiteNonConservation);
		reactor2d->SetReadjustBackup(readjust_backup);
		reactor2d->SetImpermeableWalls(impermeable_walls);
		reactor2d->SetInitialConditions(path_backup, initial_T, initial_P, initial_omega, Gamma0, initial_Z);
		reactor2d->SetGasSide(inlet_T, inlet_P, Y_gas_side);
		reactor2d->SetUniformVelocity(vx, vy);
		reactor2d->SetTimeTotal(time_total);
		reactor2d->SetDaeTimeInterval(dae_time_interval);
		reactor2d->SetOdeEndTime(ode_end_time);
		reactor2d->SetTecplotTimeInterval(tecplot_time_interval);
		reactor2d->SetDerivativeMassFractions(derivative_mass_fractions);
		reactor2d->SetDerivativeEffectiveDiffusivity(derivative_effective_diffusivity);
		reactor2d->SetDerivativeBulkDensity(derivative_bulk_density);
		if (steps_video>0)				reactor2d->SetStepsVideo(steps_video);
		if (steps_file>0)				reactor2d->SetStepsFile(steps_file);
		if (steps_update_plug_flow>0)	reactor2d->SetStepsUpdatePlugFlow(steps_update_plug_flow);
		if (on_the_fly_ropa == true)	reactor2d->SetSurfaceOnTheFlyROPA(onTheFlyROPA);

		// Solve
		{
			time_t timerStart;
			time_t timerEnd;

			time(&timerStart);
			int flag = reactor2d->SolveFromScratch(*dae_parameters, *ode_parameters);
			time(&timerEnd);

			std::cout << "Total time: " << difftime(timerEnd, timerStart) << " s" << std::endl;
		}
	}	
	
	if (problem_type == CVI_REACTOR2D && gaseous_phase == CVI::GASEOUS_PHASE_FROM_CFD)
	{
		CVI::DiskFromCFD disk(*thermodynamicsMapXML, *kineticsMapXML, *grid_x, *grid_y);
		disk.ReadFromFile(disk_file_name, initial_P);
		disk.WriteOnFile(output_path);

		// Set heterogeneous mechanism
		CVI::HeterogeneousMechanism* heterogeneous_mechanism = new CVI::HeterogeneousMechanism(*thermodynamicsMapXML, *kineticsMapXML, *transportMapXML, dictionaries(dict_name_heterogeneous_mechanism));

		// Set detailed heterogeneous mechanism
		CVI::HeterogeneousDetailedMechanism* heterogeneous_detailed_mechanism = new CVI::HeterogeneousDetailedMechanism(*thermodynamicsMapXML, *kineticsMapXML, *transportMapXML, *thermodynamicsSurfaceMapXML, *kineticsSurfaceMapXML, true, true);

		// Set porous medium
		CVI::PorousMedium* porous_medium = new CVI::PorousMedium(*thermodynamicsMapXML, *kineticsMapXML, *transportMapXML, dictionaries(dict_name_porous_medium));

		// Creates the reactor
		CVI::PlugFlowReactorCoupled* plug_flow_reactor = nullptr; // dummy
														
		// Creates the reactor
		reactor2d = new CVI::Reactor2D(	*thermodynamicsMapXML, *kineticsMapXML, *transportMapXML,
										*thermodynamicsSurfaceMapXML, *kineticsSurfaceMapXML,
										*porous_medium, *porosity_defect,
										*heterogeneous_mechanism, *heterogeneous_detailed_mechanism,
										*grid_x, *grid_y, *plug_flow_reactor,
										detailed_heterogeneous_kinetics, SiteNonConservation, gas_dae_species, surface_dae_species, output_path, porosity_treatment);

		// Initial surface fractions
		Eigen::VectorXd initial_Z(thermodynamicsSurfaceMapXML->number_of_site_species());
		for (unsigned int i = 0; i<thermodynamicsSurfaceMapXML->number_of_site_species(); i++)
			initial_Z(i) = Z0[i + 1];

		// Set options
		reactor2d->SetPlanarSymmetry(symmetry_planar);
		reactor2d->SetSiteNonConservation(SiteNonConservation);
		reactor2d->SetReadjustBackup(readjust_backup);
		reactor2d->SetImpermeableWalls(impermeable_walls);
		reactor2d->SetInitialConditions(path_backup, initial_T, initial_P, initial_omega, Gamma0, initial_Z);
		
		// Set initial/boundary conditions
		if (is_temperature_profile == false)	reactor2d->SetGasSide(disk);
		else                                    reactor2d->SetGasSide(temperature_profile, disk);

		reactor2d->SetOutputFile(disk_file_name.stem().string());
		reactor2d->SetUniformVelocity(vx, vy);
		reactor2d->SetTimeTotal(time_total);
		reactor2d->SetDaeTimeInterval(dae_time_interval);
		reactor2d->SetOdeEndTime(ode_end_time);
		reactor2d->SetTecplotTimeInterval(tecplot_time_interval);
		reactor2d->SetDerivativeMassFractions(derivative_mass_fractions);
		reactor2d->SetDerivativeEffectiveDiffusivity(derivative_effective_diffusivity);
		reactor2d->SetDerivativeBulkDensity(derivative_bulk_density);
		if (steps_video>0)				reactor2d->SetStepsVideo(steps_video);
		if (steps_file>0)				reactor2d->SetStepsFile(steps_file);
		if (on_the_fly_ropa == true)	reactor2d->SetSurfaceOnTheFlyROPA(onTheFlyROPA);
		
		// Solve
		{
			time_t timerStart;
			time_t timerEnd;

			time(&timerStart);
			int flag = reactor2d->SolveFromScratch(*dae_parameters, *ode_parameters);
			time(&timerEnd);

			std::cout << "Total time: " << difftime(timerEnd, timerStart) << " s" << std::endl;
		}
	}
	
	// Solve the 1D problem
	if (problem_type == CVI_REACTOR1D && is_temperature_profile == false)
	{
		// Plug flow ractor simulation
		CVI::PlugFlowReactorCoupled* plug_flow_reactor = new CVI::PlugFlowReactorCoupled(*thermodynamicsMapXML, *kineticsMapXML, dictionaries(dict_name_plug_flow));
		{
			// Set initial conditions
			plug_flow_reactor->SetInitialConditions(inlet_T, inlet_P, inlet_omega);

			// Solve the plug flow reactor
			plug_flow_reactor->Solve(residence_time);
		}

		// Set heterogeneous mechanism
		CVI::HeterogeneousMechanism* heterogeneous_mechanism = new CVI::HeterogeneousMechanism(*thermodynamicsMapXML, *kineticsMapXML, *transportMapXML, dictionaries(dict_name_heterogeneous_mechanism));

		// Set heterogeneous mechanism
		CVI::HeterogeneousDetailedMechanism* heterogeneous_detailed_mechanism = new CVI::HeterogeneousDetailedMechanism(*thermodynamicsMapXML, *kineticsMapXML, *transportMapXML, *thermodynamicsSurfaceMapXML, *kineticsSurfaceMapXML, true, true);

		// Set porous medium
		CVI::PorousMedium* porous_medium = new CVI::PorousMedium(*thermodynamicsMapXML, *kineticsMapXML, *transportMapXML, dictionaries(dict_name_porous_medium));

		// Creates the reactor
		reactor1d = new CVI::Reactor1D(	*thermodynamicsMapXML, *kineticsMapXML, *transportMapXML, 
										*thermodynamicsSurfaceMapXML, *kineticsSurfaceMapXML, 
										*porous_medium, 
										*heterogeneous_mechanism, *heterogeneous_detailed_mechanism, 
										*grid_x, detailed_heterogeneous_kinetics,
										SiteNonConservation, surface_dae_species, output_path);

		// Initial surface fractions
		Eigen::VectorXd initial_Z(thermodynamicsSurfaceMapXML->number_of_site_species());
		for (unsigned int i = 0; i<thermodynamicsSurfaceMapXML->number_of_site_species(); i++)
			initial_Z(i) = Z0[i + 1];

		reactor1d->SetPlanarSymmetry(symmetry_planar);
		reactor1d->SetInitialConditions(initial_T, initial_P, initial_omega, Gamma0, initial_Z);
		reactor1d->SetSiteNonConservation(SiteNonConservation);
		reactor1d->SetGasSide(inlet_T, inlet_P, plug_flow_reactor->Y());
		reactor1d->SetTimeTotal(time_total);
		reactor1d->SetDaeTimeInterval(dae_time_interval);
		reactor1d->SetOdeEndTime(ode_end_time);
		reactor1d->SetDerivativeMassFractions(derivative_mass_fractions);
		reactor1d->SetDerivativeEffectiveDiffusivity(derivative_effective_diffusivity);
		reactor1d->SetDerivativeBulkDensity(derivative_bulk_density);

		if (on_the_fly_ropa == true)
			reactor1d->SetSurfaceOnTheFlyROPA(onTheFlyROPA);

		// Solve
		{
			time_t timerStart;
			time_t timerEnd;

			time(&timerStart);
			int flag = reactor1d->SolveFromScratch(*dae_parameters, *ode_parameters);
			time(&timerEnd);

			std::cout << "Total time: " << difftime(timerEnd, timerStart) << " s" << std::endl;
		}
	}

	// Solve the 1D problem
	if (problem_type == CVI_REACTOR1D && is_temperature_profile == true)
	{
		// Plug flow ractor simulation
		Eigen::MatrixXd omega_profiles_temp(temperature_profile->x().size(), thermodynamicsMapXML->NumberOfSpecies());
		for (unsigned int j = 0; j < temperature_profile->x().size(); j++)
		{
			std::cout << "Solving plug flow @T=" << temperature_profile->y()(j) << std::endl;
			CVI::PlugFlowReactorCoupled* plug_flow_reactor = new CVI::PlugFlowReactorCoupled(*thermodynamicsMapXML, *kineticsMapXML, dictionaries(dict_name_plug_flow));

			// Set initial conditions
			plug_flow_reactor->SetInitialConditions(temperature_profile->y()(j), inlet_P, inlet_omega);
			plug_flow_reactor->SetVerboseOutput(false);

			// Solve the plug flow reactor
			plug_flow_reactor->Solve(residence_time);
			
			// Store profiles
			for (unsigned int k = 0; k < thermodynamicsMapXML->NumberOfSpecies(); k++)
				omega_profiles_temp(j,k) = plug_flow_reactor->Y()(k);
		}

		std::vector<OpenSMOKE::FixedProfile*> omega_profiles(thermodynamicsMapXML->NumberOfSpecies());
		for (unsigned int k = 0; k < thermodynamicsMapXML->NumberOfSpecies(); k++)
			omega_profiles[k] = new OpenSMOKE::FixedProfile(temperature_profile->x().size(), 
															temperature_profile->x().data(), 
															omega_profiles_temp.col(k).data() );

		// Set heterogeneous mechanism
		CVI::HeterogeneousMechanism* heterogeneous_mechanism = new CVI::HeterogeneousMechanism(*thermodynamicsMapXML, *kineticsMapXML, *transportMapXML, dictionaries(dict_name_heterogeneous_mechanism));

		// Set heterogeneous mechanism
		CVI::HeterogeneousDetailedMechanism* heterogeneous_detailed_mechanism = new CVI::HeterogeneousDetailedMechanism(*thermodynamicsMapXML, *kineticsMapXML, *transportMapXML, *thermodynamicsSurfaceMapXML, *kineticsSurfaceMapXML, true, true);

		// Set porous medium
		CVI::PorousMedium* porous_medium = new CVI::PorousMedium(*thermodynamicsMapXML, *kineticsMapXML, *transportMapXML, dictionaries(dict_name_porous_medium));

		// Creates the reactor
		reactor1d = new CVI::Reactor1D(*thermodynamicsMapXML, *kineticsMapXML, *transportMapXML,
			*thermodynamicsSurfaceMapXML, *kineticsSurfaceMapXML,
			*porous_medium,
			*heterogeneous_mechanism, *heterogeneous_detailed_mechanism,
			*grid_x, detailed_heterogeneous_kinetics,
			SiteNonConservation, surface_dae_species, output_path);

		// Initial surface fractions
		Eigen::VectorXd initial_Z(thermodynamicsSurfaceMapXML->number_of_site_species());
		for (unsigned int i = 0; i < thermodynamicsSurfaceMapXML->number_of_site_species(); i++)
			initial_Z(i) = Z0[i + 1];

		reactor1d->SetPlanarSymmetry(symmetry_planar);
		reactor1d->SetInitialConditions(initial_T, initial_P, initial_omega, Gamma0, initial_Z);
		reactor1d->SetSiteNonConservation(SiteNonConservation);
		reactor1d->SetGasSide(temperature_profile, inlet_P, omega_profiles);
		reactor1d->SetTimeTotal(time_total);
		reactor1d->SetDaeTimeInterval(dae_time_interval);
		reactor1d->SetOdeEndTime(ode_end_time);
		reactor1d->SetDerivativeMassFractions(derivative_mass_fractions);
		reactor1d->SetDerivativeEffectiveDiffusivity(derivative_effective_diffusivity);
		reactor1d->SetDerivativeBulkDensity(derivative_bulk_density);

		if (on_the_fly_ropa == true)
			reactor1d->SetSurfaceOnTheFlyROPA(onTheFlyROPA);

		// Solve
		{
			time_t timerStart;
			time_t timerEnd;

			time(&timerStart);
			int flag = reactor1d->SolveFromScratch(*dae_parameters, *ode_parameters);
			time(&timerEnd);

			std::cout << "Total time: " << difftime(timerEnd, timerStart) << " s" << std::endl;
		}
	}

	// Capillary
	if (problem_type == CVI_CAPILLARY)
	{
		// Plug flow ractor simulation
		Eigen::VectorXd Y_gas_side(thermodynamicsMapXML->NumberOfSpecies());
		CVI::PlugFlowReactorCoupled* plug_flow_reactor = new CVI::PlugFlowReactorCoupled(*thermodynamicsMapXML, *kineticsMapXML, dictionaries(dict_name_plug_flow));
		{
			// Set initial conditions
			plug_flow_reactor->SetInitialConditions(inlet_T, inlet_P, inlet_omega);

			// Solve the plug flow reactor
			plug_flow_reactor->Solve(residence_time);
		}

		// Set heterogeneous mechanism
		CVI::HeterogeneousMechanism* heterogeneous_mechanism = new CVI::HeterogeneousMechanism(*thermodynamicsMapXML, *kineticsMapXML, *transportMapXML, dictionaries(dict_name_heterogeneous_mechanism));
		
		// Set heterogeneous mechanism
		CVI::HeterogeneousDetailedMechanism* heterogeneous_detailed_mechanism = new CVI::HeterogeneousDetailedMechanism(*thermodynamicsMapXML, *kineticsMapXML, *transportMapXML, *thermodynamicsSurfaceMapXML, *kineticsSurfaceMapXML, true, true);

		// Set capillary
		CVI::Capillary* capillary = new CVI::Capillary(	*thermodynamicsMapXML, *kineticsMapXML, *transportMapXML, 
														*thermodynamicsSurfaceMapXML, *kineticsSurfaceMapXML, 
														*heterogeneous_mechanism, *heterogeneous_detailed_mechanism, 
														*grid_x, detailed_heterogeneous_kinetics, 
														SiteNonConservation, surface_dae_species);

		// Initial surface fractions
		Eigen::VectorXd initial_Z(thermodynamicsSurfaceMapXML->number_of_site_species());
		for (unsigned int i = 0; i<thermodynamicsSurfaceMapXML->number_of_site_species(); i++)
			initial_Z(i) = Z0[i + 1];
		
		capillary->SetInitialConditions(initial_T, initial_P, capillary_diameter, initial_omega, Gamma0, initial_Z);
		capillary->SetGasSide(inlet_T, inlet_P, plug_flow_reactor->Y());
		capillary->SetTimeTotal(time_total);
		capillary->SetDaeTimeInterval(dae_time_interval);
		capillary->SetOdeEndTime(ode_end_time);
		capillary->SetTecplotTimeInterval(tecplot_time_interval);

		if (on_the_fly_ropa == true)
			capillary->SetSurfaceOnTheFlyROPA(onTheFlyROPA);

		// Solve
		{
			time_t timerStart;
			time_t timerEnd;

			time(&timerStart);
			int flag = capillary->SolveFromScratch(*dae_parameters, *ode_parameters);
			time(&timerEnd);

			std::cout << "Total time: " << difftime(timerEnd, timerStart) << " s" << std::endl;
		}
	}

	
	return 0;
}
