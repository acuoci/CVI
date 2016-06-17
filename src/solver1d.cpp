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

// Grammars
#include "Grammar_CVI_Solver1D.h"
#include "Grammar_CVI_PlugFlowReactor.h"
#include "Grammar_CVI_PorousMedium.h"

// Porous medium
#include "PorousMedium.h"

// Plug Flow Reactor
#include "PlugFlowReactor.h"
#include "PlugFlowReactorProfiles.h"

// Reactor 1D
#include "Reactor1D.h"

// Reactor 2D
#include "Reactor2D.h"

// Numerical parameters
#include "math/multivalue-dae-solvers/parameters/DaeSolver_Parameters.h"

#include "math/nls-solvers/NonLinearSolver.h"
#include "math/nls-solvers/KernelSparse.h"

class BrandBrandNewProblem01tris
{
public:

	BrandBrandNewProblem01tris() {};


protected:

	unsigned int ne_;

	void MemoryAllocation()
	{
	}

	virtual void Equations(const Eigen::VectorXd &x, Eigen::VectorXd &f)
	{
		f(0) = x(1) - x(0) - 1.;
		f(1) = x(1)*x(2)*x(3) - 6.;
		f(2) = x(2) + x(1)*x(3) + x(4) - 9.;
		f(3) = x(8)*x(2) + x(3) - 19.;
		f(4) = x(4)*x(1) + x(6) - 10.;
		f(5) = x(5)*x(1)*x(9) - 45.;
		f(6) = x(6)*x(3) - 18.;
		f(7) = x(7) + x(2) - 9.;
		f(8) = x(8)*x(1) + x(9)*x(2) - 26.;
		f(9) = x(9)*x(4)*x(1) - 36.;
	}

	void Jacobian(const Eigen::VectorXd &x, Eigen::SparseMatrix<double> &J)
	{
	};

	void Print(const int call, const double phiW, const Eigen::VectorXd &x, const Eigen::VectorXd &f)
	{
		std::cout << "Solution: " << call << " " << phiW << std::endl;
	}
};

int main(int argc, char** argv)
{
	/*
	{
		// Brand New
		{
			// Problem size
			int N = 10;

			// Memory Allocation and first guess solution
			//BzzVector y0(N, -0.5, 0.1, 1.2, -4.3, -2.4, -7.5, 5.6, -4.7, 10.8, -5.9);
			BzzVector y0(N); y0 = 1.;
			Eigen::VectorXd y0_eigen(N);

			BzzVector y(N);
			Eigen::VectorXd y_eigen(N);

			BzzVector f(N);
			Eigen::VectorXd f_eigen(N);

			// Conversion
			for (int i = 0; i < N; i++)
			{
				y0_eigen(i) = y0[i + 1];
			}

			std::cout << "Problem01tris (10x10): solved as sparse: OS++ new" << std::endl;
			getchar();

			std::vector<unsigned int> rows(29);
			std::vector<unsigned int> cols(29);
			unsigned int count = 0;
			rows[count] = 1;	cols[count] = 1; count++;
			rows[count] = 1;	cols[count] = 2; count++;
			rows[count] = 2;	cols[count] = 2; count++;
			rows[count] = 2;	cols[count] = 3; count++;
			rows[count] = 2;	cols[count] = 4; count++;
			rows[count] = 3;	cols[count] = 2; count++;
			rows[count] = 3;	cols[count] = 3; count++;
			rows[count] = 3;	cols[count] = 4; count++;
			rows[count] = 3;	cols[count] = 5; count++;
			rows[count] = 4;	cols[count] = 3; count++;
			rows[count] = 4;	cols[count] = 4; count++;
			rows[count] = 4;	cols[count] = 9; count++;
			rows[count] = 5;	cols[count] = 2; count++;
			rows[count] = 5;	cols[count] = 5; count++;
			rows[count] = 5;	cols[count] = 7; count++;
			rows[count] = 6;	cols[count] = 2; count++;
			rows[count] = 6;	cols[count] = 6; count++;
			rows[count] = 6;	cols[count] = 10; count++;
			rows[count] = 7;	cols[count] = 4; count++;
			rows[count] = 7;	cols[count] = 7; count++;
			rows[count] = 8;	cols[count] = 3; count++;
			rows[count] = 8;	cols[count] = 8; count++;
			rows[count] = 9;	cols[count] = 2; count++;
			rows[count] = 9;	cols[count] = 3; count++;
			rows[count] = 9;	cols[count] = 9; count++;
			rows[count] = 9;	cols[count] = 10; count++;
			rows[count] = 10;	cols[count] = 2; count++;
			rows[count] = 10;	cols[count] = 5; count++;
			rows[count] = 10;	cols[count] = 10; count++;

			typedef NlsSMOKE::KernelSparse< BrandBrandNewProblem01tris> kernel;
			NlsSMOKE::NonLinearSolver<kernel> nls;
			nls.SetFirstGuessSolution(y0_eigen);
			//nls.SetFullPivoting(true);
			nls.SetSparsityPattern(rows, cols, false);
			nls.SetPreconditioner("ILUT");
			nls.SetAbsoluteTolerances(1e-12);
			nls.SetRelativeTolerances(1e-7);
			nls.SetPrint(true);

			const double start = BzzGetCpuTime();
			nls();
			const double end = BzzGetCpuTime();
			nls.Solution(y_eigen, f_eigen);

			std::cout << end - start << std::endl;

			nls.NlsSummary(std::cout);
		}
	}
	getchar();
	*/
	
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
	CVI::Grammar_CVI_Solver1D			grammar_cvi_solver1d;
	CVI::Grammar_CVI_PlugFlowReactor	grammar_cvi_plug_flow_reactor;
	CVI::Grammar_CVI_PorousMedium		grammar_cvi_porous_medium;
	CVI::Grammar_Defect_PorousMedium	grammar_defect_porous_medium;

	// Define the dictionaries
	OpenSMOKE::OpenSMOKE_DictionaryManager dictionaries;
	dictionaries.ReadDictionariesFromFile(input_file_name_);
	dictionaries(main_dictionary_name_).SetGrammar(grammar_cvi_solver1d);

	// Plug flow reactor
	std::string dict_name_plug_flow;
	if (dictionaries(main_dictionary_name_).CheckOption("@PlugFlowReactor") == true)
		dictionaries(main_dictionary_name_).ReadDictionary("@PlugFlowReactor", dict_name_plug_flow);

	// Porous medium
	std::string dict_name_porous_medium;
	if (dictionaries(main_dictionary_name_).CheckOption("@PorousMedium") == true)
		dictionaries(main_dictionary_name_).ReadDictionary("@PorousMedium", dict_name_porous_medium);

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
	dictionaries(dict_name_plug_flow).SetGrammar(grammar_cvi_plug_flow_reactor);
	dictionaries(dict_name_porous_medium).SetGrammar(grammar_cvi_porous_medium);

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

	if (dictionaries(main_dictionary_name_).CheckOption("@Output") == true)
	{
		boost::filesystem::path output_path;
		dictionaries(main_dictionary_name_).ReadPath("@Output", output_path);
		OpenSMOKE::CreateDirectory(output_path);
	}

	// Type of problem
	bool problem_2d = true;
	{
		std::string value;
		if (dictionaries(main_dictionary_name_).CheckOption("@Type") == true)
			dictionaries(main_dictionary_name_).ReadString("@Type", value);
		if (value == "1D")			problem_2d = false;
		else if (value == "2D")		problem_2d = true;
		else OpenSMOKE::FatalErrorMessage("Wrong @Type: 1D | 2D");
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

			grid_x = new OpenSMOKE::Grid1D(number_of_points, 0., length, stretching_factor);
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

			grid_x = new OpenSMOKE::Grid1D(number_of_points, radius_internal, radius_external, stretching_factor);
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

		grid_y = new OpenSMOKE::Grid1D(number_of_points, 0., length, stretching_factor);
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

	// Read the plug flow residence time 
	double residence_time = 1.;
	if (problem_2d == false)
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
			OpenSMOKE::FatalErrorMessage("@ResidenceTime must be specified");
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
		for (unsigned int i = 1; i <= aux.Size(); i++)
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
		for (unsigned int i = 1; i <= aux.Size(); i++)
			initial_omega(i - 1) = aux[i];
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

	// Set porous medium
	CVI::PorousMedium* porous_medium = new CVI::PorousMedium(*thermodynamicsMapXML, *kineticsMapXML, *transportMapXML, dictionaries(dict_name_porous_medium));

	// Solve the 2D problem
	if (problem_2d == true)
	{
		// Plug flow ractor simulation
		std::vector<Eigen::VectorXd> Y_gas_side;
		CVI::PlugFlowReactor* plug_flow_reactor = new CVI::PlugFlowReactor(*thermodynamicsMapXML, *kineticsMapXML, dictionaries(dict_name_plug_flow));
		{
			// Set initial conditions
			plug_flow_reactor->SetInitialConditions(inlet_T, inlet_P, inlet_omega);

			// Solve the plug flow reactor
			plug_flow_reactor->Solve(residence_time);

			// Extract the plug flow history
			Eigen::VectorXd csi(plug_flow_reactor->history_csi().size());
			for (unsigned int i = 0; i < plug_flow_reactor->history_csi().size(); i++)
				csi(i) = plug_flow_reactor->history_csi()[i];

			// Assign the boundary conditions
			PlugFlowReactorProfiles* profiles = new PlugFlowReactorProfiles(csi);
			Y_gas_side.resize(2 * grid_x->Np() + grid_y->Np());
			for (int i = 0; i < Y_gas_side.size(); i++)
			{
				Y_gas_side[i].resize(thermodynamicsMapXML->NumberOfSpecies());
				Y_gas_side[i].setZero();
			}

			// Case: only east side
			if (plug_flow_reactor->geometric_pattern() == CVI::PlugFlowReactor::ONE_SIDE)
			{
				for (int i = 0; i < grid_y->Np(); i++)
				{
					const int point = i + grid_x->Np();
					const double coordinate = grid_y->x()(i)+plug_flow_reactor->inert_length();
					profiles->Interpolate(coordinate, plug_flow_reactor->history_Y(), Y_gas_side[point]);
				}
			}
			// Case: south, east, and north sides
			else if (plug_flow_reactor->geometric_pattern() == CVI::PlugFlowReactor::THREE_SIDES)
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
		}


		// Creates the reactor
		reactor2d = new CVI::Reactor2D(*thermodynamicsMapXML, *kineticsMapXML, *transportMapXML, *porous_medium, *porosity_defect, *grid_x, *grid_y, *plug_flow_reactor);

		// Set options
		reactor2d->SetPlanarSymmetry(symmetry_planar);
		reactor2d->SetInitialConditions(initial_T, initial_P, initial_omega);
		reactor2d->SetGasSide(inlet_T, inlet_P, Y_gas_side);
		reactor2d->SetTimeTotal(time_total);
		reactor2d->SetDaeTimeInterval(dae_time_interval);
		reactor2d->SetTecplotTimeInterval(tecplot_time_interval);
		if (steps_video>0)	reactor2d->SetStepsVideo(steps_video);
		if (steps_file>0)	reactor2d->SetStepsFile(steps_file);

		// Solve
		{
			time_t timerStart;
			time_t timerEnd;

			time(&timerStart);
			int flag = reactor2d->SolveFromScratch(*dae_parameters);
			time(&timerEnd);

			std::cout << "Total time: " << difftime(timerEnd, timerStart) << " s" << std::endl;
		}
	}
	

	// Solve the 1D problem
	if (problem_2d == false)
	{
		// Plug flow ractor simulation
		Eigen::VectorXd Y_gas_side(thermodynamicsMapXML->NumberOfSpecies());
		CVI::PlugFlowReactor* plug_flow_reactor = new CVI::PlugFlowReactor(*thermodynamicsMapXML, *kineticsMapXML, dictionaries(dict_name_plug_flow));
		{
			// Set initial conditions
			plug_flow_reactor->SetInitialConditions(inlet_T, inlet_P, inlet_omega);

			// Solve the plug flow reactor
			plug_flow_reactor->Solve(residence_time);
		}

		CVI::Reactor1D* reactor1d = new CVI::Reactor1D(*thermodynamicsMapXML, *kineticsMapXML, *transportMapXML, *porous_medium, *grid_x);

		reactor1d->SetPlanarSymmetry(symmetry_planar);
		reactor1d->SetInitialConditions(initial_T, initial_P, initial_omega);
		reactor1d->SetGasSide(inlet_T, inlet_P, plug_flow_reactor->Y());

		// Solve
		{
			time_t timerStart;
			time_t timerEnd;

			time(&timerStart);
			int flag = reactor1d->SolveFromScratch(*dae_parameters);
			time(&timerEnd);

			std::cout << "Total time: " << difftime(timerEnd, timerStart) << " s" << std::endl;
		}
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
