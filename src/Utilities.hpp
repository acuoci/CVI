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
|   This file is part of OpenSMOKE++ framework.                           |
|                                                                         |
|	License                                                           |
|                                                                         |
|   Copyright(C) 2014, 2013, 2012  Alberto Cuoci                          |
|   OpenSMOKE++ is free software: you can redistribute it and/or modify   |
|   it under the terms of the GNU General Public License as published by  |
|   the Free Software Foundation, either version 3 of the License, or     |
|   (at your option) any later version.                                   |
|                                                                         |
|   OpenSMOKE++ is distributed in the hope that it will be useful,        |
|   but WITHOUT ANY WARRANTY; without even the implied warranty of        |
|   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         |
|   GNU General Public License for more details.                          |
|                                                                         |
|   You should have received a copy of the GNU General Public License     |
|   along with OpenSMOKE++. If not, see <http://www.gnu.org/licenses/>.   |
|                                                                         |
\*-----------------------------------------------------------------------*/

void ReadFromBackupFile(const boost::filesystem::path path_file,
						double& t,
						Eigen::VectorXd& x, Eigen::VectorXd& y,
						Eigen::VectorXd& T, Eigen::VectorXd& P, Eigen::VectorXd& epsilon,
						std::vector< Eigen::VectorXd >& omega,
						std::vector< Eigen::VectorXd >& Z,
						std::vector< Eigen::VectorXd >& GammaFromEqs,
						std::vector<std::string>& gas_names_species,
						std::vector<std::string>& surface_names_species)
{
	rapidxml::xml_document<> xml_main_input;
	std::vector<char> local_xml_input_string;
	OpenSMOKE::OpenInputFileXML(xml_main_input, local_xml_input_string, path_file);

	unsigned int index_T;
	unsigned int index_P;
	unsigned int index_MW;
	unsigned int number_of_gaseous_species;
	unsigned int number_of_surface_species;
	unsigned int number_of_surface_phases;
	unsigned int number_of_additional_profiles;
	unsigned int index_epsilon = 4;

	// Time
	{
		rapidxml::xml_node<>* time_node = xml_main_input.first_node("opensmoke")->first_node("time");
		if (time_node != 0)
		{
			std::stringstream values(time_node->value());
			values >> t;
		}
		else
			OpenSMOKE::FatalErrorMessage("Corrupted backup xml file: missing the time leaf");
	}

	// Indices of T, P and MW
	{
		rapidxml::xml_node<>* indices_node = xml_main_input.first_node("opensmoke")->first_node("t-p-mw");
		if (indices_node != 0)
		{
			std::stringstream values(indices_node->value());
			values >> index_T;
			values >> index_P;
			values >> index_MW;
		}
		else
			OpenSMOKE::FatalErrorMessage("Corrupted backup xml file: missing the t-p-mw leaf");
	}

	// Additional
	{
		rapidxml::xml_node<>* additional_node = xml_main_input.first_node("opensmoke")->first_node("additional");
		if (additional_node != 0)
		{
			std::stringstream values(additional_node->value());
			values >> number_of_additional_profiles;
		}
		else
			OpenSMOKE::FatalErrorMessage("Corrupted backup xml file: missing the additional leaf");
	}

	// Gaseous species
	{
		rapidxml::xml_node<>* massfractions_node = xml_main_input.first_node("opensmoke")->first_node("mass-fractions");
		if (massfractions_node != 0)
		{
			std::stringstream values(massfractions_node->value());

			values >> number_of_gaseous_species;

			gas_names_species.resize(number_of_gaseous_species);
			for (unsigned int j = 0; j<number_of_gaseous_species; j++)
			{
				std::string dummy;
				values >> dummy;
				gas_names_species[j] = dummy;

				double dummy_double;
				double dummy_unsigned_int;
				values >> dummy_double;
				values >> dummy_unsigned_int;
			}
		}
		else
			OpenSMOKE::FatalErrorMessage("Corrupted backup xml file: missing the mass-fractions leaf");
	}

	// Surface species
	{
		rapidxml::xml_node<>* surfacefractions_node = xml_main_input.first_node("opensmoke")->first_node("surface-fractions");
		if (surfacefractions_node != 0)
		{
			std::stringstream values(surfacefractions_node->value());

			values >> number_of_surface_species;

			surface_names_species.resize(number_of_surface_species);
			for (unsigned int j = 0; j<number_of_surface_species; j++)
			{
				std::string dummy;
				values >> dummy;
				surface_names_species[j] = dummy;

				double dummy_double;
				double dummy_unsigned_int;
				values >> dummy_double;
				values >> dummy_unsigned_int;
			}
		}
		else
			OpenSMOKE::FatalErrorMessage("Corrupted backup xml file: missing the surface-fractions leaf");
	}

	// Surface phases
	{
		rapidxml::xml_node<>* surfacephases_node = xml_main_input.first_node("opensmoke")->first_node("surface-phases");
		if (surfacephases_node != 0)
		{
			std::stringstream values(surfacephases_node->value());
			values >> number_of_surface_phases;
		}
		else
			OpenSMOKE::FatalErrorMessage("Corrupted backup xml file: missing the surface-phases leaf");
	}

	// Read profiles
	std::vector< Eigen::VectorXd > additional;
	{
		unsigned int number_of_abscissas;
		unsigned int number_of_ordinates;

		{
			rapidxml::xml_node<>* profiles_size_node = xml_main_input.first_node("opensmoke")->first_node("profiles-size");

			if (profiles_size_node != 0)
			{
				std::stringstream values(profiles_size_node->value());
				values >> number_of_abscissas;
				values >> number_of_ordinates;
			}
			else
				OpenSMOKE::FatalErrorMessage("Corrupted backup xml file: missing the profiles-size leaf");
		}

		omega.resize(number_of_gaseous_species);
		for (unsigned int j = 0; j<number_of_gaseous_species; j++)
			omega[j].resize(number_of_abscissas);

		Z.resize(number_of_surface_species);
		for (unsigned int j = 0; j<number_of_surface_species; j++)
			Z[j].resize(number_of_abscissas);

		GammaFromEqs.resize(number_of_surface_phases);
		for (unsigned int j = 0; j<number_of_surface_phases; j++)
			GammaFromEqs[j].resize(number_of_abscissas);

		additional.resize(number_of_additional_profiles);
		for (unsigned int j = 0; j<number_of_additional_profiles; j++)
			additional[j].resize(number_of_abscissas);

		rapidxml::xml_node<>* profiles_node = xml_main_input.first_node("opensmoke")->first_node("profiles");
		if (profiles_node != 0)
		{
			std::stringstream values(profiles_node->value());
			for (unsigned int i = 0; i<number_of_abscissas; i++)
			{
				for (unsigned int j = 0; j<number_of_additional_profiles; j++)
					values >> additional[j](i);
				for (unsigned int j = 0; j<number_of_gaseous_species; j++)
					values >> omega[j](i);
				for (unsigned int j = 0; j<number_of_surface_species; j++)
					values >> Z[j](i);
				for (unsigned int j = 0; j<number_of_surface_phases; j++)
					values >> GammaFromEqs[j](i);
			}
		}
		else
			OpenSMOKE::FatalErrorMessage("Corrupted backup xml file: missing the profiles leaf");

		// Read grid (x direction)
		{
			rapidxml::xml_node<>* grid_x_node = xml_main_input.first_node("opensmoke")->first_node("grid-x");

			if (grid_x_node != 0)
			{
				unsigned int number_x_points;
				std::stringstream values(grid_x_node->value());
				values >> number_x_points;
				x.resize(number_x_points);
				for (unsigned int i = 0; i < number_x_points; i++)
					values >> x(i);
			}
			else
				OpenSMOKE::FatalErrorMessage("Corrupted backup xml file: missing the grid-x leaf");
		}

		// Read grid (y direction)
		{
			rapidxml::xml_node<>* grid_y_node = xml_main_input.first_node("opensmoke")->first_node("grid-y");

			if (grid_y_node != 0)
			{
				unsigned int number_y_points;
				std::stringstream values(grid_y_node->value());
				values >> number_y_points;
				y.resize(number_y_points);
				for (unsigned int i = 0; i < number_y_points; i++)
					values >> y(i);
			}
			else
				OpenSMOKE::FatalErrorMessage("Corrupted backup xml file: missing the grid-y leaf");
		}

		T = additional[index_T];
		P = additional[index_P];
		epsilon = additional[index_epsilon];
	}
}

void ReadFromBackupFile(const boost::filesystem::path path_file,
						double& t,
						OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN<double>& thermodynamicsSurfaceMap,
						const Eigen::VectorXd& x, const Eigen::VectorXd& y,
						Eigen::VectorXd& T, Eigen::VectorXd& P, Eigen::VectorXd& epsilon,
						std::vector< Eigen::VectorXd >& omega,
						std::vector< Eigen::VectorXd >& Z,
						std::vector< Eigen::VectorXd >& GammaFromEqs)
{
	std::vector < Eigen::VectorXd > omega_backup;
	std::vector < Eigen::VectorXd > Z_backup;
	std::vector < Eigen::VectorXd > GammaFromEqs_backup;
	std::vector<std::string> names_gas_species_backup;
	std::vector<std::string> names_surface_species_backup;

	Eigen::VectorXd x_backup;
	Eigen::VectorXd y_backup;
	Eigen::VectorXd T_backup;
	Eigen::VectorXd P_backup;
	Eigen::VectorXd epsilon_backup;

	ReadFromBackupFile(path_file, t, x_backup, y_backup, T_backup, P_backup, epsilon_backup, omega_backup, Z_backup, GammaFromEqs_backup, names_gas_species_backup, names_surface_species_backup);

	// Check if the mesh is the same
	if (x_backup.size() != x.size())
	{
		std::cout << "Recovering mesh from backup file: Fatal error!" << std::endl;
		OpenSMOKE::FatalErrorMessage("The grid along the x direction is not the same.");
	}
	if (y_backup.size() != y.size())
	{
		std::cout << "Recovering mesh from backup file: Fatal error!" << std::endl;
		OpenSMOKE::FatalErrorMessage("The grid along the y direction is not the same.");
	}

	// Check if the gaseous kinetic mechanism is the same
	for (unsigned int i = 0; i < thermodynamicsSurfaceMap.number_of_gas_species(); i++)
		if (thermodynamicsSurfaceMap.NamesOfSpecies()[i] != names_gas_species_backup[i])
		{
			std::cout << "Recovering gaseous species names from backup file: Fatal error!" << std::endl;
			OpenSMOKE::FatalErrorMessage("The gaseous kinetic mechanism is not the same.");
		}

	// Check if the surface kinetic mechanism is the same
	for (unsigned int i = 0; i < thermodynamicsSurfaceMap.number_of_site_species(); i++)
		if (thermodynamicsSurfaceMap.NamesOfSpecies()[i + thermodynamicsSurfaceMap.number_of_gas_species()] != names_surface_species_backup[i])
		{
			std::cout << "Recovering surface species names from backup file: Fatal error!" << std::endl;
			OpenSMOKE::FatalErrorMessage("The surface kinetic mechanism is not the same.");
		}

	// Copying the fields
	for (int k = 0; k < y.size(); k++)
		for (int i = 0; i < x.size(); i++)
		{
			const int point = k*x.size() + i;

			for (unsigned int j = 0; j < thermodynamicsSurfaceMap.number_of_gas_species(); j++)
				omega[point](j) = omega_backup[j](point);

			for (unsigned int j = 0; j < thermodynamicsSurfaceMap.number_of_site_species(); j++)
				Z[point](j) = Z_backup[j](point);

			for (unsigned int j = 0; j < thermodynamicsSurfaceMap.number_of_site_phases(0); j++)
				GammaFromEqs[point](j) = GammaFromEqs_backup[j](point);
		}

	T = T_backup;
	P = P_backup;
	epsilon = epsilon_backup;
}
