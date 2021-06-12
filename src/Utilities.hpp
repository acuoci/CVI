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

#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/foreach.hpp>

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
	boost::property_tree::ptree ptree;
    	boost::property_tree::read_xml( (path_file).string(), ptree );

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
		std::cout << " * From backup file: reading time... " << std::endl;
		boost::optional< boost::property_tree::ptree& > child = ptree.get_child_optional("opensmoke.time");
		if (child)
			t = ptree.get<double>("opensmoke.time");
		else
			OpenSMOKE::FatalErrorMessage("Corrupted backup xml file: missing the time leaf");
	}

	// Indices of T, P and MW
	{
		std::cout << " * From backup file: reading t-p-mw... " << std::endl;
		boost::optional< boost::property_tree::ptree& > child = ptree.get_child_optional("opensmoke.t-p-mw");
		if (child)
		{
			std::stringstream values;
			values.str( ptree.get< std::string >("opensmoke.t-p-mw") );
			values >> index_T;
			values >> index_P;
			values >> index_MW;
		}
		else
			OpenSMOKE::FatalErrorMessage("Corrupted backup xml file: missing the t-p-mw leaf");
	}

	// Additional
	{
		std::cout << " * From backup file: reading additional... " << std::endl;
		boost::optional< boost::property_tree::ptree& > child = ptree.get_child_optional("opensmoke.additional");
		if (child)
		{
			std::stringstream values;
			values.str( ptree.get< std::string >("opensmoke.additional") );
			values >> number_of_additional_profiles;
		}
		else
			OpenSMOKE::FatalErrorMessage("Corrupted backup xml file: missing the additional leaf");
	}

	// Gaseous species
	{
		std::cout << " * From backup file: reading mass-fractions... " << std::endl;
		boost::optional< boost::property_tree::ptree& > child = ptree.get_child_optional("opensmoke.mass-fractions");
		if (child)
		{
			std::stringstream values;
			values.str( ptree.get< std::string >("opensmoke.mass-fractions") );
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
		std::cout << " * From backup file: reading surface-fractions... " << std::endl;
		boost::optional< boost::property_tree::ptree& > child = ptree.get_child_optional("opensmoke.surface-fractions");
		if (child)
		{
			std::stringstream values;
			values.str( ptree.get< std::string >("opensmoke.surface-fractions") );
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
		std::cout << " * From backup file: reading surface-phases... " << std::endl;
		boost::optional< boost::property_tree::ptree& > child = ptree.get_child_optional("opensmoke.surface-phases");
		if (child)
		{
			std::stringstream values;
			values.str( ptree.get< std::string >("opensmoke.surface-phases") );
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
			std::cout << " * From backup file: reading profiles-size... " << std::endl;
			boost::optional< boost::property_tree::ptree& > child = ptree.get_child_optional("opensmoke.profiles-size");
			if (child)
			{
				std::stringstream values;
				values.str( ptree.get< std::string >("opensmoke.profiles-size") );
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

		std::cout << " * From backup file: reading profiles... " << std::endl;
		boost::optional< boost::property_tree::ptree& > child = ptree.get_child_optional("opensmoke.profiles");
		if (child)
		{
			std::stringstream values;
			values.str( ptree.get< std::string >("opensmoke.profiles") );
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
			std::cout << " * From backup file: reading grid-x... " << std::endl;
			boost::optional< boost::property_tree::ptree& > child = ptree.get_child_optional("opensmoke.grid-x");

			if (child)
			{
				unsigned int number_x_points;
				std::stringstream values;
				values.str( ptree.get< std::string >("opensmoke.grid-x") );
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
			std::cout << " * From backup file: reading grid-y... " << std::endl;
			boost::optional< boost::property_tree::ptree& > child = ptree.get_child_optional("opensmoke.grid-y");

			if (child)
			{
				unsigned int number_y_points;
				std::stringstream values;
				values.str( ptree.get< std::string >("opensmoke.grid-y") );
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
						OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN& thermodynamicsSurfaceMap,
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

	const double tol_eps_error = 0.05;
	const double tol_eps_warning = 0.01;

	// Check and adjust the species mass fractions
	for (int k = 0; k < y.size(); k++)
		for (int i = 0; i < x.size(); i++)
		{
			const int point = k * x.size() + i;

			for (unsigned int j = 0; j < thermodynamicsSurfaceMap.number_of_gas_species(); j++)
				omega[point](j) = std::max(0., omega[point](j));

			for (unsigned int j = 0; j < thermodynamicsSurfaceMap.number_of_gas_species(); j++)
				omega[point](j) = std::min(1., omega[point](j));

			double sum = 0.;
			for (unsigned int j = 0; j < thermodynamicsSurfaceMap.number_of_gas_species(); j++)
				sum += omega[point](j);

			if (std::fabs(sum - 1.) > tol_eps_error)
			{
				std::cout << "The sum of mass fractions in point " << i << "," << k << " is outside the acceptable tolerance" << std::endl;
				std::cout << "Acceptable tolerance: " << tol_eps_error << " - Current sum: " << sum << std::endl;
				OpenSMOKE::FatalErrorMessage("Please check your backup file");
			}

			if (std::fabs(sum - 1.) > tol_eps_warning)
				std::cout << "WARNING: The sum of mass fractions in point " << i << "," << k << " is " << sum << std::endl;

			for (unsigned int j = 0; j < thermodynamicsSurfaceMap.number_of_gas_species(); j++)
				omega[point](j) /= sum;
		}

	// Check and adjust the species mass fractions
	for (int k = 0; k < y.size(); k++)
		for (int i = 0; i < x.size(); i++)
		{
			const int point = k * x.size() + i;

			for (unsigned int j = 0; j < thermodynamicsSurfaceMap.number_of_site_species(); j++)
				Z[point](j) = std::max(0., Z[point](j));

			for (unsigned int j = 0; j < thermodynamicsSurfaceMap.number_of_site_species(); j++)
				Z[point](j) = std::min(1., Z[point](j));

			double sum = 0.;
			for (unsigned int j = 0; j < thermodynamicsSurfaceMap.number_of_site_species(); j++)
				sum += Z[point](j);

			if (std::fabs(sum - 1.) > tol_eps_error)
			{
				std::cout << "The sum of surface species fractions in point " << i << "," << k << " is outside the acceptable tolerance" << std::endl;
				std::cout << "Acceptable tolerance: " << tol_eps_error << " - Current sum: " << sum << std::endl;
				OpenSMOKE::FatalErrorMessage("Please check your backup file");
			}

			if (std::fabs(sum - 1.) > tol_eps_warning)
				std::cout << "WARNING: The sum of surface species fractions in point " << i << "," << k << " is " << sum << std::endl;

			for (unsigned int j = 0; j < thermodynamicsSurfaceMap.number_of_site_species(); j++)
				Z[point](j) /= sum;
		}


	T = T_backup;
	P = P_backup;
	epsilon = epsilon_backup;
}
