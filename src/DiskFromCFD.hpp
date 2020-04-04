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

#include "Utilities.h"

namespace CVI
{
	void Interpolate(const double x_req, double& temperature_req, std::vector<double>& mass_fractions_req,
		const std::vector<double>& coordinates_from_cfd,
		const std::vector<double>& temperature_from_cfd,
		const std::vector< std::vector<double> >& mass_fractions_from_cfd);

	DiskFromCFD::DiskFromCFD(	OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap,
								OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap,
								OpenSMOKE::Grid1D& grid_x,
								OpenSMOKE::Grid1D& grid_y) :
	thermodynamicsMap_(thermodynamicsMap),
	kineticsMap_(kineticsMap),
	grid_x_(grid_x),
	grid_y_(grid_y)
	{
		ri_ = grid_x_.x()[0];
		re_ = grid_x_.x()[grid_x_.Np()-1];
		H_ = grid_y_.x()[grid_y_.Np() - 1];

		hole_ = true;
		if (std::fabs(ri_) <= 1e-12)
			hole_ = false;

		radial_coordinate_ = 0;
		axial_coordinate_ = 1;
	}

	void DiskFromCFD::ReadFromFile(const boost::filesystem::path disk_file_name)
	{
		std::cout << "Reading disk from file..." << std::endl;

		bool edge_centered_policy = false;

		rapidxml::xml_document<> doc;
		std::vector<char> xml_string;
		OpenInputFileXML(doc, xml_string, disk_file_name);

		rapidxml::xml_node<>* opensmoke_node = doc.first_node("opensmoke");
		rapidxml::xml_node<>* number_of_species_node = opensmoke_node->first_node("NumberOfSpecies");
		rapidxml::xml_node<>* names_of_species_node = opensmoke_node->first_node("NamesOfSpecies");

		std::cout << " * Check axial coordinate..." << std::endl;
		{
			rapidxml::xml_node<>* axial_coordinate_node = opensmoke_node->first_node("AxialCoordinate");
			if (axial_coordinate_node != 0)
			{
				const std::string coordinate = boost::trim_copy(std::string(axial_coordinate_node->value()));
				if (coordinate == "x")			axial_coordinate_ = 0;
				else if (coordinate == "y")		axial_coordinate_ = 1;
				else if (coordinate == "z")		axial_coordinate_ = 2;
				else OpenSMOKE::ErrorMessage("DiskFromCFD::ReadFromFile", "Wrong AxialCoordinate: x | y (default) | z");
			}
		}

		std::cout << " * Check radial coordinate..." << std::endl;
		{
			rapidxml::xml_node<>* radial_coordinate_node = opensmoke_node->first_node("RadialCoordinate");
			if (radial_coordinate_node != 0)
			{
				const std::string coordinate = boost::trim_copy(std::string(radial_coordinate_node->value()));
				if (coordinate == "x")			radial_coordinate_ = 0;
				else if (coordinate == "y")		radial_coordinate_ = 1;
				else if (coordinate == "z")		radial_coordinate_ = 2;
				else OpenSMOKE::ErrorMessage("DiskFromCFD::ReadFromFile", "Wrong RadialCoordinate: x (default) | y | z");
			}
		}

		std::cout << " * Check policy..." << std::endl;
		rapidxml::xml_node<>* policy_node = opensmoke_node->first_node("Policy");
		if (policy_node != 0)
		{
			const std::string policy = boost::trim_copy(std::string(policy_node->value()));
			if (policy == "edge-centered")		edge_centered_policy = true;
			else if (policy == "face-centered")	edge_centered_policy = false;
			else OpenSMOKE::ErrorMessage("DiskFromCFD::ReadFromFile", "Wrong policy: face-centered (default) | edge-centered)");
		}

		try
		{
			std::cout << " * Check kinetic mechanism..." << std::endl;
			const unsigned int ns = boost::lexical_cast<unsigned int>(boost::trim_copy(std::string(number_of_species_node->value())));
			
			if (ns != thermodynamicsMap_.NumberOfSpecies())
				OpenSMOKE::ErrorMessage("DiskFromCFD::ReadFromFile", "The kinetic mechanism adopted in the exported Disk data does not match with the current kinetic mechanism");

			std::vector<std::string> names_species(ns);

			std::stringstream names_of_species_string;
			names_of_species_string << names_of_species_node->value();
			for (unsigned int i = 0; i < ns; i++)
				names_of_species_string >> names_species[i];

			for (unsigned int i = 0; i<ns; i++)
				if (names_species[i] != thermodynamicsMap_.NamesOfSpecies()[i])
					OpenSMOKE::ErrorMessage("DiskFromCFD::ReadFromFile", "The kinetic mechanism adopted in the exported Disk data does not match with the current kinetic mechanism");
		}
		catch (...)
		{
			OpenSMOKE::ErrorMessage("DiskFromCFD::ReadFromFile", "Error in reading the list of species.");
		}

		// Corners
		double north_coordinate = 0.;
		double south_coordinate = 0.;
		double east_coordinate = 0.;
		double west_coordinate = 0.;

		//rapidxml::xml_node<>* data_node = opensmoke_node->first_node("Data");
		for (rapidxml::xml_node<> *child = opensmoke_node->first_node("Data"); child; child = child->next_sibling())
		{
			rapidxml::xml_attribute<> *side = child->first_attribute();

			std::cout << " * Importing disk data from side " << side->value() << std::endl;

			std::string side_name;
			std::stringstream side_name_string; side_name_string << side->value();
			side_name_string >> side_name;

			rapidxml::xml_node<>* number_points_node = child->first_node("NumberOfPoints");
			rapidxml::xml_node<>* points_node = child->first_node("Points");
			rapidxml::xml_node<>* mass_fractions_node = child->first_node("MassFractions");
			rapidxml::xml_node<>* temperature_node = child->first_node("Temperature");

			// Number of points
			const unsigned int np = boost::lexical_cast<unsigned int>(boost::trim_copy(std::string(number_points_node->value())));
			std::vector<double> coordinates_from_cfd(np);
			std::vector<double> temperature_from_cfd(np);
			std::vector< std::vector<double> > mass_fractions_from_cfd(thermodynamicsMap_.NumberOfSpecies());
			for (unsigned int i = 0; i < thermodynamicsMap_.NumberOfSpecies(); i++)
				mass_fractions_from_cfd[i].resize(np);

			// Data from CFD
			std::vector<size_t> indices_increasing(np);
		
			// Point coordinates
			{
				std::stringstream points_xml;
				points_xml << points_node->value();

				std::vector<Eigen::VectorXd> coordinates(3);
				for (unsigned int i = 0; i < 3; i++)
					coordinates[i].resize(np);

				for (unsigned int i = 0; i < np; i++)
				{
					points_xml >> coordinates[0](i);
					points_xml >> coordinates[1](i);
					points_xml >> coordinates[2](i);
				}

				Eigen::VectorXd mean(3);
				mean.setZero();
				for (unsigned int i = 0; i < 3; i++)
				{
					
					for (unsigned int j = 0; j < np; j++)
						mean(i) += coordinates[i](j);
					mean(i) /= static_cast<double>(np);
				}

				Eigen::VectorXd change(3);
				change.setZero();
				for (unsigned int i = 0; i < 3; i++)
				{
					for (unsigned int j = 0; j < np; j++)
						change(i) += std::fabs(mean(i)-coordinates[i](j));
					change(i) /= np;
				}

				for (unsigned int i = 0; i < 3; i++)
					std::cout << " * Direction " << i << ": " << " Mean(mm): " << mean(i)*1000. << " Change(mm): " << change(i)*1000. << std::endl;

				// Check alignement
				const double tolerance = 1.e-5;
				{
					unsigned int count_zeros = 0;
					for (unsigned int i = 0; i < 3; i++)
						if (std::fabs(change(i)) < tolerance)	count_zeros++;
					if (count_zeros != 2)
						OpenSMOKE::ErrorMessage("DiskFromCFD::ReadFromFile", "The current side is not aligned with none of the cartisian axes");
				}

				// Find the relevant coordinate
				int relevant_coordinate = -1;
				for (unsigned int i = 0; i < 3; i++)
					if (std::fabs(change(i)) > tolerance)
					{
						relevant_coordinate = i;
						break;
					}
				std::cout << " * Relevant direction: " << relevant_coordinate << std::endl;

				// Assign coordinates
				std::vector<double> coordinates_from_cfd_unsorted(np);
				for (unsigned int i = 0; i < np; i++)
					coordinates_from_cfd_unsorted[i] = coordinates[relevant_coordinate](i);

				indices_increasing = OpenSMOKE::SortAndTrackIndicesIncreasing(coordinates_from_cfd_unsorted);

				for (unsigned int i = 0; i < np; i++)
					coordinates_from_cfd[i] = coordinates_from_cfd_unsorted[indices_increasing[i]];

				if (edge_centered_policy == true)
				{
					coordinates_from_cfd[0]    = 0.9*coordinates_from_cfd[0] + 0.10*coordinates_from_cfd[1];
					coordinates_from_cfd[np-1] = 0.9*coordinates_from_cfd[np-1] + 0.10*coordinates_from_cfd[np-2];
				}

				if (side_name == "North")
				{
					if (relevant_coordinate != radial_coordinate_)
						OpenSMOKE::ErrorMessage("DiskFromCFD::ReadFromFile", "North side is not aligned with radial coordinate specified in XML file");

					north_coordinate = coordinates[axial_coordinate_](0);
					std::cout << " * North coordinate (mm): " << north_coordinate*1000. << std::endl;
				}
				else if (side_name == "South")
				{
					if (relevant_coordinate != radial_coordinate_)
						OpenSMOKE::ErrorMessage("DiskFromCFD::ReadFromFile", "South side is not aligned with radial coordinate specified in XML file");

					south_coordinate = coordinates[axial_coordinate_](0);
					std::cout << " * South coordinate (mm): " << south_coordinate * 1000. << std::endl;
				}
				else if (side_name == "East")
				{
					if (relevant_coordinate != axial_coordinate_)
						OpenSMOKE::ErrorMessage("DiskFromCFD::ReadFromFile", "East side is not aligned with axial coordinate specified in XML file");

					east_coordinate = coordinates[radial_coordinate_](0);
					std::cout << " * East coordinate (mm): " << east_coordinate * 1000. << std::endl;
				}
				else if (side_name == "West")
				{
					if (relevant_coordinate != axial_coordinate_)
						OpenSMOKE::ErrorMessage("DiskFromCFD::ReadFromFile", "West side is not aligned with axial coordinate specified in XML file");

					west_coordinate = coordinates[radial_coordinate_](0);
					std::cout << " * West coordinate (mm): " << west_coordinate * 1000. << std::endl;
				}
			}

			// Temperature
			{
				std::vector<double> temperature_from_cfd_unsorted(np);

				std::stringstream temperature_xml;
				temperature_xml << temperature_node->value();

				for (unsigned int i = 0; i < np; i++)
					temperature_xml >> temperature_from_cfd_unsorted[i];

				for (unsigned int i = 0; i < np; i++)
					temperature_from_cfd[i] = temperature_from_cfd_unsorted[indices_increasing[i]];
			}

			// Mass fractions
			{
				std::vector< std::vector<double> > mass_fractions_from_cfd_unsorted(thermodynamicsMap_.NumberOfSpecies());
				for (unsigned int i = 0; i < thermodynamicsMap_.NumberOfSpecies(); i++)
					mass_fractions_from_cfd_unsorted[i].resize(np);

				std::stringstream mass_fractions_xml;
				mass_fractions_xml << mass_fractions_node->value();

				for (unsigned int i = 0; i < np; i++)
					for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
						mass_fractions_xml >> mass_fractions_from_cfd_unsorted[j][i];

				for (unsigned int i = 0; i < np; i++)
					for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
						mass_fractions_from_cfd[j][i] = mass_fractions_from_cfd_unsorted[j][indices_increasing[i]];
			}

			// Extend data
			{

				if (side_name == "North" || side_name == "South")
				{
					std::cout << "Extending" << std::endl;

					if (coordinates_from_cfd[0] <= ri_ || coordinates_from_cfd[np - 1] >= re_)
					{
						std::cout << "Provided radii (mm): " << ri_*1000. << " " << re_*1000. << std::endl;
						std::cout << "Coordinates from CFD (mm): " << coordinates_from_cfd[0]*1000. << " " << coordinates_from_cfd[np - 1]*1000. << std::endl;
						OpenSMOKE::ErrorMessage("DiskFromCFD::ReadFromFile", "The size of the disk from CFD does not match the input size");
					}

					// Initial value
					coordinates_from_cfd.insert(coordinates_from_cfd.begin(), ri_*(0.9999999));
					temperature_from_cfd.insert(temperature_from_cfd.begin(), temperature_from_cfd[0]);
					for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
						mass_fractions_from_cfd[j].insert(mass_fractions_from_cfd[j].begin(), mass_fractions_from_cfd[j][0]);

					// Final value
					coordinates_from_cfd.push_back(re_*(1.0000001));
					temperature_from_cfd.push_back(temperature_from_cfd[np]);
					for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
						mass_fractions_from_cfd[j].push_back(mass_fractions_from_cfd[j][np]);
				}
			}

			if (side_name == "East" || side_name == "West")
			{
				for (unsigned int i = 0; i < np; i++)
					coordinates_from_cfd[i] -= south_coordinate;

				if (coordinates_from_cfd[0] <= 0. || coordinates_from_cfd[np - 1] >= H_)
				{
					std::cout << "Provided heights (mm): " << 0. << " " << H_*1000. << std::endl;
					std::cout << "Coordinates from CFD (mm): " << coordinates_from_cfd[0]*1000. << " " << coordinates_from_cfd[np - 1]*1000. << std::endl;
					OpenSMOKE::ErrorMessage("DiskFromCFD::ReadFromFile", "The size of the disk from CFD does not match the input size");
				}

				// Initial value
				coordinates_from_cfd.insert(coordinates_from_cfd.begin(), -0.0000001);
				temperature_from_cfd.insert(temperature_from_cfd.begin(), temperature_from_cfd[0]);
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					mass_fractions_from_cfd[j].insert(mass_fractions_from_cfd[j].begin(), mass_fractions_from_cfd[j][0]);

				// Final value
				coordinates_from_cfd.push_back(H_*(1.0000001));
				temperature_from_cfd.push_back(temperature_from_cfd[np]);
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					mass_fractions_from_cfd[j].push_back(mass_fractions_from_cfd[j][np]);
			}

			north_temperature_.resize(grid_x_.Np());
			north_mass_fractions_.resize(grid_x_.Np());

			south_temperature_.resize(grid_x_.Np());
			south_mass_fractions_.resize(grid_x_.Np());

			east_temperature_.resize(grid_y_.Np());
			east_mass_fractions_.resize(grid_y_.Np());

			if (hole_ == true)
			{
				west_temperature_.resize(grid_y_.Np());
				west_mass_fractions_.resize(grid_y_.Np());
			}
			
			for (int i = 0; i < grid_x_.Np(); i++)
			{
				north_mass_fractions_[i].resize(thermodynamicsMap_.NumberOfSpecies());
				south_mass_fractions_[i].resize(thermodynamicsMap_.NumberOfSpecies());
			}
			for (int i = 0; i < grid_y_.Np(); i++)
			{
				east_mass_fractions_[i].resize(thermodynamicsMap_.NumberOfSpecies());

				if (hole_ == true)
					west_mass_fractions_[i].resize(thermodynamicsMap_.NumberOfSpecies());
			}

			// Interpolated fields
			if (side_name == "North")
			{
				for (int i = 0; i < grid_x_.Np(); i++)
					Interpolate(grid_x_.x()[i], north_temperature_[i], north_mass_fractions_[i],
						coordinates_from_cfd, temperature_from_cfd, mass_fractions_from_cfd);
			}
			else if (side_name == "South")
			{
				for (int i = 0; i < grid_x_.Np(); i++)
					Interpolate(grid_x_.x()[i], south_temperature_[i], south_mass_fractions_[i],
						coordinates_from_cfd, temperature_from_cfd, mass_fractions_from_cfd);
			}
			else if (side_name == "East")
			{
				for (int i = 0; i < grid_y_.Np(); i++)
					Interpolate(grid_y_.x()[i], east_temperature_[i], east_mass_fractions_[i],
						coordinates_from_cfd, temperature_from_cfd, mass_fractions_from_cfd);
			}
			else if (side_name == "West")
			{
				for (int i = 0; i < grid_y_.Np(); i++)
					Interpolate(grid_y_.x()[i], west_temperature_[i], west_mass_fractions_[i],
						coordinates_from_cfd, temperature_from_cfd, mass_fractions_from_cfd);
			}
		}

		// Summary on the screen
		{
			{
				std::cout << "-------------------------------------------------------------" << std::endl;
				std::cout << "                 Gaseous phase from CFD                      " << std::endl;
				std::cout << "-------------------------------------------------------------" << std::endl;
				std::cout << " * North [mm]:           " << north_coordinate*1000. << std::endl;
				std::cout << " * South [mm]:           " << south_coordinate*1000. << std::endl;
				std::cout << " * East [mm]:            " << east_coordinate*1000. << std::endl;
				std::cout << " * West [mm]:            " << west_coordinate*1000. << std::endl;
				std::cout << " * Internal radius [mm]: " << ri_*1000. << std::endl;
				std::cout << " * External radius [mm]: " << re_*1000. << std::endl;
				std::cout << " * Height [mm]:          " << H_*1000. << std::endl;
				std::cout << "-------------------------------------------------------------" << std::endl;
			}
		}
	}

	void DiskFromCFD::WriteOnFile(const boost::filesystem::path output_folder)
	{
		// North side
		{
			boost::filesystem::path file_name = output_folder / "Disk.North.out";
			std::ofstream fNorth;
			PrepareOutputFile(file_name, fNorth);
			for (int i = 0; i < grid_x_.Np(); i++)
			{
				fNorth << std::setprecision(9) << std::setw(22) << grid_x_.x()[i] * 1000.;
				fNorth << std::setprecision(9) << std::setw(22) << north_temperature_[i];
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					fNorth << std::setprecision(9) << std::setw(22) << north_mass_fractions_[i][j];
				fNorth << std::endl;
			}
			fNorth.close();
		}

		// South side
		{
			boost::filesystem::path file_name = output_folder / "Disk.South.out";
			std::ofstream fSouth;
			PrepareOutputFile(file_name, fSouth);
			for (int i = 0; i < grid_x_.Np(); i++)
			{
				fSouth << std::setprecision(9) << std::setw(22) << grid_x_.x()[i] * 1000.;
				fSouth << std::setprecision(9) << std::setw(22) << south_temperature_[i];
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					fSouth << std::setprecision(9) << std::setw(22) << south_mass_fractions_[i][j];
				fSouth << std::endl;
			}
			fSouth.close();
		}

		// East side
		{
			boost::filesystem::path file_name = output_folder / "Disk.East.out";
			std::ofstream fEast;
			PrepareOutputFile(file_name, fEast);
			for (int i = 0; i < grid_y_.Np(); i++)
			{
				fEast << std::setprecision(9) << std::setw(22) << grid_y_.x()[i] * 1000.;
				fEast << std::setprecision(9) << std::setw(22) << east_temperature_[i];
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					fEast << std::setprecision(9) << std::setw(22) << east_mass_fractions_[i][j];
				fEast << std::endl;
			}
			fEast.close();
		}

		// West side
		if (hole_ == true)
		{
			boost::filesystem::path file_name = output_folder / "Disk.west.out";
			std::ofstream fWest;
			PrepareOutputFile(file_name, fWest);
			for (int i = 0; i < grid_y_.Np(); i++)
			{
				fWest << std::setprecision(9) << std::setw(22) << grid_y_.x()[i] * 1000.;
				fWest << std::setprecision(9) << std::setw(22) << west_temperature_[i];
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					fWest << std::setprecision(9) << std::setw(22) << west_mass_fractions_[i][j];
				fWest << std::endl;
			}
			fWest.close();
		}
	}

	void DiskFromCFD::PrepareOutputFile(const boost::filesystem::path file_name, std::ofstream& fOutput)
	{
		fOutput.open(file_name.c_str(), std::ios::out);
		fOutput.setf(std::ios::scientific);

		unsigned int count = 1;
		OpenSMOKE::PrintTagOnASCIILabel(22, fOutput, "x[mm]", count);
		OpenSMOKE::PrintTagOnASCIILabel(22, fOutput, "T[K]", count);
		for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
			OpenSMOKE::PrintTagOnASCIILabel(22, fOutput, thermodynamicsMap_.NamesOfSpecies()[j] + "_w", count);
		fOutput << std::endl;
	}

	void Interpolate(	const double x_req, double& temperature_req, std::vector<double>& mass_fractions_req,
						const std::vector<double>& coordinates_from_cfd,
						const std::vector<double>& temperature_from_cfd,
						const std::vector< std::vector<double> >& mass_fractions_from_cfd)
	{
		for (unsigned int j = 1; j < coordinates_from_cfd.size(); j++)
		{
			if (x_req <= coordinates_from_cfd[j])
			{
				const double coefficient = (x_req - coordinates_from_cfd[j - 1]) / (coordinates_from_cfd[j] - coordinates_from_cfd[j-1]);

				temperature_req = temperature_from_cfd[j - 1] + (temperature_from_cfd[j] - temperature_from_cfd[j-1])*coefficient;
				for (unsigned int k = 0; k < mass_fractions_req.size(); k++)
					mass_fractions_req[k] = mass_fractions_from_cfd[k][j - 1] + (mass_fractions_from_cfd[k][j] - mass_fractions_from_cfd[k][j-1])*coefficient;

				return;
			}

		//	OpenSMOKE::ErrorMessage("Interpolate", "Interpolation error");
		}
	}
}
