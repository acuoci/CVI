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

#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/foreach.hpp>

namespace CVI
{
	void Interpolate(const double x_req, double& temperature_req, double& pressure_req, 
		std::vector<double>& mass_fractions_req,
		const std::vector<double>& coordinates_from_cfd,
		const std::vector<double>& temperature_from_cfd,
		const std::vector<double>& pressure_from_cfd,
		const std::vector< std::vector<double> >& mass_fractions_from_cfd);

	DiskFromCFD::DiskFromCFD(	OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap,
								OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap,
								OpenSMOKE::Grid1D& grid_x,
								OpenSMOKE::Grid1D& grid_y ) :
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

	void DiskFromCFD::ReadFromFile(const boost::filesystem::path disk_file_name, const double pressure_from_input_file)
	{
		std::cout << "Reading disk from file..." << std::endl;

		bool edge_centered_policy = false;

		boost::property_tree::ptree ptree;
    		boost::property_tree::read_xml( (disk_file_name).string(), ptree );

		std::cout << " * Check axial coordinate..." << std::endl;
		{
			boost::optional< boost::property_tree::ptree& > child = ptree.get_child_optional("opensmoke.AxialCoordinate");

			if (child)
			{
				std::string coordinate = ptree.get<std::string>("opensmoke.AxialCoordinate");
				coordinate.erase(std::remove(coordinate.begin(), coordinate.end(), '\n'), coordinate.end());
				if (coordinate == "x")			axial_coordinate_ = 0;
				else if (coordinate == "y")		axial_coordinate_ = 1;
				else if (coordinate == "z")		axial_coordinate_ = 2;
				else OpenSMOKE::ErrorMessage("DiskFromCFD::ReadFromFile", "Wrong AxialCoordinate: x | y (default) | z");
			}
		}

		std::cout << " * Check radial coordinate..." << std::endl;
		{
			boost::optional< boost::property_tree::ptree& > child = ptree.get_child_optional("opensmoke.RadialCoordinate");

			if (child)
			{
				std::string coordinate = ptree.get<std::string>("opensmoke.RadialCoordinate");
				coordinate.erase(std::remove(coordinate.begin(), coordinate.end(), '\n'), coordinate.end());
				if (coordinate == "x")			radial_coordinate_ = 0;
				else if (coordinate == "y")		radial_coordinate_ = 1;
				else if (coordinate == "z")		radial_coordinate_ = 2;
				else OpenSMOKE::ErrorMessage("DiskFromCFD::ReadFromFile", "Wrong RadialCoordinate: x (default) | y | z");
			}
		}

		std::cout << " * Check policy..." << std::endl;
		{
			boost::optional< boost::property_tree::ptree& > child = ptree.get_child_optional("opensmoke.Policy");
			if (child)
			{
				std::string policy = ptree.get<std::string>("opensmoke.Policy");
				policy.erase(std::remove(policy.begin(), policy.end(), '\n'), policy.end());
				if (policy == "edge-centered")		edge_centered_policy = true;
				else if (policy == "face-centered")	edge_centered_policy = false;
				else OpenSMOKE::ErrorMessage("DiskFromCFD::ReadFromFile", "Wrong policy: face-centered (default) | edge-centered)");
			}
		}

		try
		{
			std::cout << " * Check kinetic mechanism..." << std::endl;
			const unsigned int ns = ptree.get<unsigned int>("opensmoke.NumberOfSpecies");
			
			if (ns != thermodynamicsMap_.NumberOfSpecies())
				OpenSMOKE::ErrorMessage("DiskFromCFD::ReadFromFile", "The kinetic mechanism adopted in the exported Disk data does not match with the current kinetic mechanism");

			std::vector<std::string> names_species(ns);

			std::stringstream names_of_species_string;
			names_of_species_string.str( ptree.get< std::string >("opensmoke.NamesOfSpecies") );
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

		BOOST_FOREACH( boost::property_tree::ptree::value_type const& node, ptree.get_child( "opensmoke" ) ) 
		{
			if (node.first == "Data")
			{

			boost::property_tree::ptree subtree = node.second; 
			const std::string side = subtree.get<std::string>("<xmlattr>.side");

			std::cout << " * Importing disk data from side " << side << std::endl;

			// Number of points
			const unsigned int np = subtree.get<unsigned int>("NumberOfPoints");
			std::vector<double> coordinates_from_cfd(np);
			std::vector<double> temperature_from_cfd(np);
			std::vector<double> pressure_from_cfd(np);
			std::vector< std::vector<double> > mass_fractions_from_cfd(thermodynamicsMap_.NumberOfSpecies());
			for (unsigned int i = 0; i < thermodynamicsMap_.NumberOfSpecies(); i++)
				mass_fractions_from_cfd[i].resize(np);

			// Data from CFD
			std::vector<size_t> indices_increasing(np);
		
			// Point coordinates
			{
				std::vector<Eigen::VectorXd> coordinates(3);
				for (unsigned int i = 0; i < 3; i++)
					coordinates[i].resize(np);

				std::stringstream points_xml;
				points_xml.str( subtree.get< std::string >("Points") );
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

				if (side == "North")
				{
					if (relevant_coordinate != radial_coordinate_)
						OpenSMOKE::ErrorMessage("DiskFromCFD::ReadFromFile", "North side is not aligned with radial coordinate specified in XML file");

					north_coordinate = coordinates[axial_coordinate_](0);
					std::cout << " * North coordinate (mm): " << north_coordinate*1000. << std::endl;
				}
				else if (side == "South")
				{
					if (relevant_coordinate != radial_coordinate_)
						OpenSMOKE::ErrorMessage("DiskFromCFD::ReadFromFile", "South side is not aligned with radial coordinate specified in XML file");

					south_coordinate = coordinates[axial_coordinate_](0);
					std::cout << " * South coordinate (mm): " << south_coordinate * 1000. << std::endl;
				}
				else if (side == "East")
				{
					if (relevant_coordinate != axial_coordinate_)
						OpenSMOKE::ErrorMessage("DiskFromCFD::ReadFromFile", "East side is not aligned with axial coordinate specified in XML file");

					east_coordinate = coordinates[radial_coordinate_](0);
					std::cout << " * East coordinate (mm): " << east_coordinate * 1000. << std::endl;
				}
				else if (side == "West")
				{
					if (relevant_coordinate != axial_coordinate_)
						OpenSMOKE::ErrorMessage("DiskFromCFD::ReadFromFile", "West side is not aligned with axial coordinate specified in XML file");

					west_coordinate = coordinates[radial_coordinate_](0);
					std::cout << " * West coordinate (mm): " << west_coordinate * 1000. << std::endl;
				}
			}

			// Temperature
			{
				boost::optional< boost::property_tree::ptree& > child = subtree.get_child_optional("Temperature");
				if (!child)
					OpenSMOKE::ErrorMessage("DiskFromCFD::ReadFromFile", "Missing Temperature data in XML file");

				std::stringstream temperature_xml;
				temperature_xml.str( subtree.get< std::string >("Temperature") );

				std::vector<double> temperature_from_cfd_unsorted(np);
				for (unsigned int i = 0; i < np; i++)
					temperature_xml >> temperature_from_cfd_unsorted[i];

				for (unsigned int i = 0; i < np; i++)
					temperature_from_cfd[i] = temperature_from_cfd_unsorted[indices_increasing[i]];
			}

			// Pressure
			{
				boost::optional< boost::property_tree::ptree& > child = subtree.get_child_optional("Pressure");
				if (!child)
				{
					std::cout << "DiskFromCFD::ReadFromFile: pressure from input file = " << pressure_from_input_file << " Pa" << std::endl;

					for (unsigned int i = 0; i < np; i++)
						pressure_from_cfd[i] = pressure_from_input_file;
				}
				else
				{
					std::stringstream pressure_xml;
					pressure_xml.str(subtree.get< std::string >("Pressure"));

					std::vector<double> pressure_from_cfd_unsorted(np);
					for (unsigned int i = 0; i < np; i++)
						pressure_xml >> pressure_from_cfd_unsorted[i];

					for (unsigned int i = 0; i < np; i++)
						pressure_from_cfd[i] = pressure_from_cfd_unsorted[indices_increasing[i]];
				}
			}

			// Mass fractions
			{
				boost::optional< boost::property_tree::ptree& > child = subtree.get_child_optional("MassFractions");
				if (!child)
					OpenSMOKE::ErrorMessage("DiskFromCFD::ReadFromFile", "Missing MassFractions data in XML file");

				std::stringstream mass_fractions_xml;
				mass_fractions_xml.str( subtree.get< std::string >("MassFractions") );

				std::vector< std::vector<double> > mass_fractions_from_cfd_unsorted(thermodynamicsMap_.NumberOfSpecies());
				for (unsigned int i = 0; i < thermodynamicsMap_.NumberOfSpecies(); i++)
					mass_fractions_from_cfd_unsorted[i].resize(np);

				for (unsigned int i = 0; i < np; i++)
					for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
						mass_fractions_xml >> mass_fractions_from_cfd_unsorted[j][i];

				for (unsigned int i = 0; i < np; i++)
					for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
						mass_fractions_from_cfd[j][i] = mass_fractions_from_cfd_unsorted[j][indices_increasing[i]];
			}

			// Extend data
			{

				if (side == "North" || side == "South")
				{
					std::cout << "Extending" << std::endl;

					if (coordinates_from_cfd[0] < ri_ || coordinates_from_cfd[np - 1] > re_)
					{
						std::cout << "Provided radii (mm): " << ri_*1000. << " " << re_*1000. << std::endl;
						std::cout << "Coordinates from CFD (mm): " << coordinates_from_cfd[0]*1000. << " " << coordinates_from_cfd[np - 1]*1000. << std::endl;
						OpenSMOKE::ErrorMessage("DiskFromCFD::ReadFromFile", "The size of the disk from CFD does not match the input size");
					}

					// Initial value
					coordinates_from_cfd.insert(coordinates_from_cfd.begin(), ri_*(0.9999999));
					temperature_from_cfd.insert(temperature_from_cfd.begin(), temperature_from_cfd[0]);
					pressure_from_cfd.insert(pressure_from_cfd.begin(), pressure_from_cfd[0]);
					for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
						mass_fractions_from_cfd[j].insert(mass_fractions_from_cfd[j].begin(), mass_fractions_from_cfd[j][0]);

					// Final value
					coordinates_from_cfd.push_back(re_*(1.0000001));
					temperature_from_cfd.push_back(temperature_from_cfd[np]);
					pressure_from_cfd.push_back(pressure_from_cfd[np]);
					for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
						mass_fractions_from_cfd[j].push_back(mass_fractions_from_cfd[j][np]);
				}
			}

			if (side == "East" || side == "West")
			{
				for (unsigned int i = 0; i < np; i++)
					coordinates_from_cfd[i] -= south_coordinate;

				if (coordinates_from_cfd[0] < 0. || coordinates_from_cfd[np - 1] > H_)
				{
					std::cout << "Provided heights (mm): " << 0. << " " << H_*1000. << std::endl;
					std::cout << "Coordinates from CFD (mm): " << coordinates_from_cfd[0]*1000. << " " << coordinates_from_cfd[np - 1]*1000. << std::endl;
					OpenSMOKE::ErrorMessage("DiskFromCFD::ReadFromFile", "The size of the disk from CFD does not match the input size");
				}

				// Initial value
				coordinates_from_cfd.insert(coordinates_from_cfd.begin(), -0.0000001);
				temperature_from_cfd.insert(temperature_from_cfd.begin(), temperature_from_cfd[0]);
				pressure_from_cfd.insert(pressure_from_cfd.begin(), pressure_from_cfd[0]);
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					mass_fractions_from_cfd[j].insert(mass_fractions_from_cfd[j].begin(), mass_fractions_from_cfd[j][0]);

				// Final value
				coordinates_from_cfd.push_back(H_*(1.0000001));
				temperature_from_cfd.push_back(temperature_from_cfd[np]);
				pressure_from_cfd.push_back(pressure_from_cfd[np]);
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					mass_fractions_from_cfd[j].push_back(mass_fractions_from_cfd[j][np]);
			}

			north_temperature_.resize(grid_x_.Np());
			north_pressure_.resize(grid_x_.Np());
			north_mass_fractions_.resize(grid_x_.Np());

			south_temperature_.resize(grid_x_.Np());
			south_pressure_.resize(grid_x_.Np());
			south_mass_fractions_.resize(grid_x_.Np());

			east_temperature_.resize(grid_y_.Np());
			east_pressure_.resize(grid_y_.Np());
			east_mass_fractions_.resize(grid_y_.Np());

			if (hole_ == true)
			{
				west_temperature_.resize(grid_y_.Np());
				west_pressure_.resize(grid_y_.Np());
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
			if (side == "North")
			{
				for (int i = 0; i < grid_x_.Np(); i++)
					Interpolate(grid_x_.x()[i], north_temperature_[i], north_pressure_[i], north_mass_fractions_[i],
						coordinates_from_cfd, temperature_from_cfd, pressure_from_cfd, mass_fractions_from_cfd);
			}
			else if (side == "South")
			{
				for (int i = 0; i < grid_x_.Np(); i++)
					Interpolate(grid_x_.x()[i], south_temperature_[i], south_pressure_[i], south_mass_fractions_[i],
						coordinates_from_cfd, temperature_from_cfd, pressure_from_cfd, mass_fractions_from_cfd);
			}
			else if (side == "East")
			{
				for (int i = 0; i < grid_y_.Np(); i++)
					Interpolate(grid_y_.x()[i], east_temperature_[i], east_pressure_[i], east_mass_fractions_[i],
						coordinates_from_cfd, temperature_from_cfd, pressure_from_cfd, mass_fractions_from_cfd);
			}
			else if (side == "West")
			{
				for (int i = 0; i < grid_y_.Np(); i++)
					Interpolate(grid_y_.x()[i], west_temperature_[i], west_pressure_[i], west_mass_fractions_[i],
						coordinates_from_cfd, temperature_from_cfd, pressure_from_cfd, mass_fractions_from_cfd);
			}

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

	void Interpolate(	const double x_req, double& temperature_req, double& pressure_req, 
						std::vector<double>& mass_fractions_req,
						const std::vector<double>& coordinates_from_cfd,
						const std::vector<double>& temperature_from_cfd, const std::vector<double>& pressure_from_cfd,
						const std::vector< std::vector<double> >& mass_fractions_from_cfd)
	{
		for (unsigned int j = 1; j < coordinates_from_cfd.size(); j++)
		{
			if (x_req <= coordinates_from_cfd[j])
			{
				const double coefficient = (x_req - coordinates_from_cfd[j - 1]) / (coordinates_from_cfd[j] - coordinates_from_cfd[j-1]);

				temperature_req = temperature_from_cfd[j - 1] + (temperature_from_cfd[j] - temperature_from_cfd[j-1])*coefficient;
				pressure_req = pressure_from_cfd[j - 1] + (pressure_from_cfd[j] - pressure_from_cfd[j - 1]) * coefficient;
				for (unsigned int k = 0; k < mass_fractions_req.size(); k++)
					mass_fractions_req[k] = mass_fractions_from_cfd[k][j - 1] + (mass_fractions_from_cfd[k][j] - mass_fractions_from_cfd[k][j-1])*coefficient;

				return;
			}

		//	OpenSMOKE::ErrorMessage("Interpolate", "Interpolation error");
		}
	}
}
