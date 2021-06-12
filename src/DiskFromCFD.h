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

#ifndef OpenSMOKE_DiskFromCFD_H
#define OpenSMOKE_DiskFromCFD_H

// OpenSMOKE++ Definitions
#include "OpenSMOKEpp"

// CHEMKIN maps
#include "maps/Maps_CHEMKIN"

// 1D grid
#include "utilities/grids/adaptive/Grid1D.h"

namespace CVI
{
	//!  A class to solve the reaction-diffusion equations in 1D 
	/*!
	This class provides the tools to solve the reaction-diffusion equations in 1D 
	*/

	class DiskFromCFD
	{
	public:

		/**
		*@brief Default constructor
		*@param thermodynamicsMap	reference to the thermodynamic map
		*@param kineticsMap			reference to the kinetic map
		*@param grid_x				reference to 1D grid along x
		*@param grid_y				reference to 1D grid along y
		*/
		DiskFromCFD(	OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap,
						OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap,
						OpenSMOKE::Grid1D& grid_x, 
						OpenSMOKE::Grid1D& grid_y);

		void ReadFromFile(const boost::filesystem::path disk_file_name);
		void WriteOnFile(const boost::filesystem::path output_folder);


		const std::vector<double>& north_temperature() const { return north_temperature_; }
		const std::vector<double>& south_temperature() const { return south_temperature_; }
		const std::vector<double>& east_temperature() const { return east_temperature_; }
		const std::vector<double>& west_temperature() const { return west_temperature_; }
		
		const std::vector< std::vector<double> >& north_mass_fractions() const { return north_mass_fractions_; }
		const std::vector< std::vector<double> >& south_mass_fractions() const { return south_mass_fractions_; }
		const std::vector< std::vector<double> >& east_mass_fractions() const { return east_mass_fractions_; }
		const std::vector< std::vector<double> >& west_mass_fractions() const { return west_mass_fractions_; }

	protected:

		// References
		OpenSMOKE::ThermodynamicsMap_CHEMKIN&			thermodynamicsMap_;			//!< reference to the thermodynamic map
		OpenSMOKE::KineticsMap_CHEMKIN&					kineticsMap_;				//!< reference to the kinetic map
		OpenSMOKE::Grid1D&								grid_x_;					//!< reference to the 1D grid along the x direction
		OpenSMOKE::Grid1D&								grid_y_;					//!< reference to the 1D grid along the y direction

		void PrepareOutputFile(const boost::filesystem::path file_name, std::ofstream& fOutput);

		bool	hole_;
		double	ri_;
		double	re_;
		double	H_;

		std::vector<double> north_temperature_;
		std::vector<double> south_temperature_;
		std::vector<double> east_temperature_;
		std::vector<double> west_temperature_;

		std::vector< std::vector<double> > north_mass_fractions_;
		std::vector< std::vector<double> > south_mass_fractions_;
		std::vector< std::vector<double> > east_mass_fractions_;
		std::vector< std::vector<double> > west_mass_fractions_;

		unsigned int radial_coordinate_;
		unsigned int axial_coordinate_;

	};
}

#include "DiskFromCFD.hpp"

#endif /* OpenSMOKE_DiskFromCFD_H */
