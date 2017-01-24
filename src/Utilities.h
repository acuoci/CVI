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
|   This file is part of OpenSMOKE++ framework.                           |
|                                                                         |
|	License                                                               |
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

#ifndef OpenSMOKE_CVI_Utilities_H
#define OpenSMOKE_CVI_Utilities_H

/**
*@brief Reads a solution from a backup file
*@param t current time [s]
*@param path_file path to the backup file
*@param x axial coordinate [m]
*@param y axial coordinate [m]
*@param T temperature field [K]
*@param P pressure field [Pa]
*@param epsilon porosity field [-]
*@param omega_gas gas species mass fractions field [-]
*@param Z surface species field [-]
*@param Gamma molar concentration of surface sites [kmol/m2]
*@param gas_names_species names of gaseous species
*@param surface_names_species names of surface species
*/
void ReadFromBackupFile(const boost::filesystem::path path_file,
						double& t,
						Eigen::VectorXd& x, Eigen::VectorXd& y,
						Eigen::VectorXd& T, Eigen::VectorXd& P, Eigen::VectorXd& epsilon,
						std::vector< Eigen::VectorXd >& omega,
						std::vector< Eigen::VectorXd >& Z,
						std::vector< Eigen::VectorXd >& GammaFromEqs,
						std::vector<std::string>& gas_names_species,
						std::vector<std::string>& surface_names_species);

/**
*@brief Reads a solution from a backup file
*@param t current time [s]
*@param path_file path to the backup file
*@param thermodynamicsSurfaceMap thermodynamic map from which names of species can be extracted
*@param x axial coordinate [m]
*@param y axial coordinate [m]
*@param T temperature field [K]
*@param P pressure field [Pa]
*@param epsilon porosity field [-]
*@param omega_gas gas species mass fractions field [-]
*@param Z surface species field [-]
*@param Gamma molar concentration of surface sites [kmol/m2]
*/
void ReadFromBackupFile(const boost::filesystem::path path_file,
						double& t,
						OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN<double>& thermodynamicsSurfaceMap,
						const Eigen::VectorXd& x, const Eigen::VectorXd& y,
						Eigen::VectorXd& T, Eigen::VectorXd& P, Eigen::VectorXd& epsilon,
						std::vector< Eigen::VectorXd >& omega,
						std::vector< Eigen::VectorXd >& Z,
						std::vector< Eigen::VectorXd >& GammaFromEqs);

#include "Utilities.hpp"

#endif