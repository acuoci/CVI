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
|   Copyright(C) 2015, 2014, 2013  Alberto Cuoci                          |
|   Source-code or binary products cannot be resold or distributed        |
|   Non-commercial use only                                               |
|   Cannot modify source-code for any purpose (cannot create              |
|   derivative works)                                                     |
|                                                                         |
\*-----------------------------------------------------------------------*/

#ifndef OpenSMOKE_Grammar_CVI_Solver_H
#define OpenSMOKE_Grammar_CVI_Solver_H

#include "dictionary/OpenSMOKE_DictionaryManager.h"
#include "dictionary/OpenSMOKE_DictionaryGrammar.h"
#include "dictionary/OpenSMOKE_DictionaryKeyWord.h"

namespace CVI
{
	class Grammar_CVI_Solver : public OpenSMOKE::OpenSMOKE_DictionaryGrammar
	{
	protected:

		virtual void DefineRules()
		{
			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Type",
				OpenSMOKE::SINGLE_STRING,
				"Type of problem to be solved: Capillary | 1D | 2D",
				true));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@DetailedSurfaceChemistry",
				OpenSMOKE::SINGLE_BOOL,
				"Detailed vs Global heterogeneous chemistry",
				true));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@GasDaeSpecies",
				OpenSMOKE::SINGLE_STRING,
				"Name of gas species chosen for the algebraic equation closure",
				true));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@SurfaceDaeSpecies",
				OpenSMOKE::SINGLE_STRING,
				"Name of surface species chosen for the algebraic equation closure",
				true));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Symmetry",
				OpenSMOKE::SINGLE_STRING,
				"Type of symmetry: Planar | Cylindrical",
				true));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@KineticsFolder",
				OpenSMOKE::SINGLE_PATH,
				"Name of the folder containing the kinetic scheme (XML Version)",
				true,
				"@KineticsPreProcessor",
				"none",
				"none"));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Output",
				OpenSMOKE::SINGLE_PATH,
				"Name of the folder containing the output files",
				true ));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@KineticsPreProcessor",
				OpenSMOKE::SINGLE_DICTIONARY,
				"Name of the dictionary containing the list of kinetic files to be interpreted",
				true,
				"@KineticsFolder",
				"none",
				"none"));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@PlugFlowReactor",
				OpenSMOKE::SINGLE_DICTIONARY,
				"Name of the dictionary/dictionaries defining the plug flow reactor",
				true,
				"@DiskFromCFD",
				"none",
				"none" ));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@DiskFromCFD",
				OpenSMOKE::SINGLE_PATH,
				"Name of file containing the result of the CFD simulation",
				true,
				"@PlugFlowReactor",
				"none",
				"none"));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@PorousMedium",
				OpenSMOKE::SINGLE_DICTIONARY,
				"Name of the dictionary/dictionaries defining the porous medium",
				true));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@HeterogeneousMechanism",
				OpenSMOKE::SINGLE_DICTIONARY,
				"Name of the dictionary/dictionaries defining the heterogeneous mechanism",
				true));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@PorosityDefect",
				OpenSMOKE::SINGLE_DICTIONARY,
				"Name of the dictionary/dictionaries defining the porosity defect (if any)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@InletStream",
				OpenSMOKE::SINGLE_DICTIONARY,
				"Name of the dictionary/dictionaries defining the inlet gas composition, temperature and pressure",
				true));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@InitialConditions",
				OpenSMOKE::SINGLE_DICTIONARY,
				"Name of the dictionary/dictionaries defining the initial gas composition, temperature and pressure",
				true));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@XPoints",
				OpenSMOKE::SINGLE_INT,
				"Number of grid points along the x axis",
				true));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@YPoints",
				OpenSMOKE::SINGLE_INT,
				"Number of grid points along the y axis",
				true));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@XLength",
				OpenSMOKE::SINGLE_MEASURE,
				"Length of computational domain along the x axis",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@RadiusInternal",
				OpenSMOKE::SINGLE_MEASURE,
				"Internal radius (needed only in case of cylindrical symmetry)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@RadiusExternal",
				OpenSMOKE::SINGLE_MEASURE,
				"External radius (needed only in case of cylindrical symmetry)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@YLength",
				OpenSMOKE::SINGLE_MEASURE,
				"Length of computational domain along the y axis",
				true));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@XStretchingFactor",
				OpenSMOKE::SINGLE_DOUBLE,
				"Stretching factor along the x axis",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@YStretchingFactor",
				OpenSMOKE::SINGLE_DOUBLE,
				"Stretching factor along the y axis",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@DaeParameters",
				OpenSMOKE::SINGLE_DICTIONARY,
				"Dictionary containing the numerical parameters for solving the stiff DAE system",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@DaeTimeInterval",
				OpenSMOKE::SINGLE_MEASURE,
				"Interval of time for successive DAE system solutions",
				true));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@OdeParameters",
				OpenSMOKE::SINGLE_DICTIONARY,
				"Dictionary containing the numerical parameters for solving the stiff ODE system",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@OdeEndTime",
				OpenSMOKE::SINGLE_MEASURE,
				"Time for solving the surface species equations for determining the initial conditions (default: 1 s)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@OnTheFlyROPA",
				OpenSMOKE::SINGLE_DICTIONARY,
				"Dictionary specifying the details for carrying out the ROPA (on the fly)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@TecplotTimeInterval",
				OpenSMOKE::SINGLE_MEASURE,
				"Interval of time for writing Tecplot output",
				true));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@TimeTotal",
				OpenSMOKE::SINGLE_MEASURE,
				"Total time of simulation",
				true));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@StepsVideo",
				OpenSMOKE::SINGLE_INT,
				"Number steps to update info on the screen",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@StepsFile",
				OpenSMOKE::SINGLE_INT,
				"Number steps to update info on files",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@StepsPlugFlow",
				OpenSMOKE::SINGLE_INT,
				"Number steps to update the plug flow reactor",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ResidenceTime",
				OpenSMOKE::SINGLE_MEASURE,
				"Residence time to be simulated in the gaseous phase",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@CapillaryDiameter",
				OpenSMOKE::SINGLE_MEASURE,
				"Initial diameter of capillary",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@XVelocity",
				OpenSMOKE::SINGLE_MEASURE,
				"Uniform velocity along the x axis",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@YVelocity",
				OpenSMOKE::SINGLE_MEASURE,
				"Uniform velocity along the y axis",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@DerivativeMassFractions",
				OpenSMOKE::SINGLE_STRING,
				"Derivative of mass fractions: Backward | Forward | Centered (default: Centered)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@DerivativeBulkDensity",
				OpenSMOKE::SINGLE_STRING,
				"Derivative of bulk density: Backward | Forward | Centered (default: Centered)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@DerivativeEffectiveDiffusivity",
				OpenSMOKE::SINGLE_STRING,
				"Derivative of effective diffusivity: Backward | Forward | Centered (default: Centered)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Backup",
				OpenSMOKE::SINGLE_PATH,
				"Name of backup file (XML Version)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ReadjustBackup",
				OpenSMOKE::SINGLE_BOOL,
				"Adjust surface composition after reading the backup file (default: false)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@TemperatureProfile",
				OpenSMOKE::SINGLE_DICTIONARY,
				"Dictionary defining the temperature profile",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ImpermeableWalls",
				OpenSMOKE::VECTOR_STRING,
				"List of impermeable walls (default: none)",
				false));
		}
	};
}

#endif /* OpenSMOKE_Grammar_CVI_Solver_H */
