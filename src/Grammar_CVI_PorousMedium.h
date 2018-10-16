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
|   Copyright(C) 2015, 2014, 2013  Alberto Cuoci                          |
|   Source-code or binary products cannot be resold or distributed        |
|   Non-commercial use only                                               |
|   Cannot modify source-code for any purpose (cannot create              |
|   derivative works)                                                     |
|                                                                         |
\*-----------------------------------------------------------------------*/

#ifndef OpenSMOKE_Grammar_CVI_PorousMedium_H
#define OpenSMOKE_Grammar_CVI_PorousMedium_H

#include "dictionary/OpenSMOKE_DictionaryManager.h"
#include "dictionary/OpenSMOKE_DictionaryGrammar.h"
#include "dictionary/OpenSMOKE_DictionaryKeyWord.h"

namespace CVI
{
	class Grammar_CVI_PorousMedium : public OpenSMOKE::OpenSMOKE_DictionaryGrammar
	{
	protected:

		virtual void DefineRules()
		{
			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@FiberRadius",
				OpenSMOKE::SINGLE_MEASURE,
				"Radius of the fiber",
				true));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@FiberDensity",
				OpenSMOKE::SINGLE_MEASURE,
				"Fiber density",
				true));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@InitialPorosity",
				OpenSMOKE::SINGLE_DOUBLE,
				"Initial porosity",
				true));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@PorousSubstrate",
				OpenSMOKE::SINGLE_STRING,
				"Porous substrate type: polynomial | random | random_hardcore | polynomial_onehalf | from_spheres_to_cylinders | deutschmann_correlation",
				true));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@MassDiffusionMultiplier",
				OpenSMOKE::SINGLE_DOUBLE,
				"Mass diffusion multiplier",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ThresholdPorosity",
				OpenSMOKE::SINGLE_DOUBLE,
				"Threshold porosity (default: 1e-2)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@SmoothingCoefficientPorosity",
				OpenSMOKE::SINGLE_DOUBLE,
				"Smoothing coefficient for porosity (default: 200)",
				false));
		}
	};
}

#endif /* OpenSMOKE_Grammar_CVI_PorousMedium_H */
