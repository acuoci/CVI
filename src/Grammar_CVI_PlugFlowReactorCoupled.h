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

#ifndef OpenSMOKE_Grammar_CVI_PlugFlowReactorCoupled_H
#define OpenSMOKE_Grammar_CVI_PlugFlowReactorCoupled_H

#include "dictionary/OpenSMOKE_DictionaryManager.h"
#include "dictionary/OpenSMOKE_DictionaryGrammar.h"
#include "dictionary/OpenSMOKE_DictionaryKeyWord.h"

namespace CVI
{
	class Grammar_CVI_PlugFlowReactorCoupled : public OpenSMOKE::OpenSMOKE_DictionaryGrammar
	{
	protected:

		virtual void DefineRules()
		{
			// Mandatory

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Velocity",
				OpenSMOKE::SINGLE_MEASURE,
				"Inlet velocity of plug flow reactor",
				true));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@HydraulicDiameter",
				OpenSMOKE::SINGLE_MEASURE,
				"Hydraulic diameter of the plug flow",
				true));

			// Optional

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@InletLength",
				OpenSMOKE::SINGLE_MEASURE,
				"Length of the inlet section",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ChannelWidth",
				OpenSMOKE::SINGLE_MEASURE,
				"Width of the channel",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@AsymptoticNusselt",
				OpenSMOKE::SINGLE_DOUBLE,
				"Asymptotic Nusselt number",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@InternalBoundaryLayer",
				OpenSMOKE::SINGLE_BOOL,
				"Internal boundary layer limitations",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@GeometricPattern",
				OpenSMOKE::SINGLE_STRING,
				"Geometric pattern: OneSide | ThreeSides",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Coupling",
				OpenSMOKE::SINGLE_BOOL,
				"Coupling with the carbon felt",
				true));
		}
	};
}

#endif /* OpenSMOKE_Grammar_CVI_PlugFlowReactorCoupled_H */
