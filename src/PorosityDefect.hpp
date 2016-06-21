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

namespace CVI
{
	PorosityDefect::PorosityDefect()
	{
		is_active_ = false;
		type_ = CIRCULAR;
		x_ = 0.;
		y_ = 0.;
		r_ = 0.;
		epsilon_ = 0.;
	}

	void PorosityDefect::ReadFromDictionary(OpenSMOKE::OpenSMOKE_Dictionary& dictionary)
	{
		if (dictionary.CheckOption("@X") == true)
		{
			std::string units;
			dictionary.ReadMeasure("@X", x_, units);
			if (units == "m")	x_ *= 1.;
			else if (units == "cm")	x_ *= 1.e-2;
			else if (units == "mm")	x_ *= 1.e-3;
			else if (units == "micron")	x_ *= 1.e-6;
			else OpenSMOKE::FatalErrorMessage("@X: Unknown units. Available units: m | cm | mm | micron");
		}

		if (dictionary.CheckOption("@Y") == true)
		{
			std::string units;
			dictionary.ReadMeasure("@Y", y_, units);
			if (units == "m")	y_ *= 1.;
			else if (units == "cm")	y_ *= 1.e-2;
			else if (units == "mm")	y_ *= 1.e-3;
			else if (units == "micron")	y_ *= 1.e-6;
			else OpenSMOKE::FatalErrorMessage("@Y: Unknown units. Available units: m | cm | mm | micron");
		}

		if (dictionary.CheckOption("@Radius") == true)
		{
			std::string units;
			dictionary.ReadMeasure("@Radius", r_, units);
			if (units == "m")	r_ *= 1.;
			else if (units == "cm")	r_ *= 1.e-2;
			else if (units == "mm")	r_ *= 1.e-3;
			else if (units == "micron")	r_ *= 1.e-6;
			else OpenSMOKE::FatalErrorMessage("@Radius: Unknown units. Available units: m | cm | mm | micron");
		}

		if (dictionary.CheckOption("@Porosity") == true)
			dictionary.ReadDouble("@Porosity", epsilon_);

		if (dictionary.CheckOption("@Type") == true)
		{
			std::string type;
			dictionary.ReadString("@Type", type);
			if (type == "circular") type_ = CIRCULAR;
			else  OpenSMOKE::FatalErrorMessage("Unknown @Type. Available types: circular");
		}

		is_active_ = true;
	}

	double PorosityDefect::set_porosity(const double x, const double y, const double epsilon)
	{
		if (type_ == CIRCULAR)
		{
			const double radius = std::sqrt(boost::math::pow<2>(x_ - x) + boost::math::pow<2>(y_ - y));

			if (radius <= r_)
				return epsilon_;
			else
				return epsilon;
		}
		else
		{
			return epsilon;
		}
	}
}
