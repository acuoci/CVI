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

PlugFlowReactorProfiles::PlugFlowReactorProfiles(const Eigen::VectorXd& x)
{
	n_ = x.size();
	x_ = x;
	xstart_ = x_(0);
	xend_ = x_(n_ - 1);
}

void PlugFlowReactorProfiles::Interpolate(const double x_req, const std::vector<Eigen::VectorXd>& y, Eigen::VectorXd& y_req)
{
	if (n_ != y.size())
		OpenSMOKE::FatalErrorMessage("Interpolating fixed profile: the provided y vector is not consistent with the abscissas");

	if (x_req < xstart_)
		OpenSMOKE::FatalErrorMessage("Interpolating fixed profile: the requested coordinate is smaller than the minimum available coordinate");

	if (x_req > xend_)
		OpenSMOKE::FatalErrorMessage("Interpolating fixed profile: the requested coordinate is larger than the maximum available coordinate");

	unsigned int nv = y[0].size();

	y_req.resize(nv);
	y_req.setZero();
	for (unsigned int j = 1; j < n_; j++)
	{
		if (x_req <= x_(j))
		{
			for (unsigned int k = 0; k < nv; k++)
				y_req(k) = y[j - 1](k) + (y[j](k) - y[j - 1](k)) / (x_(j) - x_(j - 1)) * (x_req - x_(j - 1));

			return;
		}
	}
}