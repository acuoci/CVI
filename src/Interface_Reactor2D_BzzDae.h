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

class OpenSMOKE_Reactor2D_BzzDaeSystem : public BzzDaeSystemObject
{
public:

	void assign(CVI::Reactor2D *reactor);

	CVI::Reactor2D *ptReactor;
	virtual void GetSystemFunctions(BzzVector &y, double t, BzzVector &dy);
	virtual void ObjectBzzPrint(void);
};

void OpenSMOKE_Reactor2D_BzzDaeSystem::ObjectBzzPrint(void)
{
}

void OpenSMOKE_Reactor2D_BzzDaeSystem::GetSystemFunctions(BzzVector &x, double t, BzzVector &f)
{
	double* ptx = x.GetHandle();
	double* ptf = f.GetHandle();

	ptReactor->Equations(t, ptx, ptf);
}

void OpenSMOKE_Reactor2D_BzzDaeSystem::assign(CVI::Reactor2D *reactor)
{
	ptReactor = reactor;
}

void DaePrint(BzzVector &y, double t)
{
	reactor2d->Print(t, y.GetHandle());
}

#include "math\multivalue-dae-solvers\interfaces\Band_BzzDae.h"
