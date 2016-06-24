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

#include "math/multivalue-ode-solvers/MultiValueSolver"

class OpenSMOKE_PlugFlowReactorCoupled_OdeSystem
{
public:

	OpenSMOKE_PlugFlowReactorCoupled_OdeSystem() {};

	void assignReactor(CVI::PlugFlowReactorCoupled *reactor);

private:

	CVI::PlugFlowReactorCoupled *ptReactor;

protected:

	unsigned int ne_;

	void MemoryAllocation()
	{
	}

	virtual void Equations(const Eigen::VectorXd& y, const double t, Eigen::VectorXd& f)
	{
		ptReactor->Equations(t, y.data(), f.data());
	}

	void Jacobian(const Eigen::VectorXd &y, const double t, Eigen::MatrixXd &J)
	{
	};

	void Print(const double t, const Eigen::VectorXd &y)
	{
		ptReactor->Print(t, y.data());
	}
};

void OpenSMOKE_PlugFlowReactorCoupled_OdeSystem::assignReactor(CVI::PlugFlowReactorCoupled *reactor)
{
	ptReactor = reactor;
}


bool PlugFlowReactorCoupled_Ode(CVI::PlugFlowReactorCoupled* reactor, const double t0, const double tf)
{
	std::cout << "Plug-flow reactor solution (OpenSMOKE++)..." << std::endl;

	const unsigned int neq = reactor->NumberOfEquations();

	// Initial conditions
	Eigen::VectorXd yInitial(neq);
	reactor->UnknownsVector(yInitial.data());

	typedef OdeSMOKE::KernelDense<OpenSMOKE_PlugFlowReactorCoupled_OdeSystem> denseOde;
	typedef OdeSMOKE::MethodGear<denseOde> methodGear;
	OdeSMOKE::MultiValueSolver<methodGear> ode_solver;

	// Set initial conditions
	ode_solver.assignReactor(reactor);
	ode_solver.SetInitialConditions(t0, yInitial);

	// Minimum/Maximum Constraints
	{
		Eigen::VectorXd xMin(neq);
		xMin.setZero();
		ode_solver.SetMinimumValues(xMin);

		Eigen::VectorXd xMax(neq);
		xMax.setConstant(1.);
		ode_solver.SetMaximumValues(xMax);
	}

	// Verbose
	ode_solver.SetPrint(true);

	// Solve the system
	double timeStart = OpenSMOKE::OpenSMOKEGetCpuTime();
	OdeSMOKE::OdeStatus status = ode_solver.Solve(tf);
	double timeEnd = OpenSMOKE::OpenSMOKEGetCpuTime();

	// Check the solution
	if (status > 0)
	{
		if (reactor->verbose_output() == true)
		{
			std::string message("OpenSMOKE++ Ode System successfully solved: ");

			if (status == 1)	message += "ODE_STATUS_TO_BE_INITIALIZED";
			else if (status == 2)	message += "ODE_STATUS_CONTINUATION";
			else if (status == 3)	message += "ODE_STATUS_STOP_INTEGRATION_FOR_SMALL_YPRIME_NORM1";

			std::cout << message << std::endl;

			std::cout << std::endl;
			std::cout << " * CPU time (s):                   " << timeEnd - timeStart << std::endl;
			std::cout << " * number of steps:                " << ode_solver.numberOfSteps() << std::endl;
			std::cout << " * number of functions:            " << ode_solver.numberOfFunctionCalls() << std::endl;
			std::cout << " * number of solutions:            " << ode_solver.numberOfLinearSystemSolutions() << std::endl;
			std::cout << " * number of Jacobians:            " << ode_solver.numberOfJacobianEvaluations() << std::endl;
			std::cout << " * number of factorizations:       " << ode_solver.numberOfMatrixFactorizations() << std::endl;
			std::cout << " * number of functions (Jacobian): " << ode_solver.numberOfFunctionCallsForJacobian() << std::endl;
			std::cout << " * last order:                     " << ode_solver.lastOrderUsed() << std::endl;
			std::cout << " * last step size:                 " << std::scientific << ode_solver.lastStepUsed() << std::endl;
			std::cout << std::endl;
		}

		Eigen::VectorXd ysol(neq);
		ode_solver.Solution(ysol);
		reactor->CorrectedUnknownsVector(ysol.data());

		return true;
	}
	else
	{
		std::string message("OpenSMOKE++ Ode Solver Error: ");
		if (status == -2)	message += "ODE_STATUS_MAX_NUMBER_OF_STEPS_REACHED";
		else if (status == -2)	message += "ODE_STATUS_TOO_STRICT_TOLERANCES";
		else if (status == -3)	message += "ODE_STATUS_ILLEGAL_CONTINUATION_REQUEST";
		else if (status == -4)	message += "ODE_STATUS_MAX_NUMBER_ERRORTEST_FAILURES";
		else if (status == -5)	message += "ODE_STATUS_MAX_NUMBER_CONVERGENCETEST_FAILURES";
		else if (status == -6)	message += "ODE_STATUS_TOO_SMALL_STEP_SIZE";
		else if (status == -7)	message += "ODE_STATUS_ILLEGAL_MAX_INDEPENDENT_VARIABLE";
		else if (status == -8)	message += "ODE_STATUS_ILLEGAL_CONSTRAINTS";
		else if (status == -10)	message += "ODE_STATUS_EXCEPTION_HANDLING_STOPP";
		else if (status == -12)	message += "ODE_STATUS_YOU_CANNOT_OVERSHOOT_TCRITIC";

		std::cout << message << std::endl;

		reactor->CorrectedUnknownsVector(yInitial.data());

		return false;
	}
}
