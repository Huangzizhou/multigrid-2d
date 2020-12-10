#pragma once
#include <Eigen/Eigen>

typedef double (*FunType)(double, double);
using namespace Eigen;

enum solverType { cg, mgV, mgW };

enum smoothType {gs, richardson, jacobi};

class FEM
{
public:
	FEM(int cellNum, solverType type, smoothType smoothSolver);
	~FEM();

public:
	void matrix_assembly(FunType f);
	void solve();
	double error(FunType u);

	void interpolation(const MatrixXd& b1, MatrixXd& b2);
	void projection(const MatrixXd& b1, MatrixXd& b2);
	void smoothing(int m, const MatrixXd& b, MatrixXd& x);
	void cycle(int m1, int m2, int p, MatrixXd& x, const MatrixXd& b);

public:
	int					nCells;
	int					iter;
	MatrixXd			source;
	MatrixXd			nodalCoef;
	solverType			linearSolver;
	smoothType			smoothSolver;
};

