#include "FEM.h"
#include <vector>
#ifdef USE_TBB
#include <tbb/tbb.h>
#endif

#define epsilon 1e-10

FEM::FEM(int cellNum, solverType type, smoothType smoothSolver)
{
	nCells = cellNum;
	linearSolver = type;
	source = MatrixXd::Zero(nCells + 1, nCells + 1);
	nodalCoef = MatrixXd::Zero(nCells + 1, nCells + 1);
	iter = 0;
	this->smoothSolver = smoothSolver;
}

FEM::~FEM()
{
}

void FEM::matrix_assembly(FunType f)
{
	//second order quadrature formula
	for (int i = 1; i < nCells; i++)
	{
		for (int j = 1; j < nCells; j++)
			source(i, j) = f((double)i / nCells, (double)j / nCells) / nCells / nCells;
	}
}

MatrixXd Amul(MatrixXd& x)
{
	Index n = x.rows();
	MatrixXd temp = MatrixXd::Zero(n, n);
#ifdef USE_TBB
    tbb::parallel_for(static_cast<Index>(1), n - 1, [&](Index i)
#else
    for (int i = 1; i < n - 1; i++)
#endif
	{
		for (int j = 1; j < n - 1; j++)
			temp(i, j) = x(i, j) * 4 - (x(i - 1, j) + x(i + 1, j) + x(i, j + 1) + x(i, j - 1));
	}
#ifdef USE_TBB
    );
#endif
	return temp;
}

double getMatrixDot(const MatrixXd& a, const MatrixXd& b)
{
	double ans = 0;
	for (int i = 0; i < a.rows(); i++)
	{
		for (int j = 0; j < a.cols(); j++)
			ans += a(i, j) * b(i, j);
	}

	return ans;
}

int conjugate_gradient(const MatrixXd& b, MatrixXd& x)
{
	Index n = b.rows();
	MatrixXd r(n + 1, n + 1), p(n + 1, n + 1), Ap(n + 1, n + 1);
	r = b - Amul(x);
	p = r;
	double alpha = 0, res = 0, res_old = getMatrixDot(r, r);
	int iter = 0;
	while (iter < n)
	{
		Ap = Amul(p);
		alpha = res_old / getMatrixDot(p, Ap);
		x += alpha * p;
		r -= alpha * Ap;
		iter++;
		res = getMatrixDot(r, r);
		if (sqrt(res) < epsilon)	break;
		p = r + (res / res_old) * p;
		res_old = res;
	}
	return iter;
}

void FEM::solve()
{
	iter = 0;
	nodalCoef = MatrixXd::Zero(nCells + 1, nCells + 1);
	switch (linearSolver)
	{
	case cg:
	{
		iter = conjugate_gradient(source, nodalCoef);
		break;
	}
	case mgV:
	{
		int p = 1, m1 = 10, m2 = 10;
		MatrixXd coarse = MatrixXd::Zero(6, 6);
		std::vector<MatrixXd> b_list;
		b_list.push_back(source);
		int depth = 0, total_depth = 0;
		while (b_list[b_list.size() - 1].rows() > 9)
		{
			MatrixXd tmp;
			projection(b_list[b_list.size() - 1], tmp);
			b_list.push_back(tmp);
		}
		total_depth = b_list.size();
		while (coarse.rows() < nodalCoef.rows())
		{
			depth++;
			while ((Amul(coarse) - b_list[total_depth - depth]).norm() > epsilon)
			{
				cycle(m1, m2, p, coarse, b_list[total_depth - depth]);
				iter++;
			}
			MatrixXd tmp;
			interpolation(coarse, tmp);
			coarse = tmp;
		}
		while ((Amul(coarse) - source).norm() > epsilon)
		{
			cycle(m1, m2, p, coarse, source);
			iter++;
		}
		nodalCoef = coarse;
		break;
	}
	case mgW:
	{
		int p = 2, m1 = 10, m2 = 0;
		MatrixXd coarse = MatrixXd::Zero(6, 6);
		std::vector<MatrixXd> b_list;
		b_list.push_back(source);
		int depth = 0, total_depth = 0;
		while (b_list[b_list.size() - 1].rows() > 9)
		{
			MatrixXd tmp;
			projection(b_list[b_list.size() - 1], tmp);
			b_list.push_back(tmp);
		}
		total_depth = b_list.size();
		while (coarse.rows() < nodalCoef.rows())
		{
			depth++;
			while ((Amul(coarse) - b_list[total_depth-depth]).norm() > epsilon)
			{
				cycle(m1, m2, p, coarse, b_list[total_depth - depth]);
				iter++;
			}
			MatrixXd tmp;
			interpolation(coarse, tmp);
			coarse = tmp;
		}
		while ((Amul(coarse) - source).norm() > epsilon)
		{
			cycle(m1, m2, p, coarse, source);
			iter++;
		}
		nodalCoef = coarse;
		break;
	}
	}
}

void FEM::interpolation(const MatrixXd& b1, MatrixXd& b2)
{
	Index n = b1.rows();//b2.rows()=2*n-1;b1.rows()=n;
	b2 = MatrixXd::Zero(2 * n - 1, 2 * n - 1);
	for (Index i = 1; i < n - 1; i++)
	{
		for (Index j = 1; j < n - 1; j++)
		{
			b2(2 * i, 2 * j) += b1(i, j);
			b2(2 * i + 1, 2 * j) += b1(i, j) / 2;
			b2(2 * i - 1, 2 * j) += b1(i, j) / 2;
			b2(2 * i, 2 * j - 1) += b1(i, j) / 2;
			b2(2 * i, 2 * j + 1) += b1(i, j) / 2;
			b2(2 * i + 1, 2 * j + 1) += b1(i, j) / 2;
			b2(2 * i - 1, 2 * j - 1) += b1(i, j) / 2;
		}
	}
}

void FEM::projection(const MatrixXd& b1, MatrixXd& b2)
{
	Index n = b1.rows() / 2 + 1;//b1.rows()=2*n-1;b2.rows()=n;
	b2 = MatrixXd::Zero(n, n);
#ifdef USE_TBB
	tbb::parallel_for(static_cast<Index>(1), n - 1, [&](Index i)
#else
    for (Index i = 1; i < n - 1; i++)
#endif
	{
		for (Index j = 1; j < n - 1; j++)
			b2(i, j) = (b1(2 * i, 2 * j) + (b1(2 * i + 1, 2 * j) + b1(2 * i - 1, 2 * j) + b1(2 * i, 2 * j + 1) + b1(2 * i, 2 * j - 1) + b1(2 * i + 1, 2 * j + 1) + b1(2 * i - 1, 2 * j - 1)) / 2) / 4;
	}
#ifdef USE_TBB
    );
#endif
}

void FEM::smoothing(int m, const MatrixXd& b, MatrixXd& x)
{
	int n = (int)x.rows();
	switch (smoothSolver)
	{
	case gs:
	{
		for (int k = 0; k < m; k++)
			for (int i = 1; i < n - 1; i++)
				for (int j = 1; j < n - 1; j++)
					x(i, j) = (x(i - 1, j) + x(i + 1, j) + x(i, j + 1) + x(i, j - 1) + b(i, j)) / 4;
		return;
	}
	case richardson:
	{
		double gamma = 10.0;
		for (int k = 0; k < m; k++)
			x -= (Amul(x) - b) / gamma;
		return;
	}
	case jacobi:
	{
		MatrixXd tmp = x;
		for (int k = 0; k < m; k++)
		{
			tmp.swap(x);
#ifdef USE_TBB
    		tbb::parallel_for(1, n - 1, 1, [&](int i)
#else
    		for (int i = 1; i < n - 1; i++)
#endif
			{
				for (int j = 1; j < n - 1; j++)
					x(i, j) = (tmp(i - 1, j) + tmp(i + 1, j) + tmp(i, j + 1) + tmp(i, j - 1) + b(i, j)) / 4;
			}
#ifdef USE_TBB
    		);
#endif
		}
		return;
	}
	}
}

void FEM::cycle(int m1, int m2, int p, MatrixXd& x, const MatrixXd& b)
{
	Index n = b.rows();
	if (n < 10)
	{
		conjugate_gradient(b, x);
	}
	else
	{
		smoothing(m1, b, x);

		MatrixXd g = b - Amul(x);
		MatrixXd g_, q = MatrixXd::Zero(n / 2 + 1, n / 2 + 1);
		projection(g, g_);

		for (int i = 0; i < p; i++)
			cycle(m1, m2, p, q, g_);

		interpolation(q, g);
		x += g;

		smoothing(m2, b, x);
	}
}

double FEM::error(FunType u)
{
	double error2 = 0;
	int M = 2;
	MatrixXd x;
	interpolation(nodalCoef, x);
	int n = x.rows() - 1;
#ifdef USE_TBB
    tbb::parallel_for(1, n, 1, [&](int i)
#else
    for (int i = 1; i < n; i++)
#endif
	{
		for (int j = 1; j < n; j++)
		{
			error2 += pow(u((double)i / n, (double)j / n) - x(i, j), 2);
		}
	}
#ifdef USE_TBB
    );
#endif
	return sqrt(error2) / n;
}