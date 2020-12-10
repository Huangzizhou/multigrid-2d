#include <iostream>
#include "FEM.h"
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <sstream>
#define Pi 3.14159265358979
using namespace std;


double f(double x, double y)
{
	return sin(Pi * y) * (x * (x - 1) * Pi * Pi - 2);
}

double u(double x, double y)
{
	return (x - 1) * x * sin(Pi * y);
}

int main()
{
	ofstream outFile;
	outFile.open("./data.csv", ios::out);
	outFile << "N,linearSolver,smoothSolver,iter,time(ms),error" << endl;
	string linearSolver[] = { "Conjugate Gradient","V Cycle","W Cycle" };
	string smoothSolver[] = { "GS","Richardson","Jacobi" };
	for (int t = 0; t < 5; t++)
	{
		int N = 10 * pow(2, t);
		FEM solver(N, cg, richardson);
		solver.matrix_assembly(&f);
		for (solverType type1 = cg; type1 <= mgW; type1 = (solverType)(type1 + 1))
		{
			solver.linearSolver = type1;
			for (smoothType type2 = gs; type2 <= jacobi; type2 = (smoothType)(type2 + 1))
			{
				solver.smoothSolver = type2;
				double start = clock();
				solver.solve();
				if (type1 == cg)
				{
					outFile << N << "," << linearSolver[type1] << "," << " " << "," << solver.iter << "," << clock() - start << "," << solver.error(&u) << endl;
					cout << N << "," << linearSolver[type1] << "," << " " << "," << solver.iter << "," << clock() - start << "," << solver.error(&u) << endl;
					break;
				}
				else
				{
					outFile << N << "," << linearSolver[type1] << "," << smoothSolver[type2] << "," << solver.iter << "," << clock() - start << "," << solver.error(&u) << endl;
					cout << N << "," << linearSolver[type1] << "," << smoothSolver[type2] << "," << solver.iter << "," << clock() - start << "," << solver.error(&u) << endl;
				}
			}
		}
		cout << "N=" << N << " over!" << endl;
	}
	outFile.close();
}

