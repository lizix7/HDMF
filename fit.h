#ifndef FIT_H
#define FIT_H

#include "global.h"
#include "parser.h"

using namespace std;

bool BasicSort(Point* lhs, Point* rhs, size_t i);
bool SortHD(Point* lhs, Point* rhs);
bool UniqueHD(Point* lhs, Point* rhs);

		
class LeastSquare {
	vector<Point*> total_points;
	vector<size_t> positions;
	vector<size_t> point_coefficients;
	vector<size_t> cube_coefficients;
	vector<vector<size_t> > combinations;
	vector<vector<size_t> > index_combinations;
	public:
		LeastSquare(Parser &);
		~LeastSquare();
		
		map<int,vector<Point*> > PointMap;
		map<int,Equation> EquationMap;
		
		void SelfCheck(int &, vector<size_t>, vector<HyperCube> &, vector<Point*> &, map<int,vector<Point*> >);
		void Combination(vector<vector<size_t> >, size_t, vector<size_t>);
		void CheckNeighbour(vector<size_t>, vector<size_t>, vector<vector<size_t> >, vector<HyperCube> &, vector<Point*> &, vector<Point*>, map<int,vector<Point*> >);
		void CheckNeighbours(vector<size_t>, vector<size_t>, vector<HyperCube> &, vector<Point*> &, vector<Point*>, map<int,vector<Point*> >);
		void Update(vector<size_t>, vector<HyperCube>, vector<Point*> &, map<int,vector<Point*> >&, map<int,Equation> &);
		void Traverse(vector<size_t>, vector<size_t>, size_t, vector<size_t>, int &, vector<HyperCube> &, vector<Point*> &, vector<Point*>, map<int,vector<Point*> >&, map<int,Equation> &);
		void DeleteDuplicate(vector<Point*> &);
		bool CheckError(vector<Point*>, Equation);
		Equation Solve_Matrix(vector<vector<double> > &, vector<double> &);
		Equation Calculate(vector<Point*>);
};

#endif 
