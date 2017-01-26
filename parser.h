#ifndef PARSER_H
#define PARSER_H

#include "global.h"

using namespace std; 

class Parser {
	HyperCube cube;
	vector<vector<size_t> > combinations;
	vector<vector<size_t> > index_combinations;
	public:
		Parser(char*,char*);
		~Parser();
		
		vector<size_t> sizes;
		vector<size_t> reverse_sizes;
		vector<size_t> point_coefficients;
		vector<size_t> cube_coefficients;
		vector<Point> points;
		vector<Point> points_total;
		vector<HyperCube> cubes;
		
		void Parse(string);
		void ParseTotal(string);
		void MeshArea(size_t, int, int);
		void Analyse(vector<Point>);
		void Combination(vector<vector<size_t> >, size_t, vector<size_t>);
		void Allocate(vector<size_t>, size_t, vector<size_t>);
		void DeleteDuplicate(vector<double> &) ;
};

#endif 
