#ifndef GLOBAL_H
#define GLOBAL_H

#include <string>
#include <vector>
#include <iostream>
#include <cmath>
#include <numeric>
#include <cstdlib>
#include <stdio.h>
#include <algorithm>
#include <valarray>
#include <fstream>
#include <sstream>
#include <ctime>
#include <map>
#include <sys/time.h>

using namespace std;

struct Point {
	// value_ stores all the data from each dimensions in a vector
	// e.g. x,y,z,t,q.
	vector<double> value_;
	vector<int> labels_;
	int label_;
	double total_capacitance_;
	double fitting_value_;
};

struct Equation {
	// parameters_ store all the coefficients of the equation in a vector
	// t=exp(ax+by+cz+dt+eq+...)
	vector<double> parameters_;
	int label_; //label 
};

struct HyperCube {
	// points_ stores all the points of a hypercube in a vector
	vector<Point*> points_;
	int label_; 
};

#endif 
