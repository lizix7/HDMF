#include "fit.h"

using namespace std;
const double kError = 0.03;

//constructor
LeastSquare::LeastSquare(Parser &MyParser) {
	vector<Point*> points;
	point_coefficients = MyParser.point_coefficients;
	cube_coefficients = MyParser.cube_coefficients;
	for (size_t i=0; i<MyParser.points.size(); i++) {
		total_points.push_back(&MyParser.points[i]);
	}
	vector<size_t> index;
	int label = 0;
	// Label all the cubes
	Traverse(MyParser.reverse_sizes,MyParser.sizes,0,index,label,MyParser.cubes,points,total_points,PointMap,EquationMap);
	// Label all the points
	// Push back the labels for each point
	for (size_t i=0; i<MyParser.cubes.size(); i++) {
		for (size_t j=0; j<MyParser.cubes[i].points_.size(); j++) {
			MyParser.cubes[i].points_[j]->labels_.push_back(MyParser.cubes[i].label_);
		}
	}
	// Find the minimum label for the point
	for (size_t i=0; i<MyParser.points.size(); i++) {
		// delete duplicate begin
		sort(MyParser.points[i].labels_.begin(),MyParser.points[i].labels_.end());
		MyParser.points[i].labels_.erase(unique(MyParser.points[i].labels_.begin(),MyParser.points[i].labels_.end()),MyParser.points[i].labels_.end());
		// delete duplicate done
		double value = MyParser.points[i].value_.back();
		vector<double> values = MyParser.points[i].value_;
		values.pop_back();
		values.push_back(1.0);
		map<double,int> ErrorMap;
		vector<double> ErrorList;
		for (size_t j=0; j<MyParser.points[i].labels_.size(); j++) {
			Equation equation = EquationMap[MyParser.points[i].labels_[j]];
			double error = abs(value-exp(inner_product(values.begin(),values.end(),equation.parameters_.begin(),0.0)));
			ErrorMap[error] = MyParser.points[i].labels_[j];
			ErrorList.push_back(error);
		}
		double min_error = *min_element(ErrorList.begin(),ErrorList.end());
		int min_label = ErrorMap[min_error];
		MyParser.points[i].label_ = min_label;
		Equation new_equation = EquationMap[min_label];
		MyParser.points[i].fitting_value_ = exp(inner_product(values.begin(),values.end(),new_equation.parameters_.begin(),0.0));
	}
}

// destructor
LeastSquare::~LeastSquare() {
}
// Check the hypercube itself
void LeastSquare::SelfCheck(int &label, vector<size_t> index, vector<HyperCube> &cubes, vector<Point*> &points, map<int,vector<Point*> > PointMap) {
	size_t position = inner_product(index.begin(),index.end(),cube_coefficients.begin(),0);
	if(cubes[position].label_==0) {
		label++;
		cubes[position].label_ = label;
		points.clear();
		points = cubes[position].points_;
	}
	else
		points = PointMap[cubes[position].label_];
}

// This function help to find the hypercubes' points
void LeastSquare::Combination(vector<vector<size_t> >indexs, size_t i, vector<size_t> combination) {
	if (i==indexs.size()) {
		combinations.push_back(combination);
	}
	else {
		vector<size_t> row = indexs[i];
		for (size_t j=0; j<row.size(); j++) {
			vector<size_t> temp_c(combination);
			temp_c.push_back(row[j]);
			Combination(indexs,i+1,temp_c);
		}
	}
	return;
}
// Solve the dense matrix by Gaussian elimination method
Equation LeastSquare::Solve_Matrix(vector<vector<double> > &A, vector<double> &B) {
	Equation X;
	size_t dimension = A.size();
	vector<double> x(dimension);
	for (size_t i=0; i<dimension; i++) {
		for (size_t j=i+1; j<dimension; j++) {
			double ratio = A[j][i]/A[i][i];
			A[j][i] = 0.0;
			for (size_t k=i+1; k<dimension; k++ ) {
				A[j][k] = A[j][k] - ratio*A[i][k];
			}
			B[j] = B[j] - ratio*B[i];
		}
	}
	for (int i=dimension-1; i>=0; i--) {
		double sum = 0.0;
		for (size_t j=i+1;j<dimension;j++) {
			sum = sum + A[i][j]*x[j];
		}
		x[i] = (B[i]-sum)/A[i][i];
	}
	X.parameters_=x;
	return X;
}
// Calculate the fitting model
Equation LeastSquare::Calculate(vector<Point*> points) {
	size_t matrix_dimension = points.back()->value_.size();
	vector<valarray<double> > B(matrix_dimension,valarray<double>(points.size()));
	vector<vector<valarray<double> > > A(matrix_dimension, B);
	for (size_t k=0; k<points.size(); k++) {
		for (size_t j=0; j<points[k]->value_.size(); j++) {
			double value = points[k]->value_.back();	
			double ln_value = log(value);
			if (j==points[k]->value_.size()-1) {
				B[j][k] = value*ln_value;
			}
			else {
				B[j][k] = points[k]->value_[j]*value*ln_value;
			}
			for (size_t i=0; i<points[k]->value_.size(); i++) {
				if(i==points[k]->value_.size()-1 && j!=points[k]->value_.size()-1) {
					A[j][i][k] = points[k]->value_[j]*value;
					
				}
				else if(i!=points[k]->value_.size()-1 && j==points[k]->value_.size()-1) {
					A[j][i][k] = points[k]->value_[i]*value;
				}
				else if (i==points[k]->value_.size()-1 && j==points[k]->value_.size()-1) {
					A[j][i][k] = value;
				}
				else {
					A[j][i][k] = points[k]->value_[i]*points[k]->value_[j]*value;
				}
			}
		}
	}
	vector<vector<double> > new_A(matrix_dimension,vector<double>(matrix_dimension));
	vector<double> new_B(matrix_dimension);
	for (size_t j=0; j<matrix_dimension; j++) {
		for (size_t i=0; i<matrix_dimension; i++) {
			new_A[j][i] = A[j][i].sum();
		}
		new_B[j] = B[j].sum();
	}
	Equation x = Solve_Matrix(new_A,new_B);
	return x;
}
// Check the single neighbour of the hypercube
void LeastSquare::CheckNeighbour(vector<size_t> index, vector<size_t> new_index, vector<vector<size_t> > index_combinations, vector<HyperCube> &cubes, vector<Point*> &points, vector<Point*> total_points, map<int,vector<Point*> > PointMap) {
	vector<size_t> combination;
	combinations.clear();
	size_t self_position = inner_product(index.begin(),index.end(),cube_coefficients.begin(),0);
	size_t neighbour_position = inner_product(new_index.begin(),new_index.end(),cube_coefficients.begin(),0);
	Combination(index_combinations,0,combination);
	if (cubes[neighbour_position].label_==0) {
		for (size_t c=0; c<combinations.size(); c++) {
			size_t point_position = inner_product(combinations[c].begin(),combinations[c].end(),point_coefficients.begin(),0);
			points.push_back(total_points[point_position]);
		}
		Equation equation = Calculate(points);
		bool flag = CheckError(points,equation);
		if (flag) {
			cubes[neighbour_position].label_=cubes[self_position].label_;
		}
		else {
			for (size_t c=0; c<combinations.size(); c++) {
				points.pop_back();
			}
		}
	}
	else if (cubes[neighbour_position].label_!=0 && cubes[neighbour_position].label_!=cubes[self_position].label_) {
		for (size_t i=0; i<PointMap[cubes[neighbour_position].label_].size(); i++) {
			points.push_back(PointMap[cubes[neighbour_position].label_][i]);
		}
		Equation equation = Calculate(points);
		bool flag = CheckError(points,equation);
		if (flag) {
			int label_1 = cubes[neighbour_position].label_;
			int label_2 = cubes[self_position].label_;
			for (size_t i=0; i<cubes.size(); i++) {
				if (cubes[i].label_==label_1) {
					cubes[i].label_=label_2;
				}
			}
		}
		else {
			for (size_t i=0; i<PointMap[cubes[neighbour_position].label_].size(); i++) {
				points.pop_back();
			}
		}
	}
}
// Check all the neighbours of the hypercube
void LeastSquare::CheckNeighbours(vector<size_t> sizes, vector<size_t> index, vector<HyperCube> &cubes, vector<Point*> &points, vector<Point*> total_points, map<int,vector<Point*> > PointMap) {
	vector<size_t> temp;
	vector<size_t> combination;
	index_combinations.clear();
	for (size_t i=0; i<index.size(); i++) {
		temp.push_back(index[i]);
		temp.push_back(index[i]+1);
		index_combinations.push_back(temp);
		temp.clear();
	}
	for (size_t i=0; i<index.size(); i++) {
		temp = index_combinations[i];
		vector<size_t> temp_index(index);
		if (index[i]!=0 && index[i]!=sizes[i]-2) {
			index_combinations[i].clear();
			index_combinations[i].push_back(index[i]-1);
			temp_index[i] = index[i]-1;
			CheckNeighbour(index,temp_index,index_combinations,cubes,points,total_points,PointMap);
			index_combinations[i].clear();
			index_combinations[i].push_back(index[i]+2);
			temp_index[i] = index[i]+1;
			CheckNeighbour(index,temp_index,index_combinations,cubes,points,total_points,PointMap);
		}
		else if (index[i]==0 && index[i]!=sizes[i]-2) {
			index_combinations[i].clear();
			index_combinations[i].push_back(index[i]+2);
			temp_index[i] = index[i]+1;
			CheckNeighbour(index,temp_index,index_combinations,cubes,points,total_points,PointMap);
		}
		else if (index[i]!=0 && index[i]==sizes[i]-2) {
			index_combinations[i].clear();
			index_combinations[i].push_back(index[i]-1);
			temp_index[i] = index[i]-1;
			CheckNeighbour(index,temp_index,index_combinations,cubes,points,total_points,PointMap);
		}
		else 
			cout<< "The data is only 1 dimension, please check your initial data!"<<endl;
		index_combinations[i] = temp;
	}
}

bool BasicSort(Point* lhs, Point* rhs, size_t i) {
	if (i==lhs->value_.size()-1) {
		return lhs->value_[i]<rhs->value_[i];
	}
	else {
		if (lhs->value_[i]==rhs->value_[i]) {
			bool result = BasicSort(lhs,rhs,i+1);
			return result;
		}
		else
			return lhs->value_[i]<rhs->value_[i];
	}
}

bool SortHD(Point* lhs, Point* rhs) {
	bool result = BasicSort(lhs,rhs,0);
	return result;
}

bool UniqueHD(Point* lhs, Point* rhs) {
	if (lhs->value_==rhs->value_)
		return true;
	else
		return false;
}
// Delete duplicate points
void LeastSquare::DeleteDuplicate(vector<Point*> &points) {
	sort(points.begin(),points.end(),SortHD);
	points.erase(unique(points.begin(),points.end(),UniqueHD),points.end());
}
// Check the error<3%
bool LeastSquare::CheckError(vector<Point*> points, Equation equation) {
	for (size_t i=0; i<points.size();i++) {
		double value = points[i]->value_.back();
		double total_capacitance = points[i]->total_capacitance_;
		vector<double> values = points[i]->value_;
		values.pop_back();
		values.push_back(1.0);
		double self_error = abs(value-exp(inner_product(values.begin(),values.end(),equation.parameters_.begin(),0.0)));
		//cout << self_error/value*100;
		double error = kError*total_capacitance;
		if (self_error>error) {
			return false;
		}
	}
	return true;
}
// Update the vector
void LeastSquare::Update(vector<size_t> index, vector<HyperCube> cubes, vector<Point*> &points, map<int,vector<Point*> > &PointMap, map<int,Equation> &EquationMap) {
	size_t position = inner_product(index.begin(),index.end(),cube_coefficients.begin(),0);
	DeleteDuplicate(points);
	Equation equation = Calculate(points);
	EquationMap[cubes[position].label_] = equation;
	PointMap[cubes[position].label_] = points;
}
// Traverse all the hypercubes
void LeastSquare::Traverse(vector<size_t> reverse_sizes, vector<size_t> sizes, size_t t, vector<size_t> index, int &label, vector<HyperCube> &cubes, vector<Point*> &points, vector<Point*> total_points, map<int,vector<Point*> > &PointMap, map<int,Equation> &EquationMap) {
	if (t==reverse_sizes.size()) {
		reverse(index.begin(),index.end());
		SelfCheck(label,index,cubes,points,PointMap);
		CheckNeighbours(sizes,index,cubes,points,total_points,PointMap);
		Update(index,cubes,points,PointMap,EquationMap);
	}
	else {
		for (size_t i=0; i<reverse_sizes[t]-1; i++) {
			vector<size_t> temp(index);
			temp.push_back(i);
			Traverse(reverse_sizes,sizes,t+1,temp,label,cubes,points,total_points,PointMap,EquationMap);
		}
	}
}
