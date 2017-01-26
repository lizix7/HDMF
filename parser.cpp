#include "parser.h"

using namespace std;

Parser::Parser(char* file_name_1, char* file_name_2) {
	vector<size_t> index;
	ifstream inFile(file_name_1);
	if (inFile.fail()) {
		cerr << "Problem occurred with input file." << endl;
		inFile.close();
		exit(1);
	}
	string line;
	while(getline(inFile, line)) {
		Parse(line);
	}
	inFile.close();
	ifstream infile(file_name_2);
	if (infile.fail()) {
		cerr << "Problem occurred with input file." << endl;
		infile.close();
		exit(1);
	}
	while(getline(infile, line)) {
		ParseTotal(line);
	}
	infile.close();
	for (size_t i=0; i<points.size(); i++) {
		points[i].total_capacitance_=points_total[i].value_.back();
	}
	Analyse(points);
	Allocate(reverse_sizes,0,index);
	return;
}

// destructor
Parser::~Parser() {
}
// This function parse the data into points vector
void Parser::Parse(string line) {
	Point point;
	double value;
	istringstream streaml (line, istringstream::in);
	while(streaml >> value) {
		point.value_.push_back(value);
	}
	streaml.clear();
	point.label_=0;
	points.push_back(point);
	return;
}
// This function parse the data into points vector
void Parser::ParseTotal(string line) {
	Point point;
	double value;
	istringstream streamt (line, istringstream::in);
	while(streamt >> value) {
		point.value_.push_back(value);
	}
	streamt.clear();
	points_total.push_back(point);
	return;
}
// This function helps to generate coefficient vector
void Parser::MeshArea(size_t position, int output_1, int output_2) {
	// base condition
	if (position ==0) {
		cube_coefficients.push_back(output_1);
		point_coefficients.push_back(output_2);
	}
	else {
		position--;
		output_1 = output_1*(sizes[position]-1);
		output_2 = output_2*(sizes[position]);
		MeshArea(position,output_1,output_2);
	}
	return;
}
// This function analyse the points and form a mesh vector and coefficient vector
void Parser::Analyse(vector<Point> points) {
	vector<double> axis;
	size_t dimension = points.back().value_.size();
	axis.reserve(points.size());
	for (size_t i=0; i<dimension-1; i++) {
		for (size_t j=0; j<points.size(); j++) {
			axis.push_back(points[j].value_[i]);
		}
		DeleteDuplicate(axis);
		sizes.push_back(axis.size());
		reverse_sizes.push_back(axis.size());
		axis.clear();
	}
	for (size_t k=0; k<sizes.size(); k++) {
		MeshArea(k,1,1);
	}
	reverse(reverse_sizes.begin(),reverse_sizes.end());
	return;
}
// This function help to find the hypercubes' points
void Parser::Combination(vector<vector<size_t> >indexs, size_t i, vector<size_t> combination) {
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
// This function allocates points to hypercubes
void Parser::Allocate(vector<size_t> reverse_sizes, size_t a, vector<size_t> index) {
	// base condition
	if (a==reverse_sizes.size()){
		reverse(index.begin(),index.end());
		vector<size_t> temp_a;
		vector<size_t> combination;
		cube.points_.clear();
		index_combinations.clear();
		combinations.clear();
		for (size_t t=0; t<index.size(); t++) {
			temp_a.push_back(index[t]);
			temp_a.push_back(index[t]+1);
			index_combinations.push_back(temp_a);
			temp_a.clear();
		}
		Combination(index_combinations,0,combination);
		for (size_t c=0; c<combinations.size(); c++) {
			size_t position = inner_product(combinations[c].begin(),combinations[c].end(),point_coefficients.begin(),0);
			cube.points_.push_back(&points[position]);
		}
		cube.label_ = 0;
		cubes.push_back(cube);
	}
	else {
		for (size_t i=0; i<reverse_sizes[a]-1; i++) {	
			vector<size_t> temp(index);
			temp.push_back(i);
			Allocate(reverse_sizes,a+1,temp);
		}	
	}
	return;
}
// This function deletes duplicate elements
void Parser::DeleteDuplicate(vector<double> &cord) {
	sort(cord.begin(),cord.end());
	cord.erase(unique(cord.begin(),cord.end()),cord.end());
	return;
}
