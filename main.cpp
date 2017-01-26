#include "parser.h"
#include "fit.h"

using namespace std;

vector<HyperCube> cubes;
vector<Point> points;

int main(int argc, char *argv[]) {
	char* file_name_1 = argv[1];
	char* file_name_2 = argv[2];
	Parser myparser = Parser(file_name_1,file_name_2);
	LeastSquare myls = LeastSquare(myparser);
	cubes = myparser.cubes;
	points = myparser.points;
	/* vector<int> temp_list;
	for (size_t i=0; i<cubes.size(); i++) {
		temp_list.push_back(cubes[i].label_);
	}
	sort(temp_list.begin(),temp_list.end());
	temp_list.erase(unique(temp_list.begin(),temp_list.end()),temp_list.end());
	cout << temp_list.size(); */
	ofstream outp("Test");
	/* // Cubes output
	for (size_t i=0; i<cubes.size(); i++) {
		outp<< cubes[i].label_;
		for (size_t j=0; j<cubes[i].points_.back()->value_.size()-1; j++) {
			valarray<double> temp(cubes[i].points_.size());
			for (size_t k=0; k<cubes[i].points_.size();k++) {
				temp[k] = cubes[i].points_[k]->value_[j];
			}
			double out_temp = temp.sum()/cubes[i].points_.size();
			outp << " "<< j+1 << ":" << out_temp;
		}
		outp <<endl;
	}
	// Points output
	for (size_t i=0; i<points.size(); i++) {
		outp<< points[i].label_;
		for (size_t j=0; j<points[i].value_.size()-1; j++) {
			outp << " "<< j+1 << ":" << points[i].value_[j];
		}
		outp <<endl;
	} */
	for (size_t i=0; i<points.size(); i++) {
		outp<< points[i].fitting_value_<<endl;
	} 
	outp.close();
	return 0;
} 
