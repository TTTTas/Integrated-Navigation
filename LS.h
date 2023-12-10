#pragma once
#include <Eigen/Dense>
#include <iostream>
using namespace std;
using namespace Eigen;

class Least_Squares
{
public:
	MatrixXd X;
	MatrixXd B;
	MatrixXd l;
	MatrixXd P;
	MatrixXd Qxx;
	double sigma;

	void LS();
};