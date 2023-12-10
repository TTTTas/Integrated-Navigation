#include"LS.h"
#include"cal.h"
#include"data.h"

void Least_Squares::LS()
{
	Qxx = (B.transpose() * P * B).inverse();
	X = Qxx * B.transpose() * P * l;
	MatrixXd v = B * X - l;
	int m = B.rows() - B.cols();
	sigma = sqrt(((v.transpose() * P * v) / m)(0, 0));
}

int DATA_SET::Set_LS()
{
	int ROWS = 0;


	return ROWS;
}