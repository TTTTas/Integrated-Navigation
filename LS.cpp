#include"LS.h"
#include"cal.h"
#include"data.h"

Least_Squares::Least_Squares()
{
	X = MatrixXd::Zero(4, 1);
}

Least_Squares::Least_Squares(Configure cfg)
{
	if (cfg.SYS_num == 1)
		X = MatrixXd::Zero(4, 1);
	else if (cfg.SYS_num == 2)
		X = MatrixXd::Zero(5, 1);
}


void Least_Squares::ELS()
{
	Qxx = (B.transpose() * P * B).inverse();
	x = Qxx * B.transpose() * P * l;
	X = X + x;
	MatrixXd v = B * x - l;
	int m = B.rows() - B.cols();
	sigma = sqrt(((v.transpose() * P * v) / m)(0, 0));
}

void Least_Squares::LS()
{
	Qxx = (B.transpose() * P * B).inverse();
	X = Qxx * B.transpose() * P * l;
	MatrixXd v = B * X - l;
	int m = B.rows() - B.cols();
	sigma = sqrt(((v.transpose() * P * v) / m)(0, 0));
}

int DATA_SET::Set_LS(Configure cfg)
{
	temp_ref.topRows(3) = LS_Pos->X.block(0, 0, 3, 1);
	double dt_G = 0;
	double dt_C = 0;
	if (cfg.SYS_num == 1)
	{
		if (cfg.GPS_Cfg.used)
			dt_G = LS_Pos->X(3, 0);
		if (cfg.BDS_Cfg.used)
			dt_C = LS_Pos->X(3, 0);
	}
	else if(cfg.SYS_num == 2)
	{
		dt_G = LS_Pos->X(3, 0);
		dt_C = LS_Pos->X(4, 0);
	}
	if (cfg.GPS_Cfg.used)
	{
		temp_ref(3, 0) = dt_G;
		setup_LS(this, cfg, SYS_GPS);
	}
	if (cfg.BDS_Cfg.used)
	{
		temp_ref(3, 0) = dt_C;
		setup_LS(this, cfg, SYS_BDS);
	}
	return LS_Pos->B.rows();
}

int Least_Squares::set_B_Pos(MatrixXd B_)
{
	if (B.rows() == 0)
	{
		B = B_;
		return B.rows();
	}
	int Last_Row = B.rows();
	int Last_Col = B.cols();
	int temp_Row = B_.rows();
	int temp_Col = B_.cols();
	B.conservativeResize(Last_Row + temp_Row, Last_Col + 1);
	B.block(Last_Row, 0, temp_Row, 3) = B_.leftCols(3);
	B.block(0, Last_Col, Last_Row, 1) = MatrixXd::Zero(Last_Row, 1);
	B.block(Last_Row, Last_Col, temp_Row, 1) = MatrixXd::Constant(temp_Row, 1, 1);
	B.block(Last_Row, 3, temp_Row, Last_Col - 3) = MatrixXd::Zero(temp_Row, Last_Col - 3);
	return B.rows();
}

int Least_Squares::set_B_Vel(MatrixXd B_)
{
	if (B.rows() == 0)
	{
		B = B_;
		return B.rows();
	}
	B.conservativeResize(l.rows() + B_.rows(), 1);
	B.bottomRows(B_.rows()) = B_;
	return B.rows();
}
int Least_Squares::set_l(MatrixXd l_)
{
	if (l.rows() == 0)
	{
		l = l_;
		return l.rows();
	}
	l.conservativeResize(l.rows() + l_.rows(), 1);
	l.bottomRows(l_.rows()) = l_;
	return l.rows();
}
int Least_Squares::set_P(MatrixXd P_)
{
	if (P.rows() == 0)
	{
		P = P_;
		return P.rows();
	}
	int Last_Row = P.rows();
	int Last_Col = P.cols();
	int temp_Row = P_.rows();
	int temp_Col = P_.cols();
	P.conservativeResize(P.rows() + P_.rows(), P.cols() + P_.cols());
	P.block(Last_Row, Last_Col, temp_Row, temp_Col) = P_;
	P.block(0, Last_Col, Last_Row, temp_Col) = MatrixXd::Zero(Last_Row, temp_Col);
	P.block(Last_Row, 0, temp_Row, Last_Col) = MatrixXd::Zero(temp_Row, Last_Col);
	return P.rows();
}

void Least_Squares::reset()
{
	B = MatrixXd::Zero(0, 0);
	l = MatrixXd::Zero(0, 0);
	P = MatrixXd::Zero(0, 0);
}