#include "KF.h"
#include "cal.h"
#include "data.h"

KalmanFilter::KalmanFilter(MatrixXd initial_state, MatrixXd initial_covariance)
{
	A_ = MatrixXd::Identity(initial_state.rows(), initial_state.rows());
	H_ = MatrixXd::Zero(0, 0);
	z_ = MatrixXd::Zero(0, 0);
	R_ = MatrixXd::Zero(0, 0);
	x_hat_ = initial_state;
	P_ = initial_covariance;
}

void KalmanFilter::predict()
{
	x_hat_minus_ = A_ * x_hat_;
	P_minus_ = A_ * P_ * A_.transpose() + Q_;
}

void KalmanFilter::update()
{
	K_ = P_minus_ * H_.transpose() * (H_ * P_minus_ * H_.transpose() + R_).inverse();
	x_hat_ = x_hat_minus_ + K_ * (z_);
	P_ = (MatrixXd::Identity(x_hat_.rows(), x_hat_.rows()) - K_ * H_) * P_minus_;
}

void KalmanFilter::set_A(MatrixXd A)
{
	A_ = A;
}

void KalmanFilter::set_H(MatrixXd H)
{
	if(H_.rows()==0)
		H_ = H;
	else
	{
		int last_row = H_.rows() / 2;
		int Last_col = H.cols();
		int temp_row = H.rows() / 2;
		int temp_col = H.cols();
		MatrixXd B_pos = H_.block(0, 0, last_row, 3);
		MatrixXd B_vel = H_.block(last_row, 3, last_row, 3);
		int sys_num = Last_col - 7;
		MatrixXd B_old_one = H.block(0, 6, last_row, sys_num);
		MatrixXd B_pos_add = H.block(0, 0, temp_row, 3);
		MatrixXd B_vel_add = H.block(temp_row, 3, temp_row, 3);
		MatrixXd B_add_one = H.block(0, 6, temp_row, 1);

		B_pos.conservativeResize(last_row + temp_row, 3);
		B_pos.bottomRows(temp_row) = B_pos_add;

		B_vel.conservativeResize(last_row + temp_row, 3);
		B_vel.bottomRows(B_vel_add.rows()) = B_vel_add;

		H_ = MatrixXd::Zero(B_pos.rows() + B_vel.rows(), 6 + sys_num + 2);
		H_.block(0, 0, B_pos.rows(), 3) = B_pos;
		H_.block(B_pos.rows(), 3, B_vel.rows(), 3) = B_vel;
		H_.block(0, 6, last_row, sys_num) = B_old_one;
		H_.block(last_row, 6 + sys_num, temp_row, 1) = B_add_one;
		H_.block(B_pos.rows(), 7 + sys_num, temp_row, 1) = B_add_one;
	}
}

void KalmanFilter::set_Z(MatrixXd z)
{
	if (z_.rows() == 0)
		z_ = z;
	else
	{
		z_.conservativeResize(z_.rows() + z.rows(), 1);
		z_.bottomRows(z.rows()) = z;
	}
}

void KalmanFilter::set_Q(MatrixXd Q)
{
	Q_ = Q;
}

void KalmanFilter::set_R(MatrixXd R)
{
	if (R_.rows() == 0)
		R_ = R;
	else
	{
		int last = R_.rows() / 2;
		int temp = R.rows() / 2;
		MatrixXd P_pos = R_.block(0, 0, last, last);
		MatrixXd P_vel = R_.block(last, last, last, last);
		MatrixXd P_pos_add = R.block(0, 0, temp, temp);
		MatrixXd P_vel_add = R.block(temp, temp, temp, temp);

		P_pos.conservativeResize(last + temp, last + temp);
		P_pos.block(last, last, temp, temp) = P_pos_add;
		P_pos.block(0, last, last, temp) = MatrixXd::Zero(last, temp);
		P_pos.block(last, 0, temp, last) = MatrixXd::Zero(temp, last);

		P_vel.conservativeResize(last + temp, last + temp);
		P_vel.block(last, last, temp, temp) = P_vel_add;
		P_vel.block(0, last, last, temp) = MatrixXd::Zero(last, temp);
		P_vel.block(last, 0, temp, last) = MatrixXd::Zero(temp, last);

		R_.conservativeResize(P_pos.rows() + P_vel.rows(), P_pos.cols() + P_vel.cols());
		R_.block(0, 0, P_pos.rows(), P_pos.cols()) = P_pos;
		R_.block(P_pos.rows(), P_pos.cols(), P_vel.rows(), P_vel.cols()) = P_vel;
		R_.block(0, P_pos.cols(), P_pos.rows(), P_vel.cols()) = MatrixXd::Zero(P_pos.rows(), P_vel.cols());
		R_.block(P_pos.rows(), 0, P_vel.rows(), P_pos.cols()) = MatrixXd::Zero(P_vel.rows(), P_pos.cols());
	}
}

void KalmanFilter::setState(MatrixXd state)
{
	x_hat_ = state;
}

void KalmanFilter::reset()
{
	H_ = MatrixXd::Zero(0, 0);
	z_ = MatrixXd::Zero(0, 0);
	R_ = MatrixXd::Zero(0, 0);
}

MatrixXd KalmanFilter::getState() const
{
	return x_hat_;
}

MatrixXd KalmanFilter::getState_minus()
{
	return x_hat_minus_;
}

MatrixXd getA(double T, Configure cfg)
{
	int row = cfg.SYS_num + 7;
	MatrixXd A = MatrixXd::Identity(row, row);
	A.block(0, 3, 3, 3) = MatrixXd::Identity(3, 3) * T;
	if (cfg.SYS_num == 1)
		A(6, 7) = T;
	else if (cfg.SYS_num == 2)
	{
		A(6, 8) = T;
		A(7, 8) = T;
	}
	return A;
}

MatrixXd getH(MatrixXd B)
{
	MatrixXd H = MatrixXd::Zero(B.rows() + B.rows(), B.cols() + B.cols());
	H.block(0, 0, B.rows(), B.cols() - 1) = B.leftCols(3);
	H.block(0, 6, B.rows(), 1) = B.rightCols(1);
	H.block(B.rows(), 3, B.rows(), 3) = B.leftCols(3);
	H.block(B.rows(), 7, B.rows(), 1) = B.rightCols(1);
	return H;
}

MatrixXd getQ(double T, Configure cfg)
{
	MatrixXd Sv = MatrixXd::Identity(3, 3) * 0.05;
	double St = 0.05;
	double Sf = 0.05;
	int row = cfg.SYS_num + 7;
	MatrixXd Q = MatrixXd::Zero(row, row);
	Q.block(0, 0, 3, 3) = Sv * pow(T, 3) / 3;
	Q.block(0, 3, 3, 3) = Sv * pow(T, 2) / 2;
	Q.block(3, 0, 3, 3) = Sv * pow(T, 2) / 2;
	Q.block(3, 3, 3, 3) = Sv * T;
	if (cfg.SYS_num == 1)
	{
		Q(6, 6) = St * T + Sf * pow(T, 3) / 3;
		Q(6, 7) = Sf * pow(T, 2) / 2;
		Q(7, 6) = Sf * pow(T, 2) / 2;
		Q(7, 7) = Sf * T;
	}
	else if (cfg.SYS_num == 2)
	{
		Q(6, 6) = St * T + Sf * pow(T, 3) / 3;
		Q(6, 8) = Sf * pow(T, 2) / 2;
		Q(8, 6) = Sf * pow(T, 2) / 2;
		Q(7, 7) = St * T + Sf * pow(T, 3) / 3;
		Q(7, 8) = Sf * pow(T, 2) / 2;
		Q(8, 7) = Sf * pow(T, 2) / 2;
		Q(8, 8) = Sf * T;
	}
	return Q;
}

MatrixXd getR(double ROW)
{
	MatrixXd R = MatrixXd::Identity(ROW + ROW, ROW + ROW);
	R.block(0, 0, ROW, ROW) = MatrixXd::Identity(ROW, ROW) * 4.5;
	R.block(ROW, ROW, ROW, ROW) = MatrixXd::Identity(ROW, ROW) * 0.18;
	return R;
}

MatrixXd getz(MatrixXd l_P, MatrixXd l_V)
{
	MatrixXd z = MatrixXd::Zero(l_P.rows() + l_V.rows(), 1);
	z.block(0, 0, l_P.rows(), 1) = l_P;
	z.block(l_P.rows(), 0, l_V.rows(), 1) = l_V;
	return z;
}