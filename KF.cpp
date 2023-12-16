#include "KF.h"
#include "cal.h"
#include "data.h"

KalmanFilter::KalmanFilter(MatrixXd initial_state, MatrixXd initial_covariance)
{
	A_ = MatrixXd::Identity(initial_state.rows(), initial_state.rows());
	H_ = MatrixXd::Zero(0, 0);
	delt_z = MatrixXd::Zero(0, 0);
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
	x_hat_ = x_hat_minus_ + K_ * delt_z;
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
		MatrixXd B_pos = H_.block(0, 0, last_row, Last_col - 4);
		MatrixXd B_vel = H_.block(last_row, Last_col - 4, last_row, 4);
		MatrixXd B_pos_add = H.block(0, 0, temp_row, temp_col - 4);
		MatrixXd B_vel_add = H.block(temp_row, temp_col - 4, temp_row, 4);

		int last_col = B_pos.cols();
		B_pos.conservativeResize(last_row + temp_row, last_col + 1);
		B_pos.block(last_row, 0, temp_row, 3) = B_pos_add.leftCols(3);
		B_pos.block(0, last_col, last_row, 1) = MatrixXd::Zero(last_row, 1);
		B_pos.block(last_row, last_col, temp_row, 1) = MatrixXd::Constant(temp_row, 1, 1);
		B_pos.block(last_row, 3, temp_row, last_col - 3) = MatrixXd::Zero(temp_row, last_col - 3);

		B_vel.conservativeResize(B_vel.rows() + B_vel_add.rows(), 4);
		B_vel.bottomRows(B_vel_add.rows()) = B_vel_add;

		H_.conservativeResize(B_pos.rows() + B_vel.rows(), B_pos.cols() + B_vel.cols());
		H_.block(0, 0, B_pos.rows(), B_pos.cols()) = B_pos;
		H_.block(B_pos.rows(), B_pos.cols(), B_vel.rows(), B_vel.cols()) = B_vel;
		H_.block(0, B_pos.cols(), B_pos.rows(), B_vel.cols()) = MatrixXd::Zero(B_pos.rows(), B_vel.cols());
		H_.block(B_pos.rows(), 0, B_vel.rows(), B_pos.cols()) = MatrixXd::Zero(B_vel.rows(), B_pos.cols());
	}
}

void KalmanFilter::set_Z(MatrixXd z)
{
	if (delt_z.rows() == 0)
		delt_z = z;
	else
	{
		delt_z.conservativeResize(delt_z.rows() + z.rows(), 1);
		delt_z.bottomRows(z.rows()) = z;
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
	x_hat_minus_ = state;
}

void KalmanFilter::reset()
{
	H_ = MatrixXd::Zero(0, 0);
	delt_z = MatrixXd::Zero(0, 0);
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
	A.block(0, row - 4, 4, 4) = MatrixXd::Identity(4, 4) * T;
	return A;
}

MatrixXd getH(MatrixXd B)
{
	MatrixXd H = MatrixXd::Zero(B.rows() + B.rows(), B.cols() + B.cols());
	H.block(0, 0, B.rows(), B.cols()) = B;
	H.block(B.rows(), B.cols(), B.rows(), B.cols()) = B;
	return H;
}

MatrixXd getQ(double T, Configure cfg)
{
	MatrixXd Q;
	MatrixXd Sv = MatrixXd::Identity(3, 3) * 10;
	double St = 1;
	double Sf = 1;
	if (cfg.SYS_num == 1)
	{
		Q = MatrixXd::Identity(8, 8);
		Q.block(0, 0, 3, 3) = Sv * pow(T, 3) / 3;
		Q.block(0, 4, 3, 3) = Sv * pow(T, 2) / 2;
		Q.block(4, 0, 3, 3) = Sv * pow(T, 2) / 2;
		Q.block(4, 4, 3, 3) = Sv * T;
		Q(3, 3) = St * T + Sf * pow(T, 3) / 3;
		Q(3, 7) = Sf * pow(T, 2) / 2;
		Q(7, 3) = Sf * pow(T, 2) / 2;
		Q(7, 7) = Sf * T;
	}
	else if (cfg.SYS_num == 2)
	{
		Q = MatrixXd::Identity(9, 9);
		Q.block(0, 0, 3, 3) = Sv * pow(T, 3) / 3;
		Q.block(0, 5, 3, 3) = Sv * pow(T, 2) / 2;
		Q.block(5, 0, 3, 3) = Sv * pow(T, 2) / 2;
		Q.block(4, 4, 3, 3) = Sv * T;
		Q(3, 3) = St * T + Sf * pow(T, 3) / 3;
		Q(3, 8) = Sf * pow(T, 2) / 2;
		Q(8, 3) = Sf * pow(T, 2) / 2;
		Q(4, 4) = St * T + Sf * pow(T, 3) / 3;
		Q(4, 8) = Sf * pow(T, 2) / 2;
		Q(8, 4) = Sf * pow(T, 2) / 2;
		Q(8, 8) = Sf * T;
	}
	return Q;
}

MatrixXd getR(double ROW)
{
	MatrixXd R = MatrixXd::Identity(ROW + ROW, ROW + ROW);
	R.block(0, 0, ROW, ROW) = MatrixXd::Identity(ROW, ROW) * 2.5;
	R.block(ROW, ROW, ROW, ROW) = MatrixXd::Identity(ROW, ROW) * 0.05;
	return R;
}

MatrixXd getz(MatrixXd l_P, MatrixXd l_V)
{
	MatrixXd z = MatrixXd::Zero(l_P.rows() + l_V.rows(), 1);
	z.block(0, 0, l_P.rows(), 1) = l_P;
	z.block(l_P.rows(), 0, l_V.rows(), 1) = l_V;
	return z;
}