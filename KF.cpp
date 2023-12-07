#include "KF.h"
#include "cal.h"
#include "data.h"

KalmanFilter::KalmanFilter(MatrixXd initial_state, MatrixXd initial_covariance)
{
	A_ = MatrixXd::Identity(initial_state.rows(), initial_state.rows());
	H_ = MatrixXd::Identity(initial_state.rows(), initial_state.rows());
	x_hat_ = initial_state;
	P_ = initial_covariance;
}

void KalmanFilter::predict()
{
	x_hat_minus_ = A_ * x_hat_;
	P_minus_ = A_ * P_ * A_.transpose() + Q_;
}

void KalmanFilter::update(MatrixXd z)
{
	K_ = P_minus_ * H_.transpose() * (H_ * P_minus_ * H_.transpose() + R_).inverse();
	x_hat_ = x_hat_minus_ + K_ * (z - H_ * x_hat_minus_);
	P_ = (MatrixXd::Identity(x_hat_.rows(), x_hat_.rows()) - K_ * H_) * P_minus_;
}

void KalmanFilter::set_A(MatrixXd A)
{
	A_ = A;
}

void KalmanFilter::set_H(MatrixXd H)
{
	H_ = H;
}

void KalmanFilter::set_Q(MatrixXd Q)
{
	Q_ = Q;
}

void KalmanFilter::set_R(MatrixXd R)
{
	R_ = R;
}

void KalmanFilter::setState(MatrixXd state)
{
	x_hat_ = state;
}

MatrixXd KalmanFilter::getState() const
{
	return x_hat_;
}

MatrixXd getA(double T)
{
	MatrixXd A = MatrixXd::Identity(8, 8);
	A.block(0, 4, 4, 4) = MatrixXd::Identity(4, 4) * T;
	return A;
}

MatrixXd getH(MatrixXd B_pos, MatrixXd B_Vel)
{
	MatrixXd H = MatrixXd::Zero(B_pos.rows() + B_Vel.rows(), B_pos.cols() + B_Vel.cols());
	H.block(0, 0, B_pos.rows(), B_pos.cols()) = B_pos;
	H.block(B_pos.rows(), B_pos.cols(), B_Vel.rows(), B_Vel.cols()) = B_Vel;
	return H;
}

MatrixXd getQ(double T)
{
	MatrixXd Q = MatrixXd::Identity(8, 8);
	MatrixXd Sv = MatrixXd::Identity(3, 3) * 0.05;
	double St = 0.05;
	double Sf = 0.05;
	Q.block(0, 0, 3, 3) = Sv * pow(T, 3) / 3;
	Q.block(0, 4, 3, 3) = Sv * pow(T, 2) / 2;
	Q.block(4, 0, 3, 3) = Sv * pow(T, 2) / 2;
	Q(3, 3) = St * T + Sf * pow(T, 3) / 3;
	Q(3, 7) = Sf * pow(T, 2) / 2;
	Q(7, 3) = Sf * pow(T, 2) / 2;
	Q(7, 7) = Sf * T;
	return Q;
}

MatrixXd getR(double ROW_P, double ROW_V)
{
	MatrixXd R = MatrixXd::Identity(ROW_P + ROW_V, ROW_P + ROW_V);
	R.block(0, 0, ROW_P, ROW_P) = MatrixXd::Identity(ROW_P, ROW_P) * 4.5;
	R.block(ROW_P, ROW_P, ROW_V, ROW_V) = MatrixXd::Identity(ROW_V, ROW_V) * 0.18;
	return R;
}

MatrixXd getz(MatrixXd l_P, MatrixXd l_V)
{
	MatrixXd z = MatrixXd::Zero(l_P.rows() + l_V.rows(), 1);
	z.block(0, 0, l_P.rows(), 1) = l_P;
	z.block(l_P.rows(), 0, l_V.rows(), 1) = l_V;
	return z;
}

MatrixXd Result_DATA::Initial_KF()
{
	MatrixXd B, l_Pos, l_Vel;
	XYZ* sate_pos0 = new XYZ();
	XYZ* sate_pos1 = new XYZ();
	double clk0 = 0;
	double clk1 = 0;
	double velocity[4] = { 0, 0, 0, 0 };
	XYZ* RcvPos = new XYZ(KF->getState()(0, 0), KF->getState()(1, 0), KF->getState()(2, 0));
	double dt_Rcv = KF->getState()(3, 0);
	MatrixXd B_new = MatrixXd::Zero(1, 4);
	MatrixXd l_new = MatrixXd::Zero(1, 1);
	for (auto Sate : range->GPS_SATE)
	{
		if (Sate->Phase_NUM < 2)
			continue;
		if (Sate->Outlier)
			continue;
		int prn = Sate->PRN;
		if (GPS_eph[prn - 1]->PRN != prn)
			continue;
		double IF = 0;
		int Index1 = 0;
		int Index2 = 0;
		while (CODE2FREQ(decode_SYN(Sate->SYS, Sate->SYG_TYPE[Index1])) != L1 && Index1 < MAXNUM)
			Index1++;
		while (CODE2FREQ(decode_SYN(Sate->SYS, Sate->SYG_TYPE[Index2])) != L2 && Index2 < MAXNUM)
			Index2++;
		if (Index1 == MAXNUM || Index2 == MAXNUM)
			continue;
		if (!(Sate->LOCK_PSE[Index1] && Sate->LOCK_PSE[Index2] && Sate->LOCK_PHA[Index1] && Sate->LOCK_PHA[Index2]))
			continue;

		// ��������λ�á��Ӳ�
		double ts = OBSTIME->SecOfWeek - Sate->PSERA[0] / velocity_c;
		double dt = abs(ts - GPS_eph[prn - 1]->toe_tow + (OBSTIME->Week - GPS_eph[prn - 1]->toe_wn) * 604800);
		if (dt > 14400)
			continue;

		for (int j = 0; j < 3; j++)
		{
			clk0 = CORRECT_CLK(ts - clk0, GPS_eph[prn - 1]);
			SAT_POS_CAL(ts - clk0, GPS_eph[prn - 1], sate_pos0, clk0, Sate->PSERA[0] / velocity_c, Sate->SYS);
			clk1 = CORRECT_CLK(ts - clk1 + 1e-3, GPS_eph[prn - 1]);
			SAT_POS_CAL(ts - clk0 + 1e-3, GPS_eph[prn - 1], sate_pos0, clk0, Sate->PSERA[0] / velocity_c, Sate->SYS);
		}
		velocity[0] = (sate_pos1->X - sate_pos0->X) / 1e-3;
		velocity[1] = (sate_pos1->Y - sate_pos0->Y) / 1e-3;
		velocity[2] = (sate_pos1->Z - sate_pos0->Z) / 1e-3;
		velocity[3] = (clk1 - clk0) * velocity_c / 1e-3;
		if (Ele_Angle(sate_pos0, RcvPos, Sate->SYS) < degree2rad(10))
			continue;
		if (Sate->SYS == SYS_GPS)
			IF = SQR(L1) * Sate->PSERA[Index1] / (SQR(L1) - SQR(L2)) - SQR(L2) * Sate->PSERA[Index2] / (SQR(L1) - SQR(L2));
		else if (Sate->SYS == SYS_BDS)
		{
			double k_1_3 = SQR(L1 / L2);
			IF = (Sate->PSERA[Index2] - k_1_3 * Sate->PSERA[Index1]) / (1 - k_1_3) + velocity_c * k_1_3 * GPS_eph[prn - 1]->T_GD1 / (1 - k_1_3);
		}
		double lamda = (1e-6 * velocity_c / L1);
		double len = sqrt(SQR(RcvPos->X - sate_pos0->X) + SQR(RcvPos->Y - sate_pos0->Y) + SQR(RcvPos->Z - sate_pos0->Z));
		double w_pos = IF - (len + dt_Rcv - velocity_c * clk0 + Hopefield(sate_pos0, RcvPos, Sate->SYS));
		if (abs(w_pos) > 10)
			continue;
		double v0 = ((sate_pos0->X - RcvPos->X) * velocity[0] + (sate_pos0->Y - RcvPos->Y) * velocity[1] + (sate_pos0->Z - RcvPos->Z) * velocity[2]) / len;
		double w_Vel = -lamda * Sate->DOPPLER[Index1] - (v0 - velocity[3]);
		if (abs(w_Vel) > 10)
			continue;
		double l = (RcvPos->X - sate_pos0->X) / len;
		double m = (RcvPos->Y - sate_pos0->Y) / len;
		double n = (RcvPos->Z - sate_pos0->Z) / len;
		B_new(0, 0) = l;
		B_new(0, 1) = m;
		B_new(0, 2) = n;
		B_new(0, 3) = 1;
		l_new(0, 0) = w_pos;
		B.conservativeResize(B.rows() + 1, B.cols());
		B.bottomRows(1) = B_new;
		l_Pos.conservativeResize(l_Pos.rows() + 1, l_Pos.cols());
		l_Pos.bottomRows(1) = l_new;
		l_new(0, 0) = w_Vel;
		l_Vel.conservativeResize(l_Vel.rows() + 1, l_Vel.cols());
		l_Vel.bottomRows(1) = l_new;
	}
	delete sate_pos0;
	delete RcvPos;
	if (B.rows() == 0)
		return MatrixXd::Identity(1, 1) * -1;
	KF->set_H(getH(B, B));
	KF->set_R(getR(B.rows(), B.rows()));
	return getz(l_Pos, l_Vel);
}

unsigned int KF_SPP(Result_DATA* data, double dt_e, bool first_flag)
{
	data->OBSTIME = data->range->OBS_TIME;
	if (data->range->GPS_SATE.size() + data->range->BDS_SATE.size() < 5)
	{
		data->solve_result = OBS_DATA_Loss;
		return -1;
	}
	for (int i = 0; i < data->range->GPS_SATE.size(); i++)
	{
		if (data->range->GPS_SATE[i]->Phase_NUM < 2)
			continue;
		int prn = data->range->GPS_SATE[i]->PRN;
		int Index1 = 0;
		int Index2 = 0;
		while (CODE2FREQ(decode_SYN(data->range->GPS_SATE[i]->SYS, data->range->GPS_SATE[i]->SYG_TYPE[Index1])) != L1 && Index1 < MAXNUM)
			Index1++;
		while (CODE2FREQ(decode_SYN(data->range->GPS_SATE[i]->SYS, data->range->GPS_SATE[i]->SYG_TYPE[Index2])) != L2 && Index2 < MAXNUM)
			Index2++;

		if (!DetectOutlier(data->range->GPS_SATE[i], data->range->GPS_SATE[i]->SYS, dt_e, Index1, Index2) || Index1 == MAXNUM || Index2 == MAXNUM)
		{
			data->range->GPS_SATE[i]->Outlier = true;
			continue;
		}
	}
	for (int i = 0; i < data->range->BDS_SATE.size(); i++)
	{
		if (data->range->BDS_SATE[i]->Phase_NUM < 2)
			continue;
		int prn = data->range->BDS_SATE[i]->PRN;
		int Index1 = 0;
		int Index2 = 0;
		while (CODE2FREQ(decode_SYN(data->range->BDS_SATE[i]->SYS, data->range->BDS_SATE[i]->SYG_TYPE[Index1])) != B1 && Index1 < MAXNUM)
			Index1++;
		while (CODE2FREQ(decode_SYN(data->range->BDS_SATE[i]->SYS, data->range->BDS_SATE[i]->SYG_TYPE[Index2])) != B3 && Index2 < MAXNUM)
			Index2++;

		if (!DetectOutlier(data->range->BDS_SATE[i], data->range->BDS_SATE[i]->SYS, dt_e, Index1, Index2) || Index1 == MAXNUM || Index2 == MAXNUM)
		{
			data->range->BDS_SATE[i]->Outlier = true;
			continue;
		}
	}
	if (first_flag)
	{
		if (!Cal_1(data, first_flag))
			return 0;
		MatrixXd state(8, 1);
		state.block(0, 0, 3, 1) = data->Pos->block(0, 0, 3, 1);
		state.block(3, 0, 3, 1) = data->Vel->block(0, 0, 3, 1);
		state(6, 0) = (*data->Pos)(3, 0);
		state(7, 0) = (*data->Vel)(3, 0);
		data->KF->setState(state);
	}

	if (KF_1(data, first_flag, dt_e))
		data->solve_result = Success;
	else
		return 0;

	// if (SYS_num == 1)
	// {
	// 	switch (User_SYS)
	// 	{
	// 	case SYS_GPS:
	// 		if (Cal_1(data, first_flag))
	// 			data->solve_result = Success;
	// 		else
	// 			return 0;
	// 		break;
	// 	case SYS_BDS:
	// 		if (Cal_1(data, first_flag))
	// 			data->solve_result = Success;
	// 		else
	// 			return 0;
	// 	default:
	// 		break;
	// 	}
	// }
	// else if (SYS_num == 2)
	// {
	// 	if (Cal_2(data, first_flag))
	// 		data->solve_result = Success;
	// 	else
	// 		return 0;
	// }

	return 1;
}

unsigned int KF_1(Result_DATA* data, bool first_flag, double T)
{
	data->KF->set_A(getA(T));
	data->KF->set_Q(getQ(T));
	data->KF->predict();
	MatrixXd z = data->Initial_KF();
	if (z(0, 0) == -1)
		return 0;
	data->KF->update(z);

	return 1;
}