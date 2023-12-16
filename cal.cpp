#include "cal.h"
#include "data.h"

double SQR(double x)
{
	return x * x;
}

double Len(XYZ* pos)
{
	return sqrt(SQR(pos->X) + SQR(pos->Y) + SQR(pos->Z));
}

double degree2rad(double degree)
{
	while (degree > 360 || degree < 0)
	{
		if (degree > 360)
			degree -= 360;
		if (degree < 0)
			degree += 360;
	}
	return degree * Pi / 180;
}

double rad2degree(double rad)
{
	return rad * 180 / Pi;
}

double Cal_PDOP(MatrixXd Qxx)
{
	if (Qxx.rows() < 3 || Qxx.cols() < 3)
		return 0;
	return sqrt(Qxx(0, 0) + Qxx(1, 1) + Qxx(2, 2));
}

unsigned int Cal_LEAST_SQR(MatrixXd B, MatrixXd l, MatrixXd P, MatrixXd* Qxx, MatrixXd& x, double* thegma, double* DOP)
{
	if (B.rows() != P.rows() || B.rows() != l.rows())
		return 0;
	if (B.rows() < B.cols())
		return 0;
	*Qxx = (B.transpose() * P * B).inverse();
	x = *Qxx * B.transpose() * P * l;
	MatrixXd v = B * x - l;
	int m = B.rows() - B.cols();
	*thegma = sqrt(((v.transpose() * P * v) / m)(0, 0));
	*DOP = Cal_PDOP(*Qxx);
	return 1;
}

unsigned int decode_SYN(int sys, int signal)
{
	switch (sys)
	{
	case SYS_GPS:
		switch (signal)
		{
		case 0:
			return CODE_L1C; /* L1C/A */
		case 5:
			return CODE_L2P; /* L2P    (OEM7) */
		case 9:
			return CODE_L2W; /* L2P(Y),semi-codeless */
		case 14:
			return CODE_L5Q; /* L5Q    (OEM6) */
		case 16:
			return CODE_L1L; /* L1C(P) (OEM7) */
		case 17:
			return CODE_L2S; /* L2C(M) (OEM7) */
		default:
			return UNKOWN;
			break;
		}
	case SYS_BDS:
		switch (signal)
		{
		case 0:
			return CODE_L2I; /* B1I with D1 (OEM6) */
		case 1:
			return CODE_L7I; /* B2I with D1 (OEM6) */
		case 2:
			return CODE_L6I; /* B3I with D1 (OEM7) */
		case 4:
			return CODE_L2I; /* B1I with D2 (OEM6) */
		case 5:
			return CODE_L7I; /* B2I with D2 (OEM6) */
		case 6:
			return CODE_L6I; /* B3I with D2 (OEM7) */
		case 7:
			return CODE_L1P; /* B1C(P) (OEM7) */
		case 9:
			return CODE_L5P; /* B2a(P) (OEM7) */
		default:
			return UNKOWN;
			break;
		}
	default:
		return UNKOWN;
		break;
	}
}

double CODE2FREQ(int code)
{
	switch (code)
	{
	case 0:
		return 0;
	case 1:
		return L1;
	case 2:
		return L2;
	case 3:
		return L2;
	case 4:
		return L5;
	case 5:
		return L1;
	case 6:
		return L2;
	case 7:
		return B1;
	case 8:
		return B2;
	case 9:
		return B3;
	case 10:
		return B1_C;
	case 11:
		return B2_a;
	default:
		return 0;
		break;
	}
}

// ÖÓ²î¸ÄÕý
double CORRECT_CLK(double t, EPHEMERIS* eph)
{
	double correct_clk = t;
	for (int i = 0; i < 10; i++)
	{
		correct_clk = t - (eph->a_f0 + eph->a_f1 * (correct_clk - eph->toc) + eph->a_f2 * (correct_clk - eph->toc) * (correct_clk - eph->toc));
	}
	return eph->a_f0 + eph->a_f1 * (correct_clk - eph->toc) + eph->a_f2 * (correct_clk - eph->toc) * (correct_clk - eph->toc);
}

double TGD(EPHEMERIS* e, double f, int sys)
{
	switch (sys)
	{
	case SYS_GPS:
		return SQR(f / L1) * e->T_GD1;
		break;
	case SYS_BDS:
		if (f == B1)
			return e->T_GD1;
		if (f == B2)
			return e->T_GD2;
		if (f == B3)
			return 0.0;
	default:
		return -1;
		break;
	}
}

// ÐÇÀúÎ»ÖÃ
unsigned int SAT_POS_CAL(double t, EPHEMERIS* eph, XYZ* Sat_Pos, double& clk, double dt, int SYS)
{
	double n, delt_t, M, E, E0, V, u_, u, r, i, dt0, F;
	switch (SYS)
	{
	case SYS_GPS:
		n = sqrt(WGS84_GM) / (eph->sqrt_A * eph->sqrt_A * eph->sqrt_A) + eph->delt_n;
		dt0 = delt_t = t - eph->toe_tow;
		F = -2 * sqrt(WGS84_GM) / (velocity_c * velocity_c);
		break;
	case SYS_BDS:
		n = sqrt(CGCS2000_GM) / (eph->sqrt_A * eph->sqrt_A * eph->sqrt_A) + eph->delt_n;
		dt0 = delt_t = t - eph->toe_tow;
		F = -2 * sqrt(CGCS2000_GM) / (velocity_c * velocity_c);
		break;
	default:
		break;
	}

	double dtr = 0;
	while (abs(delt_t) > 302400)
	{
		if (delt_t > 302400)
		{
			delt_t -= 604800;
		}
		else if (delt_t < -302400)
		{
			delt_t += 604800;
		}
	}

	for (int i = 0; i < 10; i++)
	{
		M = eph->M0 + n * delt_t;
		E = M;
		E0 = M;
		int count = 0;
		do
		{
			E0 = E;
			E = M + eph->e * sin(E0);
			count++;
		} while (abs(E - E0) > 0.000000000001 && count < 10);
		dtr = F * eph->e * eph->sqrt_A * sin(E);
		delt_t = dt0 - dtr;
	}
	clk += dtr;
	dt += (dtr + clk);
	V = atan2((sqrt(1 - eph->e * eph->e) * sin(E)), (cos(E) - eph->e));
	u_ = eph->omiga + V;
	u = u_ + eph->Cuc * cos(2 * u_) + eph->Cus * sin(2 * u_);
	r = eph->sqrt_A * eph->sqrt_A * (1 - eph->e * cos(E)) + eph->Crc * cos(2 * u_) + eph->Crs * sin(2 * u_);
	i = eph->i0 + eph->Cic * cos(2 * u_) + eph->Cis * sin(2 * u_) + eph->dot_i * delt_t;
	double x = r * cos(u);
	double y = r * sin(u);
	double z = 0;
	double L;
	double x0, y0, z0;
	switch (SYS)
	{
	case SYS_GPS:
		L = eph->Omiga0 + (eph->dot_Omiga - omiga_earth) * t - eph->dot_Omiga * eph->toe_tow;
		Sat_Pos->X = (x * cos(L) - y * cos(i) * sin(L)) + omiga_earth * dt * (x * sin(L) + y * cos(i) * cos(L));
		Sat_Pos->Y = x * sin(L) + y * cos(i) * cos(L) - omiga_earth * dt * (x * cos(L) - y * cos(i) * sin(L));
		Sat_Pos->Z = y * sin(i);
		break;
	case SYS_BDS:
		if (fabs(eph->i0 - 0.0873) < 0.1 && fabs(eph->sqrt_A - 6493) < 1) // Í¨¹ý¹ìµÀÇã½ÇºÍ¹ìµÀ¸ùÊýÅÐ¶ÏÊÇ·ñÎªGEOÎÀÐÇ  i: 5/deg sqrt_A: 6493/sqrt_meter
		{
			L = eph->Omiga0 + eph->dot_Omiga * delt_t - omiga_earth * eph->toe_tow;
			x0 = x * cos(L) - y * cos(i) * sin(L);
			y0 = x * sin(L) + y * cos(i) * cos(L);
			z0 = y * sin(i);
			MatrixXd P_GK(3, 1);
			MatrixXd R_Z(3, 3);
			MatrixXd R_X(3, 3);
			MatrixXd P(3, 1);
			P_GK << x0,
				y0,
				z0;
			R_X << 1, 0, 0,
				0, cos(-5 * Pi / 180), sin(-5 * Pi / 180),
				0, -sin(-5 * Pi / 180), cos(-5 * Pi / 180);
			R_Z << cos(omiga_earth * delt_t), sin(omiga_earth * delt_t), 0,
				-sin(omiga_earth * delt_t), cos(omiga_earth * delt_t), 0,
				0, 0, 1;
			P = R_Z * R_X * P_GK;
			Sat_Pos->X = cos(omiga_earth * dt) * P(0, 0) + sin(omiga_earth * dt) * P(1, 0);
			Sat_Pos->Y = cos(omiga_earth * dt) * P(1, 0) - sin(omiga_earth * dt) * P(0, 0);
			Sat_Pos->Z = P(2, 0);
		}
		else
		{
			L = eph->Omiga0 + (eph->dot_Omiga - omiga_earth) * t - eph->dot_Omiga * eph->toe_tow;
			Sat_Pos->X = cos(omiga_earth * dt) * (x * cos(L) - y * cos(i) * sin(L)) + sin(omiga_earth * dt) * (x * sin(L) + y * cos(i) * cos(L));
			Sat_Pos->Y = cos(omiga_earth * dt) * (x * sin(L) + y * cos(i) * cos(L)) - sin(omiga_earth * dt) * (x * cos(L) - y * cos(i) * sin(L));
			Sat_Pos->Z = y * sin(i);
		}
		break;
	default:
		break;
	}

	return 1;
}

// ÎÀÐÇ¸ß¶È½Ç¼ÆËã
double Ele_Angle(XYZ SatPos, XYZ RcvPos, int sys)
{
	XYZ Satenu = XYZ2ENU(RcvPos, SatPos, sys);

	return asin(Satenu.Z / Len(&Satenu));
}

// Hopefiled¶ÔÁ÷²ã¸ÄÕý(m)
double Hopefield(double E, double H)
{
	double Ts = T0 - 0.0065 * (H - H0);
	double hd = 40136 + 148.72 * (Ts - 273.16);
	double Ps = P0 * pow(1 - 0.000026 * (H - H0), 5.225);
	double RH = RH0 * exp(-0.0006396 * (H - H0));
	double es = RH * exp(-37.2465 + 0.213166 * Ts - 0.000256908 * Ts * Ts);
	double md = sin(degree2rad(sqrt(E * E + 6.25)));
	double mw = sin(degree2rad(sqrt(E * E + 2.25)));
	double hw = 11000;
	double ZHD = 155.2 * 1e-7 * Ps * (hd - H) / Ts;
	double ZWD = 155.2 * 1e-7 * 4810 * es * (hw - H) / (Ts * Ts);
	return ZHD / md + ZWD / mw;
}

double Hopefield(XYZ SatPos, XYZ RcvPos, int sys)
{
	if (SatPos.X == 0 && SatPos.Y == 0 && SatPos.Z == 0)
	{
		return 0;
	}
	double E = rad2degree(Ele_Angle(SatPos, RcvPos, sys));
	BLH Rcvblh;
	double H = 0;
	switch (sys)
	{
	case SYS_GPS:
		Rcvblh = XYZ2BLH(RcvPos, WGS84_e2, WGS84_a);
		H = Rcvblh.Height;
		break;
	case SYS_BDS:
		Rcvblh = XYZ2BLH(RcvPos, CGCS2000_e2, CGCS2000_a);
		H = Rcvblh.Height;
		break;
	default:
		break;
	}
	if (H < 20e3 && H > -100)
	{
		double Ts = T0 - 0.0065 * (H - H0);
		double hd = 40136 + 148.72 * (T0 - 273.16);
		double Ps = P0 * pow(1 - 0.000026 * (H - H0), 5.225);
		double RH = RH0 * exp(-0.0006396 * (H - H0));
		double es = RH * exp(-37.2465 + 0.213166 * Ts - 0.000256908 * Ts * Ts);
		double md = sin(degree2rad(sqrt(E * E + 6.25)));
		double mw = sin(degree2rad(sqrt(E * E + 2.25)));
		double hw = 11000;
		double ZHD = 155.2 * 1e-7 * Ps * (hd - H) / Ts;
		double ZWD = 155.2 * 1e-7 * 4810 * es * (hw - H) / (Ts * Ts);
		return ZHD / md + ZWD / mw;
	}
	else
		return 0;
}

double Klobuchar(XYZ RcvPos, double E, double A, double alpha[4], double beta[4], double UT, double code, int sys)
{
	if (!(alpha[0] * alpha[1] * alpha[2] * alpha[3] * beta[0] * beta[1] * beta[2] * beta[3]))
		return -1;
	BLH RcvBLH;
	double T_g = 0;
	double EA = 0;
	double B_IPP = 0;
	double L_IPP = 0;
	double B_m = 0;
	double t = 0;
	double A_I = 0;
	double P_I = 0;
	double Phase_I = 0;
	double F = 0;
	switch (sys)
	{
	case SYS_GPS:
		RcvBLH = XYZ2BLH(RcvPos, WGS84_e2, WGS84_a);
		EA = 0.0137 / (E + 0.11) - 0.022;
		B_IPP = RcvBLH.Lat + EA * cos(A);
		if (B_IPP < -0.416)
			B_IPP = -0.416;
		if (B_IPP > 0.416)
			B_IPP = 0.416;
		L_IPP = RcvBLH.Lon + EA * sin(A) / cos(B_IPP);
		B_m = B_IPP + 0.064 * cos(L_IPP - 1.617);
		t = 43200 * L_IPP + UT;
		while (t > 86400 || t < 0)
		{
			if (t > 86400)
				t -= 86400;
			if (t < 0)
				t += 86400;
		}
		A_I = 0;
		P_I = 0;
		for (int i = 0; i < 4; i++)
		{
			A_I += alpha[i] * pow(B_m, i);
			P_I += beta[i] * pow(B_m, i);
		}
		if (A_I < 0)
			A_I = 0;
		if (P_I < 72000)
			P_I = 72000;
		Phase_I = 2 * Pi * (t - 50400) / P_I;
		F = 1 + 16 * pow((0.53 - E), 3);
		if (abs(Phase_I) < 1.57)
		{
			T_g = F * (5e-9 + A_I * (1 - pow(Phase_I, 2) / 2 + pow(Phase_I, 4) / 24));
		}
		else
		{
			T_g = F * 5e-9;
		}
		return pow(L1 / code, 2) * T_g;
		break;
	case SYS_BDS:
		RcvBLH = XYZ2BLH(RcvPos, CGCS2000_e2, CGCS2000_a);
		EA = Pi / 2 - E - asin(cos(E) * 6378 / (6378 + 375));
		B_IPP = asin(sin(RcvBLH.Lat) * cos(EA) + cos(RcvBLH.Lat) * sin(EA) * cos(A));
		L_IPP = RcvBLH.Lon + asin(sin(EA) * sin(A) / cos(B_IPP));
		t = UT + L_IPP * 43200 / Pi;
		A_I = 0;
		P_I = 0;
		for (int i = 0; i < 4; i++)
		{
			A_I += alpha[i] * pow(B_IPP / Pi, i);
			P_I += beta[i] * pow(B_IPP / Pi, i);
		}
		if (A_I < 0)
			A_I = 0;
		if (P_I < 72000)
			P_I = 72000;
		if (P_I > 172800)
			P_I = 172800;
		if (abs(t - 50400) < P_I / 4)
		{
			T_g = 5e-9 + A_I * cos(2 * Pi * (t - 50400) / P_I);
		}
		else
		{
			5e-9;
		}
		return T_g / sqrt(1 - pow(cos(E) * 6378 / (6378 + 375), 2));
		break;
	default:
		return 0;
		break;
	}
}

#pragma region vector_EPOCH
unsigned int setup_Pos(GPSTIME* OBS_TIME, MatrixXd Pos, vector<Satellate*> Sates, EPOCH* eph, bool first_flag, double f1, double f2, MatrixXd* B_Pos, MatrixXd* l_Pos, MatrixXd* P_Pos, string* sate_used)
{
	XYZ sate_pos;
	double clk = 0;
	XYZ RcvPos = get_XYZ(Pos.block(0, 0, 3, 1));
	double dt_Rcv = Pos(3, 0);
	int ROWS = 0;
	MatrixXd x_Pos = MatrixXd::Zero(4, 1);
	MatrixXd B_pos_new = MatrixXd::Zero(1, 4);
	MatrixXd l_pos_new = MatrixXd::Zero(1, 1);
	for (int i = 0; i < Sates.size(); i++)
	{
		if (Sates[i]->Phase_NUM < 2)
			continue;
		if (Sates[i]->Outlier)
			continue;
		int prn = Sates[i]->PRN;
		if (eph[prn - 1].num == 0)
			continue;
		double IF = 0;
		int Index1 = 0;
		int Index2 = 0;
		while (CODE2FREQ(decode_SYN(Sates[i]->SYS, Sates[i]->SYG_TYPE[Index1])) != f1 && Index1 < MAXNUM)
			Index1++;
		while (CODE2FREQ(decode_SYN(Sates[i]->SYS, Sates[i]->SYG_TYPE[Index2])) != f2 && Index2 < MAXNUM)
			Index2++;
		if (Index1 == MAXNUM || Index2 == MAXNUM)
			continue;
		if (!(Sates[i]->LOCK_PSE[Index1] && Sates[i]->LOCK_PSE[Index2] && Sates[i]->LOCK_PHA[Index1] && Sates[i]->LOCK_PHA[Index2]))
			continue;

		// ¼ÆËãÎÀÐÇÎ»ÖÃ¡¢ÖÓ²î
		double ts = OBS_TIME->SecOfWeek - Sates[i]->PSERA[0] / velocity_c;
		if (Sates[i]->SYS == SYS_BDS)
			ts -= 14;
		double dt = abs(ts - eph[prn - 1].epoch[0]->toe_tow);
		int index = 0;
		for (int j = 1; j < eph[prn - 1].num; j++)
		{
			if (abs(ts - eph[prn - 1].epoch[j]->toe_tow))
			{
				dt = abs(ts - eph[prn - 1].epoch[j]->toe_tow);
				index = j;
			}
		}

		for (int j = 0; j < 3; j++)
		{
			clk = CORRECT_CLK(ts - clk, eph[prn - 1].epoch[index]);
			SAT_POS_CAL(ts - clk, eph[prn - 1].epoch[index], &sate_pos, clk, Sates[i]->PSERA[0] / velocity_c, Sates[i]->SYS);
		}
		if (!first_flag)
		{
			if (Ele_Angle(sate_pos, RcvPos, Sates[i]->SYS) < degree2rad(10))
				continue;
		}
		if (Sates[i]->SYS == SYS_GPS)
			IF = SQR(f1) * Sates[i]->PSERA[Index1] / (SQR(f1) - SQR(f2)) - SQR(f2) * Sates[i]->PSERA[Index2] / (SQR(f1) - SQR(f2));
		else if (Sates[i]->SYS == SYS_BDS)
		{
			double k_1_3 = SQR(f1 / f2);
			IF = (Sates[i]->PSERA[Index2] - k_1_3 * Sates[i]->PSERA[Index1]) / (1 - k_1_3) + velocity_c * k_1_3 * eph[prn - 1].epoch[index]->T_GD1 / (1 - k_1_3);
		}

		double len = sqrt(SQR(RcvPos.X - sate_pos.X) + SQR(RcvPos.Y - sate_pos.Y) + SQR(RcvPos.Z - sate_pos.Z));
		double w_pos = IF - (len + dt_Rcv - velocity_c * clk + Hopefield(sate_pos, RcvPos, Sates[i]->SYS));
		if (!first_flag)
		{
			if (abs(w_pos) > 10)
				continue;
		}
		double l = (RcvPos.X - sate_pos.X) / len;
		double m = (RcvPos.Y - sate_pos.Y) / len;
		double n = (RcvPos.Z - sate_pos.Z) / len;
		B_pos_new(0, 0) = l;
		B_pos_new(0, 1) = m;
		B_pos_new(0, 2) = n;
		B_pos_new(0, 3) = 1;
		l_pos_new(0, 0) = w_pos;
		B_Pos->conservativeResize(B_Pos->rows() + 1, B_Pos->cols());
		B_Pos->bottomRows(1) = B_pos_new;
		l_Pos->conservativeResize(l_Pos->rows() + 1, l_Pos->cols());
		l_Pos->bottomRows(1) = l_pos_new;
		switch (Sates[i]->SYS)
		{
		case SYS_GPS:
			*sate_used += "G" + to_string(Sates[i]->PRN);
			break;
		case SYS_BDS:
			*sate_used += "C" + to_string(Sates[i]->PRN);
			break;
		default:
			break;
		}
		ROWS++;
	}
	*P_Pos = MatrixXd::Identity(ROWS, ROWS);
	return ROWS;
}

unsigned int setup_Pos(GPSTIME* OBS_TIME, MatrixXd Pos, vector<Satellate*> Sates, EPOCH* eph, bool first_flag, double f, MatrixXd* B_Pos, MatrixXd* l_Pos, MatrixXd* P_Pos, string* sate_used)
{
	XYZ sate_pos;
	double clk = 0;
	XYZ RcvPos = get_XYZ(Pos.block(0, 0, 3, 1));
	double dt_Rcv = Pos(3, 0);
	int ROWS = 0;
	MatrixXd x_Pos = MatrixXd::Zero(4, 1);
	MatrixXd B_pos_new = MatrixXd::Zero(1, 4);
	MatrixXd l_pos_new = MatrixXd::Zero(1, 1);
	for (int i = 0; i < Sates.size(); i++)
	{
		if (Sates[i]->Phase_NUM < 2)
			continue;
		if (Sates[i]->Outlier)
			continue;
		int prn = Sates[i]->PRN;
		if (eph[prn - 1].num == 0)
			continue;
		int Index = 0;
		while (CODE2FREQ(decode_SYN(Sates[i]->SYS, Sates[i]->SYG_TYPE[Index])) != f && Index < MAXNUM)
			Index++;
		if (Index == MAXNUM)
			continue;
		if (!(Sates[i]->LOCK_PSE[Index] && Sates[i]->LOCK_PHA[Index]))
			continue;

		// ¼ÆËãÎÀÐÇÎ»ÖÃ¡¢ÖÓ²î
		double ts = OBS_TIME->SecOfWeek - Sates[i]->PSERA[0] / velocity_c;
		if (Sates[i]->SYS == SYS_BDS)
			ts -= 14;
		double dt = abs(ts - eph[prn - 1].epoch[0]->toe_tow);
		int index = 0;
		for (int j = 1; j < eph[prn - 1].num; j++)
		{
			if (abs(ts - eph[prn - 1].epoch[j]->toe_tow))
			{
				dt = abs(ts - eph[prn - 1].epoch[j]->toe_tow);
				index = j;
			}
		}
		double tgd = 0;
		tgd = TGD(eph[prn - 1].epoch[index], f, Sates[i]->SYS);

		for (int j = 0; j < 3; j++)
		{
			clk = CORRECT_CLK(ts - clk, eph[prn - 1].epoch[index]);
			SAT_POS_CAL(ts - clk + tgd, eph[prn - 1].epoch[index], &sate_pos, clk, Sates[i]->PSERA[0] / velocity_c, Sates[i]->SYS);
		}
		if (!first_flag)
		{
			if (Ele_Angle(sate_pos, RcvPos, Sates[i]->SYS) < degree2rad(10))
				continue;
		}

		double len = sqrt(SQR(RcvPos.X - sate_pos.X) + SQR(RcvPos.Y - sate_pos.Y) + SQR(RcvPos.Z - sate_pos.Z));
		double w_pos = Sates[i]->PSERA[Index] - (len + dt_Rcv - velocity_c * clk + Hopefield(sate_pos, RcvPos, Sates[i]->SYS));
		if (!first_flag)
		{
			if (abs(w_pos) > 10)
				continue;
		}
		double l = (RcvPos.X - sate_pos.X) / len;
		double m = (RcvPos.Y - sate_pos.Y) / len;
		double n = (RcvPos.Z - sate_pos.Z) / len;
		B_pos_new(0, 0) = l;
		B_pos_new(0, 1) = m;
		B_pos_new(0, 2) = n;
		B_pos_new(0, 3) = 1;
		l_pos_new(0, 0) = w_pos;
		B_Pos->conservativeResize(B_Pos->rows() + 1, B_Pos->cols());
		B_Pos->bottomRows(1) = B_pos_new;
		l_Pos->conservativeResize(l_Pos->rows() + 1, l_Pos->cols());
		l_Pos->bottomRows(1) = l_pos_new;
		switch (Sates[i]->SYS)
		{
		case SYS_GPS:
			*sate_used += "G" + to_string(Sates[i]->PRN);
			break;
		case SYS_BDS:
			*sate_used += "C" + to_string(Sates[i]->PRN);
			break;
		default:
			break;
		}
		ROWS++;
	}
	*P_Pos = MatrixXd::Identity(ROWS, ROWS);
	return ROWS;
}

unsigned int setup_Vel(GPSTIME* OBS_TIME, MatrixXd Pos, vector<Satellate*> Sates, EPOCH* eph, MatrixXd* B_Vel, MatrixXd* l_Vel, MatrixXd* P_Vel)
{
	XYZ* sate_pos0 = new XYZ();
	XYZ* sate_pos1 = new XYZ();
	double clk0 = 0;
	double clk1 = 0;
	double velocity[4] = { 0, 0, 0, 0 };
	XYZ* RcvPos = new XYZ(Pos(0, 0), Pos(1, 0), Pos(2, 0));
	MatrixXd x_vel = MatrixXd::Zero(4, 1);
	MatrixXd B_vel_new = MatrixXd::Zero(1, 4);
	MatrixXd l_vel_new = MatrixXd::Zero(1, 1);
	int ROWS = 0;
	for (int i = 0; i < Sates.size(); i++)
	{
		if (Sates[i]->Phase_NUM < 2)
			continue;
		if (Sates[i]->Outlier)
			continue;
		int prn = Sates[i]->PRN;
		if (eph[prn - 1].num == 0)
			continue;
		double ts = OBS_TIME->SecOfWeek - Sates[i]->PSERA[0] / velocity_c;
		if (Sates[i]->SYS == SYS_BDS)
			ts -= 14;
		double dt = abs(ts - eph[prn - 1].epoch[0]->toe_tow);
		int index = 0;
		for (int j = 1; j < eph[prn - 1].num; j++)
		{
			if (abs(ts - eph[prn - 1].epoch[j]->toe_tow) < dt)
			{
				dt = abs(ts - eph[prn - 1].epoch[j]->toe_tow);
				index = j;
			}
		}

		for (int j = 0; j < 3; j++)
		{
			clk0 = CORRECT_CLK(ts - clk0, eph[prn - 1].epoch[index]);
			clk1 = CORRECT_CLK(ts + 1e-3 - clk1, eph[prn - 1].epoch[index]);
			SAT_POS_CAL(ts - clk0, eph[prn - 1].epoch[index], sate_pos0, clk0, Sates[i]->PSERA[0] / velocity_c, Sates[i]->SYS);
			SAT_POS_CAL(ts + 1e-3 - clk1, eph[prn - 1].epoch[index], sate_pos1, clk1, Sates[i]->PSERA[0] / velocity_c, Sates[i]->SYS);
		}
		velocity[0] = (sate_pos1->X - sate_pos0->X) / 1e-3;
		velocity[1] = (sate_pos1->Y - sate_pos0->Y) / 1e-3;
		velocity[2] = (sate_pos1->Z - sate_pos0->Z) / 1e-3;
		velocity[3] = (clk1 - clk0) * velocity_c / 1e-3;
		double f1 = CODE2FREQ(decode_SYN(Sates[i]->SYS, Sates[i]->SYG_TYPE[0]));
		double lamda = (1e-6 * velocity_c / f1);
		double len = sqrt(SQR(RcvPos->X - sate_pos0->X) + SQR(RcvPos->Y - sate_pos0->Y) + SQR(RcvPos->Z - sate_pos0->Z));
		double l = (RcvPos->X - sate_pos0->X) / len;
		double m = (RcvPos->Y - sate_pos0->Y) / len;
		double n = (RcvPos->Z - sate_pos0->Z) / len;
		double v0 = ((sate_pos0->X - RcvPos->X) * velocity[0] + (sate_pos0->Y - RcvPos->Y) * velocity[1] + (sate_pos0->Z - RcvPos->Z) * velocity[2]) / len;
		double w_Vel = -lamda * Sates[i]->DOPPLER[0] - (v0 - velocity[3]);
		if (abs(w_Vel) > 10)
			continue;
		B_vel_new(0, 0) = l;
		B_vel_new(0, 1) = m;
		B_vel_new(0, 2) = n;
		B_vel_new(0, 3) = 1;
		l_vel_new(0, 0) = w_Vel;
		B_Vel->conservativeResize(B_Vel->rows() + 1, B_Vel->cols());
		B_Vel->bottomRows(1) = B_vel_new;
		l_Vel->conservativeResize(l_Vel->rows() + 1, l_Vel->cols());
		l_Vel->bottomRows(1) = l_vel_new;
		ROWS++;
	}
	*P_Vel = MatrixXd::Identity(ROWS, ROWS);
	delete sate_pos0;
	delete sate_pos1;
	delete RcvPos;
	return ROWS;
}

#pragma endregion

unsigned int setup_Pos(GPSTIME* OBS_TIME, MatrixXd Pos, vector<Satellate*> Sates, EPHEMERIS** eph, bool first_flag, double f1, double f2, MatrixXd* B_Pos, MatrixXd* l_Pos, MatrixXd* P_Pos, string* sate_used)
{
	XYZ sate_pos;
	double clk = 0;
	XYZ RcvPos = get_XYZ(Pos.block(0, 0, 3, 1));
	double dt_Rcv = Pos(3, 0);
	int ROWS = 0;
	MatrixXd x_Pos = MatrixXd::Zero(4, 1);
	MatrixXd B_pos_new = MatrixXd::Zero(1, 4);
	MatrixXd l_pos_new = MatrixXd::Zero(1, 1);
	for (int i = 0; i < Sates.size(); i++)
	{
		if (Sates[i]->Phase_NUM < 2)
			continue;
		if (Sates[i]->Outlier)
			continue;
		int prn = Sates[i]->PRN;
		if (eph[prn - 1]->PRN != prn)
			continue;
		double IF = 0;
		int Index1 = 0;
		int Index2 = 0;
		while (CODE2FREQ(decode_SYN(Sates[i]->SYS, Sates[i]->SYG_TYPE[Index1])) != f1 && Index1 < MAXNUM)
			Index1++;
		while (CODE2FREQ(decode_SYN(Sates[i]->SYS, Sates[i]->SYG_TYPE[Index2])) != f2 && Index2 < MAXNUM)
			Index2++;
		if (Index1 == MAXNUM || Index2 == MAXNUM)
			continue;
		if (!(Sates[i]->LOCK_PSE[Index1] && Sates[i]->LOCK_PSE[Index2] && Sates[i]->LOCK_PHA[Index1] && Sates[i]->LOCK_PHA[Index2]))
			continue;

		// ¼ÆËãÎÀÐÇÎ»ÖÃ¡¢ÖÓ²î
		double ts = OBS_TIME->SecOfWeek - Sates[i]->PSERA[0] / velocity_c;
		double dt = abs(ts - eph[prn - 1]->toe_tow + (OBS_TIME->Week - eph[prn - 1]->toe_wn) * 604800);
		if (dt > 14400)
			continue;

		for (int j = 0; j < 3; j++)
		{
			clk = CORRECT_CLK(ts - clk, eph[prn - 1]);
			SAT_POS_CAL(ts - clk, eph[prn - 1], &sate_pos, clk, Sates[i]->PSERA[0] / velocity_c, Sates[i]->SYS);
		}
		if (!first_flag)
		{
			if (Ele_Angle(sate_pos, RcvPos, Sates[i]->SYS) < degree2rad(10))
				continue;
		}
		if (Sates[i]->SYS == SYS_GPS)
			IF = SQR(f1) * Sates[i]->PSERA[Index1] / (SQR(f1) - SQR(f2)) - SQR(f2) * Sates[i]->PSERA[Index2] / (SQR(f1) - SQR(f2));
		else if (Sates[i]->SYS == SYS_BDS)
		{
			double k_1_3 = SQR(f1 / f2);
			IF = (Sates[i]->PSERA[Index2] - k_1_3 * Sates[i]->PSERA[Index1]) / (1 - k_1_3) + velocity_c * k_1_3 * eph[prn - 1]->T_GD1 / (1 - k_1_3);
		}

		double len = sqrt(SQR(RcvPos.X - sate_pos.X) + SQR(RcvPos.Y - sate_pos.Y) + SQR(RcvPos.Z - sate_pos.Z));
		double w_pos = IF - (len + dt_Rcv - velocity_c * clk + Hopefield(sate_pos, RcvPos, Sates[i]->SYS));
		if (!first_flag)
		{
			if (abs(w_pos) > 10)
				continue;
		}
		double l = (RcvPos.X - sate_pos.X) / len;
		double m = (RcvPos.Y - sate_pos.Y) / len;
		double n = (RcvPos.Z - sate_pos.Z) / len;
		B_pos_new(0, 0) = l;
		B_pos_new(0, 1) = m;
		B_pos_new(0, 2) = n;
		B_pos_new(0, 3) = 1;
		l_pos_new(0, 0) = w_pos;
		B_Pos->conservativeResize(B_Pos->rows() + 1, B_Pos->cols());
		B_Pos->bottomRows(1) = B_pos_new;
		l_Pos->conservativeResize(l_Pos->rows() + 1, l_Pos->cols());
		l_Pos->bottomRows(1) = l_pos_new;
		switch (Sates[i]->SYS)
		{
		case SYS_GPS:
			*sate_used += "G" + to_string(Sates[i]->PRN);
			break;
		case SYS_BDS:
			*sate_used += "C" + to_string(Sates[i]->PRN);
			break;
		default:
			break;
		}
		ROWS++;
	}
	*P_Pos = MatrixXd::Identity(ROWS, ROWS);
	return ROWS;
}

unsigned int setup_Pos(GPSTIME* OBS_TIME, MatrixXd Pos, vector<Satellate*> Sates, EPHEMERIS** eph, bool first_flag, double f, MatrixXd* B_Pos, MatrixXd* l_Pos, MatrixXd* P_Pos, string* sate_used)
{
	XYZ sate_pos;
	double clk = 0;
	XYZ RcvPos = get_XYZ(Pos.block(0, 0, 3, 1));
	double dt_Rcv = Pos(3, 0);
	int ROWS = 0;
	MatrixXd x_Pos = MatrixXd::Zero(4, 1);
	MatrixXd B_pos_new = MatrixXd::Zero(1, 4);
	MatrixXd l_pos_new = MatrixXd::Zero(1, 1);
	for (int i = 0; i < Sates.size(); i++)
	{
		if (Sates[i]->Phase_NUM < 2)
			continue;
		if (Sates[i]->Outlier)
			continue;
		int prn = Sates[i]->PRN;
		if (eph[prn - 1]->PRN != prn)
			continue;
		int Index = 0;
		while (CODE2FREQ(decode_SYN(Sates[i]->SYS, Sates[i]->SYG_TYPE[Index])) != f && Index < MAXNUM)
			Index++;
		if (Index == MAXNUM)
			continue;
		if (!(Sates[i]->LOCK_PSE[Index] && Sates[i]->LOCK_PHA[Index]))
			continue;

		// ¼ÆËãÎÀÐÇÎ»ÖÃ¡¢ÖÓ²î
		double ts = OBS_TIME->SecOfWeek - Sates[i]->PSERA[0] / velocity_c;
		double dt = abs(ts - eph[prn - 1]->toe_tow + (OBS_TIME->Week - eph[prn - 1]->toe_wn) * 604800);
		if (dt > 14400)
			continue;

		double tgd = 0;
		tgd = TGD(eph[prn - 1], f, Sates[i]->SYS);

		for (int j = 0; j < 3; j++)
		{
			clk = CORRECT_CLK(ts - clk, eph[prn - 1]);
			SAT_POS_CAL(ts - clk + tgd, eph[prn - 1], &sate_pos, clk, Sates[i]->PSERA[0] / velocity_c, Sates[i]->SYS);
		}
		if (!first_flag)
		{
			if (Ele_Angle(sate_pos, RcvPos, Sates[i]->SYS) < degree2rad(10))
				continue;
		}

		double len = sqrt(SQR(RcvPos.X - sate_pos.X) + SQR(RcvPos.Y - sate_pos.Y) + SQR(RcvPos.Z - sate_pos.Z));
		double w_pos = Sates[i]->PSERA[Index] - (len + dt_Rcv - velocity_c * clk + Hopefield(sate_pos, RcvPos, Sates[i]->SYS));
		if (!first_flag)
		{
			if (abs(w_pos) > 10)
				continue;
		}
		double l = (RcvPos.X - sate_pos.X) / len;
		double m = (RcvPos.Y - sate_pos.Y) / len;
		double n = (RcvPos.Z - sate_pos.Z) / len;
		B_pos_new(0, 0) = l;
		B_pos_new(0, 1) = m;
		B_pos_new(0, 2) = n;
		B_pos_new(0, 3) = 1;
		l_pos_new(0, 0) = w_pos;
		B_Pos->conservativeResize(B_Pos->rows() + 1, B_Pos->cols());
		B_Pos->bottomRows(1) = B_pos_new;
		l_Pos->conservativeResize(l_Pos->rows() + 1, l_Pos->cols());
		l_Pos->bottomRows(1) = l_pos_new;
		switch (Sates[i]->SYS)
		{
		case SYS_GPS:
			*sate_used += "G" + to_string(Sates[i]->PRN);
			break;
		case SYS_BDS:
			*sate_used += "C" + to_string(Sates[i]->PRN);
			break;
		default:
			break;
		}
		ROWS++;
	}
	*P_Pos = MatrixXd::Identity(ROWS, ROWS);
	return ROWS;
}

unsigned int setup_Vel(GPSTIME* OBS_TIME, MatrixXd Pos, vector<Satellate*> Sates, EPHEMERIS** eph, MatrixXd* B_Vel, MatrixXd* l_Vel, MatrixXd* P_Vel)
{
	XYZ* sate_pos0 = new XYZ();
	XYZ* sate_pos1 = new XYZ();
	double clk0 = 0;
	double clk1 = 0;
	double velocity[4] = { 0, 0, 0, 0 };
	XYZ* RcvPos = new XYZ(Pos(0, 0), Pos(1, 0), Pos(2, 0));
	MatrixXd x_vel = MatrixXd::Zero(4, 1);
	MatrixXd B_vel_new = MatrixXd::Zero(1, 4);
	MatrixXd l_vel_new = MatrixXd::Zero(1, 1);
	int ROWS = 0;
	for (int i = 0; i < Sates.size(); i++)
	{
		if (Sates[i]->Phase_NUM < 2)
			continue;
		if (Sates[i]->Outlier)
			continue;
		int prn = Sates[i]->PRN;
		if (eph[prn - 1]->PRN != prn)
			continue;
		double ts = OBS_TIME->SecOfWeek - Sates[i]->PSERA[0] / velocity_c;

		double dt = abs(ts - eph[prn - 1]->toe_tow + (OBS_TIME->Week - eph[prn - 1]->toe_wn) * 604800);
		if (dt > 14400)
			continue;

		for (int j = 0; j < 3; j++)
		{
			clk0 = CORRECT_CLK(ts - clk0, eph[prn - 1]);
			clk1 = CORRECT_CLK(ts + 1e-3 - clk1, eph[prn - 1]);
			SAT_POS_CAL(ts - clk0, eph[prn - 1], sate_pos0, clk0, Sates[i]->PSERA[0] / velocity_c, Sates[i]->SYS);
			SAT_POS_CAL(ts + 1e-3 - clk1, eph[prn - 1], sate_pos1, clk1, Sates[i]->PSERA[0] / velocity_c, Sates[i]->SYS);
		}
		velocity[0] = (sate_pos1->X - sate_pos0->X) / 1e-3;
		velocity[1] = (sate_pos1->Y - sate_pos0->Y) / 1e-3;
		velocity[2] = (sate_pos1->Z - sate_pos0->Z) / 1e-3;
		velocity[3] = (clk1 - clk0) * velocity_c / 1e-3;
		double f1 = CODE2FREQ(decode_SYN(Sates[i]->SYS, Sates[i]->SYG_TYPE[0]));
		double lamda = (1e-6 * velocity_c / f1);
		double len = sqrt(SQR(RcvPos->X - sate_pos0->X) + SQR(RcvPos->Y - sate_pos0->Y) + SQR(RcvPos->Z - sate_pos0->Z));
		double l = (RcvPos->X - sate_pos0->X) / len;
		double m = (RcvPos->Y - sate_pos0->Y) / len;
		double n = (RcvPos->Z - sate_pos0->Z) / len;
		double v0 = ((sate_pos0->X - RcvPos->X) * velocity[0] + (sate_pos0->Y - RcvPos->Y) * velocity[1] + (sate_pos0->Z - RcvPos->Z) * velocity[2]) / len;
		double w_Vel = -lamda * Sates[i]->DOPPLER[0] - (v0 - velocity[3]);
		if (abs(w_Vel) > 10)
			continue;
		B_vel_new(0, 0) = l;
		B_vel_new(0, 1) = m;
		B_vel_new(0, 2) = n;
		B_vel_new(0, 3) = 1;
		l_vel_new(0, 0) = w_Vel;
		B_Vel->conservativeResize(B_Vel->rows() + 1, B_Vel->cols());
		B_Vel->bottomRows(1) = B_vel_new;
		l_Vel->conservativeResize(l_Vel->rows() + 1, l_Vel->cols());
		l_Vel->bottomRows(1) = l_vel_new;
		ROWS++;
	}
	*P_Vel = MatrixXd::Identity(ROWS, ROWS);
	delete sate_pos0;
	delete sate_pos1;
	delete RcvPos;
	return ROWS;
}



unsigned int setup_LS(DATA_SET* data, Configure cfg, int sys)
{
	int ROWS = 0;
	vector<Satellate*> Sates;
	EPHEMERIS** eph;
	double f = 0;
	switch (sys)
	{
	case SYS_GPS:
		Sates = data->range->GPS_SATE;
		eph = data->GPS_eph;
		f = cfg.GPS_Cfg.f1;
		break;
	case SYS_BDS:
		Sates = data->range->BDS_SATE;
		eph = data->BDS_eph;
		f = cfg.BDS_Cfg.f1;
		break;
	default:
		return 0;
		break;
	}

	XYZ RcvPos = get_XYZ(data->temp_ref.block(0, 0, 3, 1));
	double dt_Rcv = data->temp_ref(3, 0);
	XYZ sate_pos, sate_pos1;
	double clk = 0;
	double clk1 = 0;
	double velocity[4] = { 0, 0, 0, 0 };
	MatrixXd B, l_Pos, P_Pos;
	MatrixXd l_Vel, P_Vel;
	B = MatrixXd::Zero(0, 4);
	l_Vel = l_Pos = MatrixXd::Zero(0, 1);
	MatrixXd B_pos_new = MatrixXd::Zero(1, 4);
	MatrixXd l_pos_new = MatrixXd::Zero(1, 1);
	MatrixXd l_vel_new = MatrixXd::Zero(1, 1);
	for (int i = 0; i < Sates.size(); i++)
	{
		int prn = Sates[i]->PRN;
		if (eph[prn - 1]->PRN != prn)
			continue;
		int Index = 0;
		while (CODE2FREQ(decode_SYN(Sates[i]->SYS, Sates[i]->SYG_TYPE[Index])) != f && Index < MAXNUM)
			Index++;
		if (Index == MAXNUM)
			return 0;

		double measure = get_measure(Sates[i], cfg, eph[prn - 1]);
		if (!measure)
			continue;
		// ¼ÆËãÎÀÐÇÎ»ÖÃ¡¢ÖÓ²î
		double ts = data->OBSTIME->SecOfWeek - Sates[i]->PSERA[0] / velocity_c;
		int week = data->OBSTIME->Week;
		if (sys == SYS_BDS)
		{
			ts -= 14;
			week -= 1356;
		}
		double dt = abs(ts - eph[prn - 1]->toe_tow + (week - eph[prn - 1]->toe_wn) * 604800);
		if (dt > 14400)
			continue;

		double tgd = 0;
		if (cfg.phase_num == 1)
			tgd = TGD(eph[prn - 1], f, Sates[i]->SYS);

		for (int j = 0; j < 3; j++)
		{
			clk = CORRECT_CLK(ts - clk, eph[prn - 1]);
			SAT_POS_CAL(ts - clk + tgd, eph[prn - 1], &sate_pos, clk, Sates[i]->PSERA[0] / velocity_c, Sates[i]->SYS);
			clk1 = CORRECT_CLK(ts - clk1 + 1e-3, eph[prn - 1]);
			SAT_POS_CAL(ts - clk1 + tgd + 1e-3, eph[prn - 1], &sate_pos1, clk1, Sates[i]->PSERA[0] / velocity_c, Sates[i]->SYS);
		}
		if (!data->LS_first)
		{
			if (Ele_Angle(sate_pos, RcvPos, Sates[i]->SYS) < degree2rad(10))
				continue;
		}
		velocity[0] = (sate_pos1.X - sate_pos.X) / 1e-3;
		velocity[1] = (sate_pos1.Y - sate_pos.Y) / 1e-3;
		velocity[2] = (sate_pos1.Z - sate_pos.Z) / 1e-3;
		velocity[3] = (clk1 - clk) * velocity_c / 1e-3;
		double lamda = (1e-6 * velocity_c / f);
		double len = sqrt(SQR(RcvPos.X - sate_pos.X) + SQR(RcvPos.Y - sate_pos.Y) + SQR(RcvPos.Z - sate_pos.Z));
		double w_pos = measure - (len + dt_Rcv - velocity_c * clk + Hopefield(sate_pos, RcvPos, Sates[i]->SYS));
		double v0 = ((sate_pos.X - RcvPos.X) * velocity[0] + (sate_pos.Y - RcvPos.Y) * velocity[1] + (sate_pos.Z - RcvPos.Z) * velocity[2]) / len;
		double w_vel = -lamda * Sates[i]->DOPPLER[Index] - (v0 - velocity[3]);
		if (!data->LS_first)
		{
			if (abs(w_pos) > 10)
				continue;
		}
		double l = (RcvPos.X - sate_pos.X) / len;
		double m = (RcvPos.Y - sate_pos.Y) / len;
		double n = (RcvPos.Z - sate_pos.Z) / len;
		B_pos_new(0, 0) = l;
		B_pos_new(0, 1) = m;
		B_pos_new(0, 2) = n;
		B_pos_new(0, 3) = 1;
		l_pos_new(0, 0) = w_pos;
		l_vel_new(0, 0) = w_vel;
		B.conservativeResize(B.rows() + 1, B.cols());
		B.bottomRows(1) = B_pos_new;
		l_Pos.conservativeResize(l_Pos.rows() + 1, l_Pos.cols());
		l_Pos.bottomRows(1) = l_pos_new;
		l_Vel.conservativeResize(l_Vel.rows() + 1, l_Vel.cols());
		l_Vel.bottomRows(1) = l_vel_new;
		switch (sys)
		{
		case SYS_GPS:
			*(data->LS_SATES) += "G" + to_string(Sates[i]->PRN);
			break;
		case SYS_BDS:
			*(data->LS_SATES) += "C" + to_string(Sates[i]->PRN);
			break;
		default:
			break;
		}
		ROWS++;
	}
	P_Pos = MatrixXd::Identity(ROWS, ROWS);
	P_Vel = MatrixXd::Identity(ROWS, ROWS);
	if (ROWS != 0)
	{
		data->LS_Pos->set_B_Pos(B);
		data->LS_Pos->set_l(l_Pos);
		data->LS_Pos->set_P(P_Pos);
		data->LS_Vel->set_B_Vel(B);
		data->LS_Vel->set_l(l_Vel);
		data->LS_Vel->set_P(P_Vel);
	}
	switch (sys)
	{
	case SYS_GPS:
		data->LS_GPS_num = ROWS;
		break;
	case SYS_BDS:
		data->LS_BDS_num = ROWS;
		break;
	default:
		break;
	}
	return ROWS;
}

unsigned int setup_KF(DATA_SET* data, Configure cfg, int sys)
{
	int ROWS = 0;
	vector<Satellate*> Sates;
	EPHEMERIS** eph;
	double f = 0;
	switch (sys)
	{
	case SYS_GPS:
		Sates = data->range->GPS_SATE;
		eph = data->GPS_eph;
		f = cfg.GPS_Cfg.f1;
		break;
	case SYS_BDS:
		Sates = data->range->BDS_SATE;
		eph = data->BDS_eph;
		f = cfg.BDS_Cfg.f1;
		break;
	default:
		return 0;
		break;
	}

	XYZ RcvPos = get_XYZ(data->temp_ref.block(0, 0, 3, 1));
	double dt_Rcv = data->temp_ref(3, 0);
	XYZ sate_pos, sate_pos1;
	double clk = 0;
	double clk1 = 0;
	double velocity[4] = { 0, 0, 0, 0 };
	MatrixXd B, l_Pos;
	MatrixXd l_Vel;
	B = MatrixXd::Zero(0, 4);
	l_Vel = l_Pos = MatrixXd::Zero(0, 1);
	MatrixXd B_pos_new = MatrixXd::Zero(1, 4);
	MatrixXd l_pos_new = MatrixXd::Zero(1, 1);
	MatrixXd l_vel_new = MatrixXd::Zero(1, 1);
	for (int i = 0; i < Sates.size(); i++)
	{
		int prn = Sates[i]->PRN;
		if (eph[prn - 1]->PRN != prn)
			continue;
		int Index = 0;
		while (CODE2FREQ(decode_SYN(Sates[i]->SYS, Sates[i]->SYG_TYPE[Index])) != f && Index < MAXNUM)
			Index++;
		if (Index == MAXNUM)
			return 0;

		double measure = get_measure(Sates[i], cfg, eph[prn - 1]);
		if (!measure)
			continue;
		// ¼ÆËãÎÀÐÇÎ»ÖÃ¡¢ÖÓ²î
		double ts = data->OBSTIME->SecOfWeek - Sates[i]->PSERA[0] / velocity_c;
		int week = data->OBSTIME->Week;
		if (sys == SYS_BDS)
		{
			ts -= 14;
			week -= 1356;
		}
		double dt = abs(ts - eph[prn - 1]->toe_tow + (week - eph[prn - 1]->toe_wn) * 604800);
		if (dt > 14400)
			continue;

		double tgd = 0;
		if (cfg.phase_num == 1)
			tgd = TGD(eph[prn - 1], f, Sates[i]->SYS);

		for (int j = 0; j < 3; j++)
		{
			clk = CORRECT_CLK(ts - clk, eph[prn - 1]);
			SAT_POS_CAL(ts - clk + tgd, eph[prn - 1], &sate_pos, clk, Sates[i]->PSERA[0] / velocity_c, Sates[i]->SYS);
			clk1 = CORRECT_CLK(ts - clk1 + 1e-3, eph[prn - 1]);
			SAT_POS_CAL(ts - clk1 + tgd + 1e-3, eph[prn - 1], &sate_pos1, clk1, Sates[i]->PSERA[0] / velocity_c, Sates[i]->SYS);
		}
		if (Ele_Angle(sate_pos, RcvPos, Sates[i]->SYS) < degree2rad(10))
			continue;
		velocity[0] = (sate_pos1.X - sate_pos.X) / 1e-3;
		velocity[1] = (sate_pos1.Y - sate_pos.Y) / 1e-3;
		velocity[2] = (sate_pos1.Z - sate_pos.Z) / 1e-3;
		velocity[3] = (clk1 - clk) * velocity_c / 1e-3;
		double lamda = (1e-6 * velocity_c / f);
		double len = sqrt(SQR(RcvPos.X - sate_pos.X) + SQR(RcvPos.Y - sate_pos.Y) + SQR(RcvPos.Z - sate_pos.Z));
		double w_pos = measure - (len + dt_Rcv - velocity_c * clk + Hopefield(sate_pos, RcvPos, Sates[i]->SYS));
		double v0 = ((sate_pos.X - RcvPos.X) * velocity[0] + (sate_pos.Y - RcvPos.Y) * velocity[1] + (sate_pos.Z - RcvPos.Z) * velocity[2]) / len;
		double w_vel = -lamda * Sates[i]->DOPPLER[Index] - (v0 - velocity[3]);
		if (abs(w_pos) > 10)
			continue;
		if (abs(w_vel) > 10)
			continue;
		double l = (RcvPos.X - sate_pos.X) / len;
		double m = (RcvPos.Y - sate_pos.Y) / len;
		double n = (RcvPos.Z - sate_pos.Z) / len;
		B_pos_new(0, 0) = l;
		B_pos_new(0, 1) = m;
		B_pos_new(0, 2) = n;
		B_pos_new(0, 3) = 1;
		l_pos_new(0, 0) = w_pos;
		l_vel_new(0, 0) = w_vel;
		B.conservativeResize(B.rows() + 1, B.cols());
		B.bottomRows(1) = B_pos_new;
		l_Pos.conservativeResize(l_Pos.rows() + 1, l_Pos.cols());
		l_Pos.bottomRows(1) = l_pos_new;
		l_Vel.conservativeResize(l_Vel.rows() + 1, l_Vel.cols());
		l_Vel.bottomRows(1) = l_vel_new;
		switch (sys)
		{
		case SYS_GPS:
			*(data->KF_SATES) += "G" + to_string(Sates[i]->PRN);
			break;
		case SYS_BDS:
			*(data->KF_SATES) += "C" + to_string(Sates[i]->PRN);
			break;
		default:
			break;
		}
		ROWS++;
	}
	if (ROWS == 0)
		return -0;
	data->KF->set_H(getH(B));
	data->KF->set_R(getR(ROWS));
	data->KF->set_Z(getz(l_Pos, l_Vel));
	switch (sys)
	{
	case SYS_GPS:
		data->KF_GPS_num = ROWS;
		break;
	case SYS_BDS:
		data->KF_BDS_num = ROWS;
		break;
	default:
		break;
	}
	return ROWS;
}

double get_measure(Satellate* Sate, Configure cfg, EPHEMERIS* eph)
{
	double measure = 0;
	double f1 = 0;
	double f2 = 0;
	switch (Sate->SYS)
	{
	case SYS_GPS:
		f1 = cfg.GPS_Cfg.f1;
		f2 = cfg.GPS_Cfg.f2;
		break;
	case SYS_BDS:
		f1 = cfg.BDS_Cfg.f1;
		f2 = cfg.BDS_Cfg.f2;
		break;
	default:
		break;
	}

	if (Sate->Outlier)
		return 0;
	int prn = Sate->PRN;
	int Index1 = 0;
	int Index2 = 0;
	while (CODE2FREQ(decode_SYN(Sate->SYS, Sate->SYG_TYPE[Index1])) != f1 && Index1 < MAXNUM)
		Index1++;
	while (CODE2FREQ(decode_SYN(Sate->SYS, Sate->SYG_TYPE[Index2])) != f2 && Index2 < MAXNUM)
		Index2++;
	if (Index1 == MAXNUM)
		return 0;
	if (!(Sate->LOCK_PSE[Index1] && Sate->LOCK_PHA[Index1]))
		return 0;

	if (cfg.SYS_num == 1)
		measure = Sate->PSERA[Index1];
	if (cfg.SYS_num = 2)
	{
		if (Sate->Phase_NUM < 2)
			return 0;
		if (Index2 == MAXNUM)
			return 0;
		if (!(Sate->LOCK_PSE[Index2] && Sate->LOCK_PHA[Index2]))
			return 0;
		measure = SQR(f1) * Sate->PSERA[Index1] / (SQR(f1) - SQR(f2)) - SQR(f2) * Sate->PSERA[Index2] / (SQR(f1) - SQR(f2));
		if (Sate->SYS == SYS_BDS)
			measure += velocity_c * SQR(f1 / f2) * eph->T_GD1 / (1 - SQR(f1 / f2));
	}

	return measure;
}

unsigned int LS_SPV(DATA_SET* data, Configure cfg)
{
	int val = 0;
	for (int i = 0; i < 4; i++)
	{
		*data->LS_SATES = "";
		data->LS_Pos->reset();
		data->LS_Vel->reset();
		val = data->Set_LS(cfg);
		if (val)
		{
			data->LS_Pos->ELS();
			data->LS_Vel->LS();
			data->LS_result = Success_Solve;
		}
	}

	return val;
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
	else if (cfg.SYS_num == 2)
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

unsigned int KF_SPV(DATA_SET* data, double dt_e, Configure cfg)
{
	int val = 0;
	if (data->KF_first && data->LS_result == Success_Solve)
	{
		MatrixXd state(8, 1);
		state.block(0, 0, 3, 1) = data->LS_Pos->X.block(0, 0, 3, 1);
		if (cfg.SYS_num == 1)
		{
			state.block(4, 0, 3, 1) = data->LS_Vel->X.block(0, 0, 3, 1);
			state(3, 0) = data->LS_Pos->X(3, 0);
			state(7, 0) = data->LS_Vel->X(3, 0);
		}
		else if (cfg.SYS_num == 2)
		{
			state.conservativeResize(9, 1);
			state.block(5, 0, 3, 1) = data->LS_Vel->X.block(0, 0, 3, 1);
			state(3, 0) = data->LS_Pos->X(3, 0);
			state(4, 0) = data->LS_Pos->X(4, 0);
			state(8, 0) = data->LS_Vel->X(3, 0);
		}
		data->KF->setState(state);
	}

	*data->KF_SATES = "";
	data->KF->reset();
	data->KF->set_A(getA(dt_e, cfg));
	data->KF->set_Q(getQ(dt_e, cfg));
	data->KF->predict();
	val = data->Set_KF(cfg);
	if (val)
	{
		data->KF->update();
		data->KF_result = Success_Solve;
	}
	return val;
}

int DATA_SET::Set_KF(Configure cfg)
{
	int val = 0;
	temp_ref.topRows(3) = KF->getState_minus().block(0, 0, 3, 1);
	double dt_G = 0;
	double dt_C = 0;
	if (cfg.SYS_num == 1)
	{
		if (cfg.GPS_Cfg.used)
			dt_G = KF->getState_minus()(3, 0);
		if (cfg.BDS_Cfg.used)
			dt_C = KF->getState_minus()(3, 0);
	}
	else if (cfg.SYS_num == 2)
	{
		dt_G = KF->getState_minus()(3, 0);
		dt_C = KF->getState_minus()(4, 0);
	}
	if (cfg.GPS_Cfg.used)
	{
		temp_ref(3, 0) = dt_G;
		val = setup_KF(this, cfg, SYS_GPS);
	}
	if (cfg.BDS_Cfg.used)
	{
		temp_ref(3, 0) = dt_C;
		val = setup_KF(this, cfg, SYS_BDS);
	}
	return val;
}