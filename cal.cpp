#include "cal.h"
#include "data.h"

/* 函数名: SQR
	输入参数:
	 - x: 要计算平方的值
	返回值: 输入值的平方
*/
double SQR(double x)
{
	return x * x;
}

// 函数名: Len
// 输入参数:
//   - pos: 包含三维坐标(X, Y, Z)的结构体
// 返回值: 坐标点的长度（模）
double Len(XYZ* pos)
{
	return sqrt(SQR(pos->X) + SQR(pos->Y) + SQR(pos->Z));
}

// 函数名: degree2rad
// 输入参数:
//   - degree: 角度值（度）
// 返回值: 角度值对应的弧度值
double degree2rad(double degree)
{
	// 将角度限制在[0, 360)范围内
	while (degree > 360 || degree < 0)
	{
		if (degree > 360)
			degree -= 360;
		if (degree < 0)
			degree += 360;
	}
	return degree * Pi / 180; // 转换为弧度
}

// 函数名: rad2degree
// 输入参数:
//   - rad: 弧度值
// 返回值: 弧度值对应的角度值（度）
double rad2degree(double rad)
{
	return rad * 180 / Pi; // 转换为角度
}

// 函数名: Cal_PDOP
// 输入参数:
//   - Qxx: 协方差矩阵
// 返回值: PDOP值
double Cal_PDOP(MatrixXd Qxx)
{
	if (Qxx.rows() < 3 || Qxx.cols() < 3)
		return 0;
	return sqrt(Qxx(0, 0) + Qxx(1, 1) + Qxx(2, 2));
}

// 函数名: Cal_LEAST_SQR
// 输入参数:
//   - B: 设计矩阵
//   - l: 观测残差
//   - P: 观测权重矩阵
//   - Qxx: 输出参数，协方差矩阵
//   - x: 输出参数，最小二乘解
//   - thegma: 输出参数，单位权中误差
//   - DOP: 输出参数，PDOP值
// 返回值: 操作结果，成功返回1，否则返回0
unsigned int Cal_LEAST_SQR(MatrixXd B, MatrixXd l, MatrixXd P, MatrixXd* Qxx, MatrixXd& x, double* thegma, double* DOP)
{
	// 检查输入矩阵的维度是否合法
	if (B.rows() != P.rows() || B.rows() != l.rows())
		return 0;
	if (B.rows() < B.cols())
		return 0;

	// 计算协方差矩阵
	*Qxx = (B.transpose() * P * B).inverse();

	// 计算最小二乘解
	x = *Qxx * B.transpose() * P * l;

	// 计算观测残差
	MatrixXd v = B * x - l;

	// 计算单位权中误差
	int m = B.rows() - B.cols();
	*thegma = sqrt(((v.transpose() * P * v) / m)(0, 0));

	// 计算PDOP值
	*DOP = Cal_PDOP(*Qxx);

	return 1; // 成功
}

// 函数名: decode_SYN
// 输入参数:
//   - sys: 系统编号
//   - signal: 信号编号
// 返回值: 解码后的信号编号，未知信号返回UNKOWN
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

// 函数名: CODE2FREQ
// 输入参数:
//   - code: 信号编号
// 返回值: 对应信号编号的频率值
double CODE2FREQ(int code)
{
	switch (code)
	{
	case 0:
		return 0;
	case 1:
		return L1;
	case 2:
	case 3:
		return L2;
	case 4:
		return L5;
	case 5:
	case 6:
		return L1;
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

// 函数名: CORRECT_CLK
// 输入参数:
//   - t: 当前时刻
//   - eph: 指向EPHEMERIS结构的指针
// 返回值: 改正后的钟差
double CORRECT_CLK(double t, EPHEMERIS* eph)
{
	double correct_clk = t;
	for (int i = 0; i < 10; i++)
	{
		correct_clk = t - (eph->a_f0 + eph->a_f1 * (correct_clk - eph->toc) + eph->a_f2 * (correct_clk - eph->toc) * (correct_clk - eph->toc));
	}
	return eph->a_f0 + eph->a_f1 * (correct_clk - eph->toc) + eph->a_f2 * (correct_clk - eph->toc) * (correct_clk - eph->toc);
}

// 函数名: TGD
// 输入参数:
//   - e: 指向EPHEMERIS结构的指针
//   - f: 频率值
//   - sys: 卫星系统
// 返回值: TGD值
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

// 函数名: SAT_POS_CAL
// 输入参数:
//   - t: 观测时刻
//   - eph: 卫星星历数据结构指针
//   - Sat_Pos: 卫星位置结果结构指针
//   - clk: 卫星钟差
//   - dt: 传播时间修正值
//   - SYS: 卫星系统类型 (SYS_GPS 或 SYS_BDS)
// 返回值: 1 表示计算成功

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
	dt += clk;
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
		if (fabs(eph->i0 - 0.0873) < 0.1 && fabs(eph->sqrt_A - 6493) < 1) // 通过轨道倾角和轨道根数判断是否为GEO卫星  i: 5/deg sqrt_A: 6493/sqrt_meter
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

// 函数名: Ele_Angle
// 输入参数:
//   - SatPos: 卫星直角坐标系下的位置
//   - RcvPos: 接收机直角坐标系下的位置
//   - sys: 卫星系统类型 (SYS_GPS 或 SYS_BDS)
// 返回值: 卫星高度角

double Ele_Angle(XYZ SatPos, XYZ RcvPos, int sys)
{
	// 将卫星坐标转换到ENU坐标系
	XYZ Satenu = XYZ2ENU(RcvPos, SatPos, sys);

	// 计算卫星高度角
	return asin(Satenu.Z / Len(&Satenu));
}

// 函数名: Hopefield
// 输入参数:
//   - E: 卫星高度角
//   - H: 接收机海拔高度
// 返回值: Hopefield对流层改正值

double Hopefield(double E, double H)
{
	// 计算标准大气压力和温度
	double Ts = T0 - 0.0065 * (H - H0);
	double hd = 40136 + 148.72 * (Ts - 273.16);

	// 计算标准大气压和相对湿度
	double Ps = P0 * pow(1 - 0.0000226 * (H - H0), 5.225);
	double RH = RH0 * exp(-0.0006396 * (H - H0));

	// 计算饱和水汽压力
	double es = RH * exp(-37.2465 + 0.213166 * Ts - 0.000256908 * Ts * Ts);

	// 计算方位角对应的方向余弦
	double md = sin(degree2rad(sqrt(E * E + 6.25)));
	double mw = sin(degree2rad(sqrt(E * E + 2.25)));

	// 计算水汽分压和水汽刻度高
	double hw = 11000;
	double ZHD = 155.2 * 1e-7 * Ps * (hd - H) / Ts;
	double ZWD = 155.2 * 1e-7 * 4810 * es * (hw - H) / (Ts * Ts);

	// 返回Hopefield对流层改正值
	return ZHD / md + ZWD / mw;
}

// 函数名: Hopefield
// 输入参数:
//   - SatPos: 卫星直角坐标系下的位置
//   - RcvPos: 接收机直角坐标系下的位置
//   - sys: 卫星系统类型 (SYS_GPS 或 SYS_BDS)
// 返回值: Hopefield对流层改正值

double Hopefield(XYZ SatPos, XYZ RcvPos, int sys)
{
	// 若卫星坐标为零向量，返回0
	if (SatPos.X == 0 && SatPos.Y == 0 && SatPos.Z == 0)
	{
		return 0;
	}

	// 计算卫星高度角
	double E = rad2degree(Ele_Angle(SatPos, RcvPos, sys));

	// 获取接收机海拔高度
	double H = 0;
	switch (sys)
	{
	case SYS_GPS:
		H = XYZ2BLH(RcvPos, WGS84_e2, WGS84_a).Height;
		break;
	case SYS_BDS:
		H = XYZ2BLH(RcvPos, CGCS2000_e2, CGCS2000_a).Height;
		break;
	default:
		break;
	}

	// 判断海拔高度是否在有效范围内，计算Hopefield对流层改正值
	if (H < 20e3 && H > -100)
	{
		return Hopefield(E, H);
	}
	else
	{
		return 0;
	}
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

// 函数名: setup_LS
// 输入参数:
//   - data: 数据集结构体指针
//   - cfg: 配置信息结构体
//   - sys: 卫星系统类型 (SYS_GPS 或 SYS_BDS)
// 返回值: LS定位矩阵行数

unsigned int setup_LS(DATA_SET* data, Configure cfg, int sys)
{
	// 声明变量
	int ROWS = 0;
	vector<Satellate*> Sates;
	EPHEMERIS** eph;
	double f = 0;

	// 根据卫星系统类型选择相应的卫星数据和星历数据
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
	if (Sates.size() < cfg.SYS_num + 3)
	{
		data->LS_result = OBS_DATA_Loss;
		return 0;
	}
	// 获取接收机位置和接收机钟差
	XYZ RcvPos = get_XYZ(data->temp_ref.block(0, 0, 3, 1));
	double dt_Rcv = data->temp_ref(3, 0);

	// 初始化变量
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

	double p_pos = 0;

	// 循环遍历卫星
	for (int i = 0; i < Sates.size(); i++)
	{
		int prn = Sates[i]->PRN;

		// 检查星历数据是否匹配
		if (eph[prn - 1]->PRN != prn)
			continue;

		// 查找卫星对应的频率索引
		int Index = 0;
		while (CODE2FREQ(decode_SYN(Sates[i]->SYS, Sates[i]->SYG_TYPE[Index])) != f && Index < MAXNUM)
			Index++;
		if (Index == MAXNUM)
			continue;

		// 获取测量值
		double measure = get_measure(Sates[i], cfg, eph[prn - 1], p_pos);
		if (!measure)
			continue;

		// 计算卫星位置、钟差
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

		// 计算卫星位置和速度
		for (int j = 0; j < 3; j++)
		{
			clk = CORRECT_CLK(ts - clk, eph[prn - 1]);
			SAT_POS_CAL(ts - clk + tgd, eph[prn - 1], &sate_pos, clk, Sates[i]->PSERA[0] / velocity_c, Sates[i]->SYS);
			clk1 = CORRECT_CLK(ts - clk1 + 1e-3, eph[prn - 1]);
			SAT_POS_CAL(ts - clk1 + tgd + 1e-3, eph[prn - 1], &sate_pos1, clk1, Sates[i]->PSERA[0] / velocity_c, Sates[i]->SYS);
		}

		// 检查是否为首次定位
		if (!data->LS_first)
		{
			if (Ele_Angle(sate_pos, RcvPos, Sates[i]->SYS) < degree2rad(cfg.Ele_Mask))
				continue;
		}

		// 计算卫星速度
		velocity[0] = (sate_pos1.X - sate_pos.X) / 1e-3;
		velocity[1] = (sate_pos1.Y - sate_pos.Y) / 1e-3;
		velocity[2] = (sate_pos1.Z - sate_pos.Z) / 1e-3;
		velocity[3] = (clk1 - clk) * velocity_c / 1e-3;

		// 计算方位角对应的波长和长度
		double lamda = (1e-6 * velocity_c / f);
		double len = sqrt(SQR(RcvPos.X - sate_pos.X) + SQR(RcvPos.Y - sate_pos.Y) + SQR(RcvPos.Z - sate_pos.Z));

		// 计算预测值
		double w_pos = measure - (len + dt_Rcv - velocity_c * clk + Hopefield(sate_pos, RcvPos, Sates[i]->SYS));
		double v0 = ((sate_pos.X - RcvPos.X) * velocity[0] + (sate_pos.Y - RcvPos.Y) * velocity[1] + (sate_pos.Z - RcvPos.Z) * velocity[2]) / len;
		double w_vel = -lamda * Sates[i]->DOPPLER[Index] - (v0 - velocity[3]);

		// 检查是否为首次定位
		if (!data->LS_first)
		{
			if (abs(w_pos) > cfg.w_thresh)
				continue;
		}

		// 构造最小二乘问题的设计矩阵B和观测矩阵l
		double l = (RcvPos.X - sate_pos.X) / len;
		double m = (RcvPos.Y - sate_pos.Y) / len;
		double n = (RcvPos.Z - sate_pos.Z) / len;
		B_pos_new(0, 0) = l;
		B_pos_new(0, 1) = m;
		B_pos_new(0, 2) = n;
		B_pos_new(0, 3) = 1;
		l_pos_new(0, 0) = w_pos;
		l_vel_new(0, 0) = w_vel;
		ROWS++;

		// 更新B、l、P矩阵
		B.conservativeResize(B.rows() + 1, B.cols());
		B.bottomRows(1) = B_pos_new;
		l_Pos.conservativeResize(l_Pos.rows() + 1, l_Pos.cols());
		l_Pos.bottomRows(1) = l_pos_new;
		l_Vel.conservativeResize(l_Vel.rows() + 1, l_Vel.cols());
		l_Vel.bottomRows(1) = l_vel_new;
		if (ROWS == 1)
		{
			P_Pos = MatrixXd::Identity(1, 1);
			P_Pos(0, 0) = p_pos;
		}
		else
		{
			MatrixXd temp_P = MatrixXd::Identity(ROWS, ROWS);
			temp_P.block(0, 0, P_Pos.rows(), P_Pos.cols()) = P_Pos;
			temp_P(ROWS - 1, ROWS - 1) = p_pos;
			P_Pos = temp_P;
		}

		// 根据卫星系统类型记录卫星PRN
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
	}

	// 初始化权阵矩阵
	P_Pos = MatrixXd::Identity(ROWS, ROWS);
	P_Vel = MatrixXd::Identity(ROWS, ROWS);

	// 若有卫星数据，则更新LS_Pos和LS_Vel矩阵
	if (ROWS >= cfg.SYS_num + 3)
	{
		data->LS_Pos->set_B_Pos(B);
		data->LS_Pos->set_l(l_Pos);
		data->LS_Pos->set_P(P_Pos);
		data->LS_Vel->set_B_Vel(B);
		data->LS_Vel->set_l(l_Vel);
		data->LS_Vel->set_P(P_Vel);
	}
	else
	{
		data->LS_result = OBS_DATA_Loss;
		return 0;
	}

	// 根据卫星系统类型记录卫星数量
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

	// 返回LS定位矩阵行数
	return ROWS;
}


// 函数名: setup_KF
// 输入参数:
//   - data: 数据集结构体指针
//   - cfg: 配置信息结构体
//   - sys: 卫星系统类型 (SYS_GPS 或 SYS_BDS)
// 返回值: KF定位矩阵行数

unsigned int setup_KF(DATA_SET* data, Configure cfg, int sys)
{
	// 声明变量
	int ROWS = 0;
	vector<Satellate*> Sates;
	EPHEMERIS** eph;
	double f = 0;

	// 根据卫星系统类型选择相应的卫星数据和星历数据
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

	// 获取接收机位置和接收机钟差
	XYZ RcvPos = get_XYZ(data->temp_ref.block(0, 0, 3, 1));
	MatrixXd RcvVel(4, 4);
	RcvVel.block(0, 0, 3, 1) = data->KF->getState_minus().block(3, 0, 3, 1);
	RcvVel(3,0)= data->KF->getState_minus()(cfg.SYS_num + 6, 0);
	double dt_Rcv = data->temp_ref(3, 0);
	XYZ sate_pos, sate_pos1;
	double clk = 0;
	double clk1 = 0;
	double velocity[4] = { 0, 0, 0, 0 };
	MatrixXd B, l_Pos, l_Vel;
	B = MatrixXd::Zero(0, 4);
	l_Vel = l_Pos = MatrixXd::Zero(0, 1);
	MatrixXd B_pos_new = MatrixXd::Zero(1, 4);
	MatrixXd l_pos_new = MatrixXd::Zero(1, 1);
	MatrixXd l_vel_new = MatrixXd::Zero(1, 1);
	MatrixXd P_Pos;
	double p_pos = 0;
	// 循环遍历卫星
	for (int i = 0; i < Sates.size(); i++)
	{
		int prn = Sates[i]->PRN;

		// 检查星历数据是否匹配
		if (eph[prn - 1]->PRN != prn)
			continue;

		// 查找卫星对应的频率索引
		int Index = 0;
		while (CODE2FREQ(decode_SYN(Sates[i]->SYS, Sates[i]->SYG_TYPE[Index])) != f && Index < MAXNUM)
			Index++;
		if (Index == MAXNUM)
			continue;

		// 获取测量值
		double measure = get_measure(Sates[i], cfg, eph[prn - 1], p_pos);
		if (!measure)
			continue;

		// 计算卫星位置、钟差
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

		// 计算卫星位置和速度
		for (int j = 0; j < 3; j++)
		{
			clk = CORRECT_CLK(ts - clk, eph[prn - 1]);
			SAT_POS_CAL(ts - clk + tgd, eph[prn - 1], &sate_pos, clk, Sates[i]->PSERA[0] / velocity_c, Sates[i]->SYS);
			clk1 = CORRECT_CLK(ts - clk1 + 1e-3, eph[prn - 1]);
			SAT_POS_CAL(ts - clk1 + tgd + 1e-3, eph[prn - 1], &sate_pos1, clk1, Sates[i]->PSERA[0] / velocity_c, Sates[i]->SYS);
		}

		// 检查仰角是否小于阈值
		if (Ele_Angle(sate_pos, RcvPos, Sates[i]->SYS) < degree2rad(cfg.Ele_Mask))
			continue;

		// 计算卫星速度
		velocity[0] = (sate_pos1.X - sate_pos.X) / 1e-3;
		velocity[1] = (sate_pos1.Y - sate_pos.Y) / 1e-3;
		velocity[2] = (sate_pos1.Z - sate_pos.Z) / 1e-3;
		velocity[3] = (clk1 - clk) * velocity_c / 1e-3;

		// 计算方位角对应的波长和长度
		double lamda = (1e-6 * velocity_c / f);
		double len = sqrt(SQR(RcvPos.X - sate_pos.X) + SQR(RcvPos.Y - sate_pos.Y) + SQR(RcvPos.Z - sate_pos.Z));
		// 构造卡尔曼滤波问题的设计矩阵B和观测矩阵l
		double l = (RcvPos.X - sate_pos.X) / len;
		double m = (RcvPos.Y - sate_pos.Y) / len;
		double n = (RcvPos.Z - sate_pos.Z) / len;
		B_pos_new(0, 0) = l;
		B_pos_new(0, 1) = m;
		B_pos_new(0, 2) = n;
		B_pos_new(0, 3) = 1;
		// 计算预测值
		double w_pos = measure - (len + dt_Rcv - velocity_c * clk + Hopefield(sate_pos, RcvPos, Sates[i]->SYS));
		double v0 = ((sate_pos.X - RcvPos.X) * velocity[0] + (sate_pos.Y - RcvPos.Y) * velocity[1] + (sate_pos.Z - RcvPos.Z) * velocity[2]) / len;
		double w_vel = -lamda * Sates[i]->DOPPLER[Index] - (v0 - velocity[3]) - (B_pos_new * RcvVel)(0, 0);

		// 检查预测值是否大于阈值
		if (abs(w_pos) > cfg.w_thresh)
			continue;

		// 检查预测速度是否大于阈值
		if (abs(w_vel) > cfg.w_thresh)
			continue;

		l_pos_new(0, 0) = w_pos;
		l_vel_new(0, 0) = w_vel;
		ROWS++;

		// 更新B、l矩阵
		B.conservativeResize(B.rows() + 1, B.cols());
		B.bottomRows(1) = B_pos_new;
		l_Pos.conservativeResize(l_Pos.rows() + 1, l_Pos.cols());
		l_Pos.bottomRows(1) = l_pos_new;
		l_Vel.conservativeResize(l_Vel.rows() + 1, l_Vel.cols());
		l_Vel.bottomRows(1) = l_vel_new;
		if (ROWS == 1)
		{
			P_Pos = MatrixXd::Identity(1, 1);
			P_Pos(0, 0) = p_pos;
		}
		else
		{
			MatrixXd temp_P = MatrixXd::Identity(ROWS, ROWS);
			temp_P.block(0, 0, P_Pos.rows(), P_Pos.cols()) = P_Pos;
			temp_P(ROWS - 1, ROWS - 1) = p_pos;
			P_Pos = temp_P;
		}

		// 根据卫星系统类型记录卫星PRN
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
	}

	// 如果卫星数为0，返回0
	if (ROWS == 0)
		return 0;

	MatrixXd R = MatrixXd::Identity(ROWS + ROWS, ROWS + ROWS);
	R.block(0, 0, ROWS, ROWS) = P_Pos;
	R.block(ROWS, ROWS, ROWS, ROWS) = MatrixXd::Identity(ROWS, ROWS) * 0.18;

	// 更新Kalman Filter对象的H、R、Z矩阵
	data->KF->set_H(getH(B));
	data->KF->set_R(getR(ROWS));
	data->KF->set_Z(getz(l_Pos, l_Vel));

	// 根据卫星系统类型记录卫星数量
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

	// 返回KF定位矩阵行数
	return ROWS;
}


	// 函数名: get_measure
	// 输入参数:
	//   - Sate: 卫星信息结构体指针
	//   - cfg: 配置信息结构体
	//   - eph: 星历数据结构体指针
	//   - p: 权重
	// 返回值: 测量值

	double get_measure(Satellate* Sate, Configure cfg, EPHEMERIS* eph, double &p)
	{
		// 声明变量
		double measure = 0;
		double f1 = 0;
		double f2 = 0;

		// 根据卫星系统类型选择相应的频率
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

		// 如果卫星标记为异常值，返回测量值为0
		if (Sate->Outlier)
			return 0;

		int prn = Sate->PRN;
		int Index1 = 0;
		int Index2 = 0;

		// 查找卫星对应的频率索引
		while (CODE2FREQ(decode_SYN(Sate->SYS, Sate->SYG_TYPE[Index1])) != f1 && Index1 < MAXNUM)
			Index1++;
		while (CODE2FREQ(decode_SYN(Sate->SYS, Sate->SYG_TYPE[Index2])) != f2 && Index2 < MAXNUM)
			Index2++;

		// 如果未找到对应频率索引，返回测量值为0
		if (Index1 == MAXNUM)
			return 0;

		// 如果指定的频率1没有相应的相位和伪距观测值，返回测量值为0
		if (!(Sate->LOCK_PSE[Index1] && Sate->LOCK_PHA[Index1]))
			return 0;

		// 根据配置的系统数目计算测量值
		if (cfg.phase_num == 1)
		{
			measure = Sate->PSERA[Index1];
			p = SQR(Sate->PSE_PREC[Index1]);
		}
		else if (cfg.phase_num == 2)
		{
			// 如果相位观测值数量小于2，返回测量值为0
			if (Sate->Phase_NUM < 2)
				return 0;

			// 如果未找到对应频率索引，返回测量值为0
			if (Index2 == MAXNUM)
				return 0;

			// 如果指定的频率2没有相应的相位和伪距观测值，返回测量值为0
			if (!(Sate->LOCK_PSE[Index2] && Sate->LOCK_PHA[Index2]))
				return 0;

			// 计算测量值（无电离层组合）
			measure = SQR(f1) * Sate->PSERA[Index1] / (SQR(f1) - SQR(f2)) - SQR(f2) * Sate->PSERA[Index2] / (SQR(f1) - SQR(f2));
			p = SQR(SQR(f1) * Sate->PSE_PREC[Index1] / (SQR(f1) - SQR(f2))) + SQR(SQR(f2) * Sate->PSE_PREC[Index2] / (SQR(f1) - SQR(f2)));

			// 对于BDS系统，加上TGD1对应的钟差改正
			if (Sate->SYS == SYS_BDS)
				measure += velocity_c * SQR(f2 / f1) * eph->T_GD1 / (1 - SQR(f2 / f1));
		}

		// 返回计算得到的测量值
		return measure;
	}


// 函数名: LS_SPV
// 输入参数:
//   - data: 数据集，包含卡尔曼滤波器、最小二乘解等数据
//   - cfg: 配置信息，包括系统编号和传感器配置
// 返回值: 操作结果，成功返回0，否则返回非零值
unsigned int LS_SPV(DATA_SET* data, Configure cfg)
{
	int val = 0;
	

	// 迭代4次
	for (int i = 0; i < 4; i++)
	{
		// 清空卫星信息
		*data->LS_SATES = "";
		// 重置最小二乘解和速度
		data->LS_Pos->reset();
		data->LS_Vel->reset();
		// 设置最小二乘
		val = data->Set_LS(cfg);

		// 如果成功设置最小二乘，则执行最小二乘解算
		if (val)
		{
			data->LS_Pos->ELS();
			data->LS_Vel->LS();
			data->LS_result = Success_Solve;
		}
	}

	return val;
}

// 函数名: Set_LS
// 输入参数:
//   - cfg: 配置信息，包括系统数和传感器配置
// 返回值: 最小二乘问题的行数
int DATA_SET::Set_LS(Configure cfg)
{
	// 将最小二乘解的前三行赋值给临时参考数据
	temp_ref.topRows(3) = LS_Pos->X.block(0, 0, 3, 1);

	// 初始化GPS和BDS的接收机钟差
	double dt_G = 0;
	double dt_C = 0;

	// 根据系统数量设置接收机钟差
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

	// 设置GPS的接收机钟差
	if (cfg.GPS_Cfg.used)
	{
		temp_ref(3, 0) = dt_G;
		setup_LS(this, cfg, SYS_GPS);
	}

	// 设置BDS的接收机钟差
	if (cfg.BDS_Cfg.used)
	{
		temp_ref(3, 0) = dt_C;
		setup_LS(this, cfg, SYS_BDS);
	}

	// 返回最小二乘的行数
	return LS_Pos->B.rows();
}


// 函数名: KF_SPV
// 输入参数:
//   - data: 数据集，包含卡尔曼滤波器、最小二乘解等数据
//   - dt_e: 时间步长
//   - cfg: 配置信息，包括系统编号和传感器配置
// 返回值: 操作结果，成功返回0，否则返回非零值
unsigned int KF_SPV(DATA_SET* data, double dt_e, Configure cfg) {
	int val = 0;

	// 如果是第一次滤波且最小二乘解成功
	if (data->KF_first && data->LS_result == Success_Solve) {
		// 从最小二乘解中提取状态向量
		MatrixXd state(8, 1);
		state.block(0, 0, 3, 1) = data->LS_Pos->X.block(0, 0, 3, 1);
		state.block(3, 0, 3, 1) = data->LS_Vel->X.block(0, 0, 3, 1);
		// 根据系统编号配置状态向量
		if (cfg.SYS_num == 1) {
			state(6, 0) = data->LS_Pos->X(3, 0);
			state(7, 0) = data->LS_Vel->X(3, 0);
		}
		else if (cfg.SYS_num == 2) {
			state.conservativeResize(9, 1);
			state(6, 0) = data->LS_Pos->X(3, 0);
			state(7, 0) = data->LS_Pos->X(4, 0);
			state(8, 0) = data->LS_Vel->X(3, 0);
		}

		// 设置卡尔曼滤波器状态
		data->KF->setState(state);
		data->KF_GPS_num = data->LS_GPS_num;
		data->KF_BDS_num = data->LS_BDS_num;
		*data->KF_SATES = *(data->LS_SATES);
		data->KF_result = Success_Solve;
		return 1;
	}

	// 重置卡尔曼滤波器并进行预测
	*data->KF_SATES = "";
	data->KF->reset();
	data->KF->set_A(getA(dt_e, cfg));
	data->KF->set_Q(getQ(dt_e, cfg));
	data->KF->predict();

	// 设置卡尔曼滤波器的状态
	val = data->Set_KF(cfg);

	// 如果设置成功，执行更新步骤
	if (val) {
		data->KF->update();
		data->KF_result = Success_Solve;
	}

	return val;
}

// 函数名: DATA_SET::Set_KF
// 输入参数:
//   - cfg: 配置信息，包括系统编号和传感器配置
// 返回值: 操作结果，成功返回0，否则返回非零值
int DATA_SET::Set_KF(Configure cfg) {
	int val = 0;
	temp_ref.topRows(3) = KF->getState_minus().block(0, 0, 3, 1);
	double dt_G = 0;
	double dt_C = 0;

	// 根据系统编号配置GPS和BDS的接收机钟差
	if (cfg.SYS_num == 1) {
		if (cfg.GPS_Cfg.used)
			dt_G = KF->getState_minus()(6, 0);
		if (cfg.BDS_Cfg.used)
			dt_C = KF->getState_minus()(6, 0);
	}
	else if (cfg.SYS_num == 2) {
		dt_G = KF->getState_minus()(6, 0);
		dt_C = KF->getState_minus()(7, 0);
	}

	// 根据配置信息设置GPS的状态向量
	if (cfg.GPS_Cfg.used) {
		temp_ref(3, 0) = dt_G;
		val = setup_KF(this, cfg, SYS_GPS);
	}

	// 根据配置信息设置BDS的状态向量
	if (cfg.BDS_Cfg.used) {
		temp_ref(3, 0) = dt_C;
		val = setup_KF(this, cfg, SYS_BDS);
	}

	return val;
}