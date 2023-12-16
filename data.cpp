#include "data.h"
#include "cal.h"
#include <io.h>
#include <direct.h>

DATA_SET::DATA_SET(Configure cfg)
{
	OBSTIME = new GPSTIME();
	Pos = new MatrixXd();
	Vel = new MatrixXd();
	Q_Pos = new MatrixXd();
	Q_Vel = new MatrixXd();
	thegma_Pos = new double;
	thegma_Vel = new double;
	PDOP = new double;
	VDOP = new double;
	LS_SATES = new string();
	KF_SATES = new string();
	*Pos = MatrixXd::Zero(5, 1);
	*Vel = MatrixXd::Zero(4, 1);
	*Q_Pos = MatrixXd::Zero(5, 5);
	*Q_Vel = MatrixXd::Zero(4, 4);
	*thegma_Pos = 0;
	*thegma_Vel = 0;
	*PDOP = 0;
	*VDOP = 0;
	LS_GPS_num = 0;
	LS_BDS_num = 0;
	KF_GPS_num = 0;
	KF_BDS_num = 0;
	*LS_SATES = "";
	*KF_SATES = "";

	range = new OBS_DATA();
	range->OBS_TIME = new GPSTIME();
	for (int i = 0; i < GPS_SAT_QUAN; i++)
	{
		GPS_GF[i] = 0;
		GPS_MW[i] = 0;
		GPS_eph[i] = new EPHEMERIS();
	}
	for (int i = 0; i < BDS_SAT_QUAN; i++)
	{
		BDS_GF[i] = 0;
		BDS_MW[i] = 0;
		BDS_eph[i] = new EPHEMERIS();
	}

	LS_result = UN_Solve;
	KF_result = UN_Solve;

	LS_first = true;
	KF_first = true;
	temp_ref = MatrixXd::Zero(4, 1);

	LS_Pos = new Least_Squares(cfg);
	LS_Vel = new Least_Squares();

	int row = 7 + cfg.SYS_num;
	MatrixXd P_ = MatrixXd::Zero(row, row);
	P_.block(0, 0, row - 4, row - 4) = MatrixXd::Identity(row - 4, row - 4) * 10 * 10;
	P_.block(row - 4, row - 4, 4, 4) = MatrixXd::Identity(4, 4) * 0.1 * 0.1;
	KF = new KalmanFilter(MatrixXd::Zero(7 + cfg.SYS_num, 1), P_);
	Real_Pos = new XYZ(-2267807.853, 5009320.431, 3221020.875);
}

void DATA_SET::reset()
{
	for (vector<Satellate*>::iterator it = range->GPS_SATE.begin(); it != range->GPS_SATE.end(); it++)
	{
		if (NULL != *it)
		{
			delete* it;
			*it = NULL;
		}
	}
	range->GPS_SATE.clear();
	for (vector<Satellate*>::iterator it = range->BDS_SATE.begin(); it != range->BDS_SATE.end(); it++)
	{
		if (NULL != *it)
		{
			delete* it;
			*it = NULL;
		}
	}
	range->BDS_SATE.clear();
}

int DATA_SET::LS_print(Configure cfg)
{
	//XYZ xyz = get_XYZ(( * Pos).block(0, 0, 3, 1));
	//BLH blh = XYZ2BLH(xyz, WGS84_e2, WGS84_a);
	//XYZ enu = XYZ2ENU(*Real_Pos, xyz, SYS_GPS);
	//MatrixXd R = get_Rot(degree2rad(blh.Lat), degree2rad(blh.Lon));
	//MatrixXd Q_local;
	//double m_H, m_V, B, L;
	//Q_local = R * ((*Q_Pos).block(0, 0, 3, 3)) * R.transpose();
	//m_H = (*thegma_Pos) * sqrt(Q_local(2, 2));
	//m_V = (*thegma_Pos) * sqrt(Q_local(0, 0) + Q_local(1, 1));

	XYZ xyz = get_XYZ(LS_Pos->X.block(0, 0, 3, 1));
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
	BLH blh = XYZ2BLH(xyz, WGS84_e2, WGS84_a);
	XYZ enu = XYZ2ENU(*Real_Pos, xyz, SYS_GPS);
	MatrixXd R = get_Rot(degree2rad(blh.Lat), degree2rad(blh.Lon));
	MatrixXd Q_local;
	double m_H, m_V, B, L;
	Q_local = R * ((*Q_Pos).block(0, 0, 3, 3)) * R.transpose();
	m_H = (LS_Pos->sigma) * sqrt(Q_local(2, 2));
	m_V = (LS_Pos->sigma) * sqrt(Q_local(0, 0) + Q_local(1, 1));
	switch (LS_result)
	{
	case UN_Solve:
		printf("GPSTIME: %d\t%.3f\tUN_Solve\n", OBSTIME->Week, OBSTIME->SecOfWeek);
		return 0;
	case Success_Solve:
		printf("GPSTIME: %d\t%.3f\tXYZ: %.4f\t%.4f\t%.4f\tBLH: %8.4f\t%8.4f\t%7.4f\tENU: %7.4f\t%7.4f\t%7.4f\t",
			OBSTIME->Week, OBSTIME->SecOfWeek,
			xyz.X, xyz.Y, xyz.Z,
			blh.Lat, blh.Lon, blh.Height,
			enu.X, enu.Y, enu.Z);
		if (cfg.GPS_Cfg.used)
			printf("GPS Clk : % 7.4f\t", dt_G);
		if (cfg.BDS_Cfg.used)
			printf("BDS Clk : % 7.4f\t", dt_C);
		printf("m_H : % 7.4f\tm_V : % 7.4f\tVelocity : % 8.4f\t % 8.4f\t % 8.4f\t % 8.4f\tthegma_P: % 6.4f\tthegma_V : % 6.4f\tPDOP : % 6.4f\tGPS : % 2d\tBDS : % 2d\t % s\n",
			m_H, m_V,
			LS_Vel->X(0, 0), LS_Vel->X(1, 0), LS_Vel->X(2, 0), LS_Vel->X(3, 0),
			LS_Pos->sigma, LS_Vel->sigma, Cal_PDOP(LS_Pos->Qxx),
			LS_GPS_num, LS_BDS_num,
			LS_SATES->c_str());
		return 1;
		break;
	case OBS_DATA_Loss:
		printf("GPSTIME: %d\t%.3f\tOBS_DATA_Loss\n", OBSTIME->Week, OBSTIME->SecOfWeek);
		return 0;
	case Epoch_Loss:
		printf("GPSTIME: %d\t%.3f\tEpoch_Loss\n", OBSTIME->Week, OBSTIME->SecOfWeek);
		return 0;
	case Set_UP_B_fail:
		printf("GPSTIME: %d\t%.3f\tSet_UP_B_fail\n", OBSTIME->Week, OBSTIME->SecOfWeek);
		return 0;
	default:
		break;
	}
	return 0;
}

int DATA_SET::LS_Filewrite(FILE* fpr, Configure cfg)
{
	XYZ xyz = get_XYZ(LS_Pos->X.block(0, 0, 3, 1));
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
	BLH blh = XYZ2BLH(xyz, WGS84_e2, WGS84_a);
	XYZ enu = XYZ2ENU(*Real_Pos, xyz, SYS_GPS);
	MatrixXd R = get_Rot(degree2rad(blh.Lat), degree2rad(blh.Lon));
	MatrixXd Q_local;
	double m_H, m_V, B, L;
	Q_local = R * ((*Q_Pos).block(0, 0, 3, 3)) * R.transpose();
	m_H = (LS_Pos->sigma) * sqrt(Q_local(2, 2));
	m_V = (LS_Pos->sigma) * sqrt(Q_local(0, 0) + Q_local(1, 1));
	switch (LS_result)
	{
	case UN_Solve:
		fprintf(fpr, "GPSTIME: %d\t%.3f\tUN_Solve\n", OBSTIME->Week, OBSTIME->SecOfWeek);
		return 0;
	case Success_Solve:
		fprintf(fpr, "GPSTIME: %d\t%.3f\tXYZ: %.4f\t%.4f\t%.4f\tBLH: %8.4f\t%8.4f\t%7.4f\tENU: %7.4f\t%7.4f\t%7.4f\t",
			OBSTIME->Week, OBSTIME->SecOfWeek,
			xyz.X, xyz.Y, xyz.Z,
			blh.Lat, blh.Lon, blh.Height,
			enu.X, enu.Y, enu.Z);
		if (cfg.GPS_Cfg.used)
			fprintf(fpr, "GPS Clk : % 7.4f\t", dt_G);
		if (cfg.BDS_Cfg.used)
			fprintf(fpr, "BDS Clk : % 7.4f\t", dt_C);
		fprintf(fpr, "m_H : % 7.4f\tm_V : % 7.4f\tVelocity : % 8.4f\t % 8.4f\t % 8.4f\t % 8.4f\tthegma_P: % 6.4f\tthegma_V : % 6.4f\tPDOP : % 6.4f\tGPS : % 2d\tBDS : % 2d\t % s\n",
			m_H, m_V,
			LS_Vel->X(0, 0), LS_Vel->X(1, 0), LS_Vel->X(2, 0), LS_Vel->X(3, 0),
			LS_Pos->sigma, LS_Vel->sigma, Cal_PDOP(LS_Pos->Qxx),
			LS_GPS_num, LS_BDS_num,
			LS_SATES->c_str());
		return 1;
		break;
	case OBS_DATA_Loss:
		fprintf(fpr, "GPSTIME: %d\t%.3f\tOBS_DATA_Loss\n", OBSTIME->Week, OBSTIME->SecOfWeek);
		return 0;
	case Epoch_Loss:
		fprintf(fpr, "GPSTIME: %d\t%.3f\tEpoch_Loss\n", OBSTIME->Week, OBSTIME->SecOfWeek);
		return 0;
	case Set_UP_B_fail:
		fprintf(fpr, "GPSTIME: %d\t%.3f\tSet_UP_B_fail\n", OBSTIME->Week, OBSTIME->SecOfWeek);
		return 0;
	default:
		break;
	}
	return 0;
}

void DATA_SET::KF_Print(FILE* fpr, Configure cfg)
{
	XYZ xyz = get_XYZ(KF->getState().block(0, 0, 3, 1));
	BLH blh = XYZ2BLH(xyz, WGS84_e2, WGS84_a);
	XYZ enu = XYZ2ENU(*Real_Pos, xyz, SYS_GPS);
	double dt_G = KF->getState()(3, 0);
	double dt_C = 0;
	XYZ Vel;
	double Rcv_t_v = 0;
	if (cfg.SYS_num == 1)
	{
		if (cfg.GPS_Cfg.used)
			dt_G = KF->getState()(3, 0);
		if (cfg.BDS_Cfg.used)
			dt_C = KF->getState()(3, 0);
		Vel = get_XYZ(KF->getState().block(4, 0, 3, 1));
		Rcv_t_v = KF->getState()(7, 0);
	}
	else if (cfg.SYS_num == 2)
	{
		dt_G = KF->getState()(3, 0);
		dt_C = KF->getState()(4, 0);
		Vel = get_XYZ(KF->getState().block(5, 0, 3, 1));
		Rcv_t_v = KF->getState()(8, 0);
	}
	switch (KF_result)
	{
	case UN_Solve:
		printf("GPSTIME: %d\t%.3f\tUN_Solve\n", OBSTIME->Week, OBSTIME->SecOfWeek);
		fprintf(fpr, "GPSTIME: %d\t%.3f\tUN_Solve\n", OBSTIME->Week, OBSTIME->SecOfWeek);
		break;
	case Success_Solve:
		printf("GPSTIME: %d\t%.3f\tXYZ: %.4f\t%.4f\t%.4f\tBLH: %8.4f\t%8.4f\t%7.4f\tENU: %7.4f\t%7.4f\t%7.4f\t",
			OBSTIME->Week, OBSTIME->SecOfWeek,
			xyz.X, xyz.Y, xyz.Z,
			blh.Lat, blh.Lon, blh.Height,
			enu.X, enu.Y, enu.Z);
		if (cfg.GPS_Cfg.used)
			printf("GPS Clk : % 7.4f\t", dt_G);
		if (cfg.BDS_Cfg.used)
			printf("BDS Clk : % 7.4f\t", dt_C);
		printf("Velocity : % 8.4f\t % 8.4f\t % 8.4f\t % 8.4f\tGPS : % 2d\tBDS : % 2d\t % s\n",
			Vel.X, Vel.Y, Vel.Z, Rcv_t_v,
			KF_GPS_num, KF_BDS_num,
			KF_SATES->c_str());

		fprintf(fpr, "GPSTIME: %d\t%.3f\tXYZ: %.4f\t%.4f\t%.4f\tBLH: %8.4f\t%8.4f\t%7.4f\tENU: %7.4f\t%7.4f\t%7.4f\t",
			OBSTIME->Week, OBSTIME->SecOfWeek,
			xyz.X, xyz.Y, xyz.Z,
			blh.Lat, blh.Lon, blh.Height,
			enu.X, enu.Y, enu.Z);
		if (cfg.GPS_Cfg.used)
			fprintf(fpr, "GPS Clk : % 7.4f\t", dt_G);
		if (cfg.BDS_Cfg.used)
			fprintf(fpr, "BDS Clk : % 7.4f\t", dt_C);
		fprintf(fpr, "Velocity : % 8.4f\t % 8.4f\t % 8.4f\t % 8.4f\tGPS : % 2d\tBDS : % 2d\t % s\n",
			Vel.X, Vel.Y, Vel.Z, Rcv_t_v,
			KF_GPS_num, KF_BDS_num,
			KF_SATES->c_str());
		break;
	case OBS_DATA_Loss:
		printf("GPSTIME: %d\t%.3f\tOBS_DATA_Loss\n", OBSTIME->Week, OBSTIME->SecOfWeek);
		fprintf(fpr, "GPSTIME: %d\t%.3f\tOBS_DATA_Loss\n", OBSTIME->Week, OBSTIME->SecOfWeek);
	case Epoch_Loss:
		printf("GPSTIME: %d\t%.3f\tEpoch_Loss\n", OBSTIME->Week, OBSTIME->SecOfWeek);
		fprintf(fpr, "GPSTIME: %d\t%.3f\tEpoch_Loss\n", OBSTIME->Week, OBSTIME->SecOfWeek);
	case Set_UP_B_fail:
		printf("GPSTIME: %d\t%.3f\tSet_UP_B_fail\n", OBSTIME->Week, OBSTIME->SecOfWeek);
		fprintf(fpr, "GPSTIME: %d\t%.3f\tSet_UP_B_fail\n", OBSTIME->Week, OBSTIME->SecOfWeek);
	default:
		break;
	}
}


int createDirectory(string path)
{
	int len = path.length();
	char tmpDirPath[256] = { 0 };
	for (int i = 0; i < len; i++)
	{
		tmpDirPath[i] = path[i];
		if (tmpDirPath[i] == '\\' || tmpDirPath[i] == '/')
		{
			if (_access(tmpDirPath, 0) == -1)
			{
				int ret = _mkdir(tmpDirPath);
				if (ret == -1)
					return ret;
			}
		}
	}
	return 0;
}

int DATA_SET::CheckOBSConsist(Satellate* sate, int sys, double t, int index, bool& PSE_flag, bool& PHA_flag)
{
	//函数可以优化
	int prn = sate->PRN;
	int Index = decode_SYN(sys, sate->SYG_TYPE[index]);
	double f = CODE2FREQ(Index);
	double lamda = 1e-6 * velocity_c / f;
	double C_D = 0;
	double mean_D = 0;
	double C_P = 0;
	double C_L = 0;
	if (sate->DOPPLER[0] == 0)
	{
		cout << "DOPPLER DATA LOSS!" << endl;
		return 0;
	}
	switch (sys)
	{
	case SYS_GPS:
		Index -= 1;
		if (GPS_PHA[Index][prn - 1] == 0 && GPS_PSE[Index][prn - 1] == 0 && GPS_DOP[Index][prn - 1] == 0)
		{
			GPS_PHA[Index][prn - 1] = sate->PHASE[index];
			GPS_PSE[Index][prn - 1] = sate->PSERA[index];
			GPS_DOP[Index][prn - 1] = sate->DOPPLER[index];
			return 1;
		}
		C_D = lamda * abs(GPS_DOP[Index][prn - 1] - sate->DOPPLER[index]);
		if (C_D > 20)
			return -1;
		mean_D = (GPS_DOP[Index][prn - 1] + sate->DOPPLER[index]) / 2;
		C_P = abs((sate->PSERA[index] - GPS_PSE[Index][prn - 1]) + lamda * mean_D * t);
		C_L = abs(lamda * (sate->PHASE[index] - GPS_PHA[Index][prn - 1]) + lamda * mean_D * t);
		if (C_P > 8)
			PSE_flag = false;
		if (C_L > 0.5)
			PHA_flag = false;
		GPS_PSE[Index][prn - 1] = sate->PSERA[index];
		GPS_PHA[Index][prn - 1] = sate->PHASE[index];
		GPS_DOP[Index][prn - 1] = sate->DOPPLER[index];

		return 1;
	case SYS_BDS:
		Index -= 7;
		if (BDS_PHA[Index][prn - 1] == 0 && BDS_PSE[Index][prn - 1] == 0 && BDS_DOP[Index][prn - 1] == 0)
		{
			BDS_PHA[Index][prn - 1] = sate->PHASE[index];
			BDS_PSE[Index][prn - 1] = sate->PSERA[index];
			BDS_DOP[Index][prn - 1] = sate->DOPPLER[index];
			return 1;
		}
		C_D = lamda * abs(BDS_DOP[Index][prn - 1] - sate->DOPPLER[index]);
		if (C_D > 20)
			return -1;
		mean_D = (BDS_DOP[Index][prn - 1] + sate->DOPPLER[index]) / 2;
		C_P = abs((sate->PSERA[index] - BDS_PSE[Index][prn - 1]) + lamda * mean_D * t);
		C_L = abs(lamda * (sate->PHASE[index] - BDS_PHA[Index][prn - 1]) + lamda * mean_D * t);
		if (C_P > 8)
			PSE_flag = false;
		if (C_L > 0.5)
			PHA_flag = false;
		BDS_PSE[Index][prn - 1] = sate->PSERA[index];
		BDS_PHA[Index][prn - 1] = sate->PHASE[index];
		BDS_DOP[Index][prn - 1] = sate->DOPPLER[index];

		return 1;
	default:
		return 0;
		break;
	}
}

int DATA_SET::DetectOutlier(Satellate* sate, int sys, double t, int index1, int index2)
{
	int prn = sate->PRN;
	double f1 = CODE2FREQ(decode_SYN(sys, sate->SYG_TYPE[index1]));
	double f2 = CODE2FREQ(decode_SYN(sys, sate->SYG_TYPE[index2]));
	double lamda1 = 1e-6 * velocity_c / f1;
	double lamda2 = 1e-6 * velocity_c / f2;
	bool PSE_flag1 = true;
	bool PHA_flag1 = true;
	bool PSE_flag2 = true;
	bool PHA_flag2 = true;
	double GF = 0;
	double MW = 0;
	double dGF = 0;
	double dMW = 0;
	switch (sys)
	{
	case SYS_GPS:
		if (CheckOBSConsist(sate, sys, t, index1, PSE_flag1, PHA_flag1) && CheckOBSConsist(sate, sys, t, index2, PSE_flag2, PHA_flag2) && PSE_flag1 && PHA_flag1 && PSE_flag2 && PHA_flag2)
		{
			GF = sate->PSERA[index1] - sate->PSERA[index2];
			MW = (f1 - f2) * (sate->PSERA[index1] / lamda1 + sate->PSERA[index2] / lamda2) / (f1 + f2) - (sate->PHASE[index1] - sate->PHASE[index2]);
			if (GPS_GF[prn - 1] == 0 && GPS_MW[prn - 1] == 0)
			{
				GPS_GF[prn - 1] = GF;
				GPS_MW[prn - 1] = MW;
				GPS_COUNT[prn - 1]++;
				return 1;
			}
			dGF = abs(GF - GPS_GF[prn - 1]);
			dMW = abs(MW - GPS_MW[prn - 1]);
			GPS_GF[prn - 1] = GF;
			GPS_MW[prn - 1] = (GPS_MW[prn - 1] * (GPS_COUNT[prn - 1]++) + MW);
			GPS_MW[prn - 1] /= GPS_COUNT[prn - 1];
			if (dGF > GF_THRESH || dMW > MW_THRESH)
				return 0;

			return 1;
		}
		else
		{
			return 0;
		}
	case SYS_BDS:
		if (CheckOBSConsist(sate, sys, t, index1, PSE_flag1, PHA_flag1) && CheckOBSConsist(sate, sys, t, index2, PSE_flag2, PHA_flag2) && PSE_flag1 && PHA_flag1 && PSE_flag2 && PHA_flag2)
		{
			GF = sate->PSERA[index1] - sate->PSERA[index2];
			MW = (f1 - f2) * (sate->PSERA[index1] / lamda1 + sate->PSERA[index2] / lamda2) / (f1 + f2) - (sate->PHASE[index1] - sate->PHASE[index2]);
			if (BDS_GF[prn - 1] == 0 && BDS_MW[prn - 1] == 0)
			{
				BDS_GF[prn - 1] = GF;
				BDS_MW[prn - 1] = MW;
				BDS_COUNT[prn - 1]++;
				return 1;
			}
			dGF = abs(GF - BDS_GF[prn - 1]);
			dMW = abs(MW - BDS_MW[prn - 1]);
			BDS_GF[prn - 1] = GF;
			BDS_MW[prn - 1] = (BDS_MW[prn - 1] * (BDS_COUNT[prn - 1]++) + MW);
			BDS_MW[prn - 1] /= BDS_COUNT[prn - 1];
			if (dGF > GF_THRESH || dMW > MW_THRESH)
				return 0;
			return 1;
		}
		else
		{
			return 0;
		}
	default:
		return 0;
		break;
	}
}

void DATA_SET::DetectOut(Configure cfg, double dt_e)
{
	for (int i = 0; i < range->GPS_SATE.size(); i++)
	{
		if (range->GPS_SATE[i]->Phase_NUM < 2)					//需修改
			continue;
		int prn = range->GPS_SATE[i]->PRN;
		int Index1 = 0;
		int Index2 = 0;
		while (CODE2FREQ(decode_SYN(range->GPS_SATE[i]->SYS, range->GPS_SATE[i]->SYG_TYPE[Index1])) != cfg.GPS_Cfg.f1 && Index1 < MAXNUM)
			Index1++;
		while (CODE2FREQ(decode_SYN(range->GPS_SATE[i]->SYS, range->GPS_SATE[i]->SYG_TYPE[Index2])) != cfg.GPS_Cfg.f2 && Index2 < MAXNUM)
			Index2++;

		if (!DetectOutlier(range->GPS_SATE[i], range->GPS_SATE[i]->SYS, dt_e, Index1, Index2) || Index1 == MAXNUM || Index2 == MAXNUM)
		{
			range->GPS_SATE[i]->Outlier = true;
			continue;
		}
	}
	for (int i = 0; i < range->BDS_SATE.size(); i++)
	{
		if (range->BDS_SATE[i]->Phase_NUM < 2)
			continue;
		int prn = range->BDS_SATE[i]->PRN;
		int Index1 = 0;
		int Index2 = 0;
		while (CODE2FREQ(decode_SYN(range->BDS_SATE[i]->SYS, range->BDS_SATE[i]->SYG_TYPE[Index1])) != cfg.BDS_Cfg.f1 && Index1 < MAXNUM)
			Index1++;
		while (CODE2FREQ(decode_SYN(range->BDS_SATE[i]->SYS, range->BDS_SATE[i]->SYG_TYPE[Index2])) != cfg.BDS_Cfg.f2 && Index2 < MAXNUM)
			Index2++;

		if (!DetectOutlier(range->BDS_SATE[i], range->BDS_SATE[i]->SYS, dt_e, Index1, Index2) || Index1 == MAXNUM || Index2 == MAXNUM)
		{
			range->BDS_SATE[i]->Outlier = true;
			continue;
		}
	}
}