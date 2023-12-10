#include "data.h"
#include "cal.h"
#include <io.h>
#include <direct.h>

DATA_SET::DATA_SET()
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
	SATES = new string();
	*Pos = MatrixXd::Zero(5, 1);
	*Vel = MatrixXd::Zero(4, 1);
	*Q_Pos = MatrixXd::Zero(5, 5);
	*Q_Vel = MatrixXd::Zero(4, 4);
	*thegma_Pos = 0;
	*thegma_Vel = 0;
	*PDOP = 0;
	*VDOP = 0;
	GPS_num = 0;
	BDS_num = 0;
	*SATES = "";

	range = new OBS_DATA();
	range->OBS_TIME = new GPSTIME();
	for (int i = 0; i < GPS_SAT_QUAN; i++)
	{
		GPS_eph[i] = new EPHEMERIS();
	}
	for (int i = 0; i < BDS_SAT_QUAN; i++)
	{
		BDS_eph[i] = new EPHEMERIS();
	}

	solve_result = UN_Solve;

	MatrixXd P_ = MatrixXd::Zero(8, 8);
	P_.block(0, 0, 4, 4) = MatrixXd::Identity(4, 4) * 10 * 10;
	P_.block(4, 4, 4, 4) = MatrixXd::Identity(4, 4) * 0.1 * 0.1;
	KF = new KalmanFilter(MatrixXd::Zero(8, 1), P_);
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

int DATA_SET::LS_print()
{
	XYZ xyz = get_XYZ(( * Pos).block(0, 0, 3, 1));
	BLH blh = XYZ2BLH(xyz, WGS84_e2, WGS84_a);
	XYZ enu = XYZ2ENU(*Real_Pos, xyz, SYS_GPS);
	MatrixXd R = get_Rot(degree2rad(blh.Lat), degree2rad(blh.Lon));
	MatrixXd Q_local;
	double m_H, m_V, B, L;
	Q_local = R * ((*Q_Pos).block(0, 0, 3, 3)) * R.transpose();
	m_H = (*thegma_Pos) * sqrt(Q_local(2, 2));
	m_V = (*thegma_Pos) * sqrt(Q_local(0, 0) + Q_local(1, 1));
	switch (solve_result)
	{
	case Success:
		printf("GPSTIME: %d\t%.3f\tXYZ: %.4f\t%.4f\t%.4f\tBLH: %8.4f\t%8.4f\t%7.4f\tENU: %7.4f\t%7.4f\t%7.4f\tGPS Clk: %7.4f\tBDS Clk: %7.4f\tm_H: %7.4f\tm_V: %7.4f\tVelocity: %8.4f\t%8.4f\t%8.4f\t%8.4f\tthegma_P: %6.4f\tthegma_V: %6.4f\tPDOP: %6.4f\tGPS: %2d\tBDS: %2d\t%s\n",
			OBSTIME->Week, OBSTIME->SecOfWeek,
			xyz.X,xyz.Y,xyz.Z,
			blh.Lat, blh.Lon, blh.Height,
			enu.X, enu.Y, enu.Z,
			(*Pos)(3, 0), (*Pos)(4, 0),
			m_H, m_V,
			(*Vel)(0, 0), (*Vel)(1, 0), (*Vel)(2, 0), (*Vel)(3, 0),
			*thegma_Pos, *thegma_Vel, *PDOP,
			GPS_num, BDS_num,
			SATES->c_str());
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

int DATA_SET::LS_Filewrite(FILE* fpr)
{
	XYZ xyz = get_XYZ((*Pos).block(0, 0, 3, 1));
	BLH blh = XYZ2BLH(xyz, WGS84_e2, WGS84_a);
	XYZ enu = XYZ2ENU(*Real_Pos, xyz, SYS_GPS);
	MatrixXd R = get_Rot(degree2rad(blh.Lat), degree2rad(blh.Lon));
	MatrixXd Q_local;
	double m_H, m_V;
	Q_local = R * ((*Q_Pos).block(0, 0, 3, 3)) * R.transpose();
	m_H = (*thegma_Pos) * sqrt(Q_local(2, 2));
	m_V = (*thegma_Pos) * sqrt(Q_local(0, 0) + Q_local(1, 1));
	switch (solve_result)
	{
	case Success:
		fprintf(fpr, "GPSTIME: %d\t%.3f\tXYZ: %.4f\t%.4f\t%.4f\tBLH: %8.4f\t%8.4f\t%7.4f\tENU: %7.4f\t%7.4f\t%7.4f\tGPS Clk: %7.4f\tBDS Clk: %7.4f\tm_H: %7.4f\tm_V: %7.4f\tVelocity: %8.4f\t%8.4f\t%8.4f\t%8.4f\tthegma_P: %6.4f\tthegma_V: %6.4f\tPDOP: %6.4f\tGPS: %2d\tBDS: %2d\t%s\n",
			OBSTIME->Week, OBSTIME->SecOfWeek,
			xyz.X, xyz.Y, xyz.Z,
			blh.Lat, blh.Lon, blh.Height,
			enu.X, enu.Y, enu.Z,
			(*Pos)(3, 0), (*Pos)(4, 0),
			m_H, m_V,
			(*Vel)(0, 0), (*Vel)(1, 0), (*Vel)(2, 0), (*Vel)(3, 0),
			*thegma_Pos, *thegma_Vel, *PDOP,
			GPS_num, BDS_num,
			SATES->c_str());
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

void DATA_SET::KF_Print(FILE* fpr)
{	
	XYZ xyz = get_XYZ(KF->getState().block(0, 0, 3, 1));
	BLH blh = XYZ2BLH(xyz, WGS84_e2, WGS84_a);
	XYZ enu = XYZ2ENU(*Real_Pos, xyz, SYS_GPS);
	double Rcv_t = KF->getState()(3, 0);
	XYZ Vel = get_XYZ(KF->getState().block(4, 0, 3, 1));
	double Rcv_t_v = KF->getState()(7, 0);
	printf("GPSTIME: %d\t%.3f\tXYZ: %.4f\t%.4f\t%.4f\tBLH: %8.4f\t%8.4f\t%7.4f\tENU: %7.4f\t%7.4f\t%7.4f\tGPS Clk: %7.4f\tVelocity: %8.4f\t%8.4f\t%8.4f\t%8.4f\n",
		OBSTIME->Week, OBSTIME->SecOfWeek,
		xyz.X, xyz.Y, xyz.Z,
		blh.Lat, blh.Lon, blh.Height,
		enu.X, enu.Y, enu.Z,
		Rcv_t,
		Vel.X, Vel.Y, Vel.Z,
		Rcv_t_v);
	fprintf(fpr, "GPSTIME: %d\t%.3f\tXYZ: %.4f\t%.4f\t%.4f\tBLH: %8.4f\t%8.4f\t%7.4f\tENU: %7.4f\t%7.4f\t%7.4f\tGPS Clk: %7.4f\tVelocity: %8.4f\t%8.4f\t%8.4f\t%8.4f\n",
		OBSTIME->Week, OBSTIME->SecOfWeek,
		xyz.X, xyz.Y, xyz.Z,
		blh.Lat, blh.Lon, blh.Height,
		enu.X, enu.Y, enu.Z,
		Rcv_t,
		Vel.X, Vel.Y, Vel.Z,
		Rcv_t_v);
}

const char* Configure::NetIP = "8.140.46.126";
const unsigned short Configure::NetPort = 5002;

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

int decodestream(DATA_SET* result, unsigned char Buff[], int& d)
{
	unsigned char TempBuff[MAXRAWLEN];
	int len;
	int msgID, msgTYPE;
	GPSTIME* gpstime;
	int key = 0;
	int i, j;
	int val;
	i = 0;
	val = 0;
	while (1)
	{
		/*文件预处理*/
		for (; i < d - 2; i++) // 同步
		{
			if (Buff[i] == OEM4SYNC1 && Buff[i + 1] == OEM4SYNC2 && Buff[i + 2] == OEM4SYNC3)
			{
				break;
			}
		}
		key++;
		if (i + OEM4HLEN >= d) // 消息不完整，跳出
			break;

		for (j = 0; j < OEM4HLEN; j++)
			TempBuff[j] = Buff[i + j]; // 拷贝消息头到待解码缓存中

		len = U2(TempBuff + 8) + OEM4HLEN; // 消息头+消息体长度
		// if (OEM4HLEN != U1(TempBuff + 3))
		//{
		//	i += len + 4;
		//	continue;
		// }

		if ((len + 4 + i) > d || len > MAXRAWLEN) // 消息不完整，跳出
			break;

		for (j = OEM4HLEN; j < len + 4; j++) // 拷贝消息体到待解码缓存中
			TempBuff[j] = Buff[i + j];

		msgID = U2(TempBuff + 4);

		if (CRC32(TempBuff, len) != UCRC32(TempBuff + len, 4)) // 检验CRC32
		{
			i += len + 4;
			continue;
		}

		/*录入文件头信息*/
		msgTYPE = (U1(TempBuff + 6) >> 4) & 0X3;
		gpstime = new GPSTIME(U2(TempBuff + 14), (double)(U4(TempBuff + 16) * 1e-3));
		if (msgTYPE != 0)
			continue;
		int prn = 0;
		/*处理不同数据信息*/
		switch (msgID)
		{
		case ID_RANGE:
			result->range->OBS_TIME->Week = gpstime->Week;
			result->range->OBS_TIME->SecOfWeek = gpstime->SecOfWeek;
			result->range->Sate_Num = U4(TempBuff + OEM4HLEN);
			val = decode_RANGE(TempBuff + OEM4HLEN + 4, result->range->Sate_Num, result->range);
			break;
		case ID_GPSEPHEMERIS:
			prn = U4(TempBuff + OEM4HLEN);
			decode_GPSEPH_STAT(TempBuff + OEM4HLEN, result->GPS_eph[prn - 1]);
			break;
		case ID_BDSEPHEMERIS:
			prn = U4(TempBuff + OEM4HLEN);
			decode_BDSEPH_STAT(TempBuff + OEM4HLEN, result->BDS_eph[prn - 1]);
			break;
		default:
			break;
		}
		delete gpstime;
		i += len + 4;

		if (val == 1)  //观测数据解码成功
			break;
	}
	//---------------解码后，缓存的处理-------------------//
	for (j = 0; j < d - i; j++)
		Buff[j] = Buff[i + j];

	d = j; // 解码后，缓存中剩余的尚未解码的字节数
	//---------------解码后，缓存的处理-------------------//
	return val;
}