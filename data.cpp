#include "data.h"
#include "cal.h"
#include <io.h>
#include <direct.h>

Result_DATA::Result_DATA()
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

	KF = new KalmanFilter(MatrixXd::Zero(8, 1), MatrixXd::Identity(8, 8));
	Real_Pos = new XYZ(-2267807.853, 5009320.431, 3221020.875);
}

void Result_DATA::reset()
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

int Result_DATA::OUTPUT()
{
	BLH* blh = new BLH();
	XYZ* xyz = new XYZ();
	XYZ* enu = new XYZ();
	MatrixXd R(3, 3);
	MatrixXd Q_local;
	double m_H, m_V, B, L;
	switch (solve_result)
	{
	case Success:
		xyz->X = (*Pos)(0, 0);
		xyz->Y = (*Pos)(1, 0);
		xyz->Z = (*Pos)(2, 0);
		XYZ2BLH(blh, xyz, WGS84_e2, WGS84_a);
		XYZ2ENU(Real_Pos, xyz, enu, SYS_GPS);
		B = degree2rad(blh->Lat);
		L = degree2rad(blh->Lon);
		R(0, 0) = -sin(L);
		R(0, 1) = cos(L);
		R(0, 2) = 0;
		R(1, 0) = -sin(B) * cos(L);
		R(1, 1) = -sin(B) * sin(L);
		R(1, 2) = cos(B);
		R(2, 0) = cos(B) * cos(L);
		R(2, 1) = cos(B) * sin(L);
		R(2, 2) = sin(B);
		Q_local = R * ((*Q_Pos).block(0, 0, 3, 3)) * R.transpose();
		m_H = (*thegma_Pos) * sqrt(Q_local(2, 2));
		m_V = (*thegma_Pos) * sqrt(Q_local(0, 0) + Q_local(1, 1));
		printf("GPSTIME: %d\t%.3f\tXYZ: %.4f\t%.4f\t%.4f\tBLH: %8.4f\t%8.4f\t%7.4f\tENU: %7.4f\t%7.4f\t%7.4f\tGPS Clk: %7.4f\tBDS Clk: %7.4f\tm_H: %7.4f\tm_V: %7.4f\tVelocity: %8.4f\t%8.4f\t%8.4f\t%8.4f\tthegma_P: %6.4f\tthegma_V: %6.4f\tPDOP: %6.4f\tGPS: %2d\tBDS: %2d\t%s\n",
			OBSTIME->Week, OBSTIME->SecOfWeek,
			(*Pos)(0, 0), (*Pos)(1, 0), (*Pos)(2, 0),
			blh->Lat, blh->Lon, blh->Height,
			enu->X, enu->Y, enu->Z,
			(*Pos)(3, 0), (*Pos)(4, 0),
			m_H, m_V,
			(*Vel)(0, 0), (*Vel)(1, 0), (*Vel)(2, 0), (*Vel)(3, 0),
			*thegma_Pos, *thegma_Vel, *PDOP,
			GPS_num, BDS_num,
			SATES->c_str());
		delete blh;
		delete xyz;
		delete enu;
		return 1;
		break;
	case OBS_DATA_Loss:
		printf("GPSTIME: %d\t%.3f\tOBS_DATA_Loss\n", OBSTIME->Week, OBSTIME->SecOfWeek);
		delete blh;
		delete xyz;
		delete enu;
		return 0;
	case Epoch_Loss:
		printf("GPSTIME: %d\t%.3f\tEpoch_Loss\n", OBSTIME->Week, OBSTIME->SecOfWeek);
		delete blh;
		delete xyz;
		delete enu;
		return 0;
	case Set_UP_B_fail:
		printf("GPSTIME: %d\t%.3f\tSet_UP_B_fail\n", OBSTIME->Week, OBSTIME->SecOfWeek);
		delete blh;
		delete xyz;
		delete enu;
		return 0;
	default:
		break;
	}
	return 0;
}

int Result_DATA::WRITEOUTPUT(FILE* fpr)
{
	BLH* blh = new BLH();
	XYZ* xyz = new XYZ();
	XYZ* enu = new XYZ();
	MatrixXd R(3, 3);
	MatrixXd Q_local;
	double m_H, m_V, B, L;
	switch (solve_result)
	{
	case Success:
		xyz->X = (*Pos)(0, 0);
		xyz->Y = (*Pos)(1, 0);
		xyz->Z = (*Pos)(2, 0);
		XYZ2BLH(blh, xyz, WGS84_e2, WGS84_a);
		XYZ2ENU(Real_Pos, xyz, enu, SYS_GPS);
		B = degree2rad(blh->Lat);
		L = degree2rad(blh->Lon);
		R(0, 0) = -sin(L);
		R(0, 1) = cos(L);
		R(0, 2) = 0;
		R(1, 0) = -sin(B) * cos(L);
		R(1, 1) = -sin(B) * sin(L);
		R(1, 2) = cos(B);
		R(2, 0) = cos(B) * cos(L);
		R(2, 1) = cos(B) * sin(L);
		R(2, 2) = sin(B);
		Q_local = R * ((*Q_Pos).block(0, 0, 3, 3)) * R.transpose();
		m_H = (*thegma_Pos) * sqrt(Q_local(2, 2));
		m_V = (*thegma_Pos) * sqrt(Q_local(0, 0) + Q_local(1, 1));
		fprintf(fpr, "GPSTIME:%d\t%.3f\tXYZ:%.4f\t%.4f\t%.4f\tBLH:%8.4f\t%8.4f\t%7.4f\tENU:%7.4f\t%7.4f\t%7.4f\tGPS Clk:%7.4f\tBDS Clk:%7.4f\tm_H:%7.4f\tm_V:%7.4f\tVelocity:%8.4f\t%8.4f\t%8.4f\t%8.4f\tthegma_P:%6.4f\tthegma_V:%6.4f\tPDOP:%6.4f\tGPS:%2d\tBDS:%2d\t%s\n",
			OBSTIME->Week, OBSTIME->SecOfWeek,
			(*Pos)(0, 0), (*Pos)(1, 0), (*Pos)(2, 0),
			blh->Lat, blh->Lon, blh->Height,
			enu->X, enu->Y, enu->Z,
			(*Pos)(3, 0), (*Pos)(4, 0),
			m_H, m_V,
			(*Vel)(0, 0), (*Vel)(1, 0), (*Vel)(2, 0), (*Vel)(3, 0),
			*thegma_Pos, *thegma_Vel, *PDOP,
			GPS_num, BDS_num,
			SATES->c_str());
		delete blh;
		delete xyz;
		delete enu;
		return 1;
		break;
	case OBS_DATA_Loss:
		fprintf(fpr, "GPSTIME: %d\t%.3f\tOBS_DATA_Loss\n", OBSTIME->Week, OBSTIME->SecOfWeek);
		delete blh;
		delete xyz;
		delete enu;
		return 0;
	case Epoch_Loss:
		fprintf(fpr, "GPSTIME: %d\t%.3f\tEpoch_Loss\n", OBSTIME->Week, OBSTIME->SecOfWeek);
		delete blh;
		delete xyz;
		delete enu;
		return 0;
	case Set_UP_B_fail:
		fprintf(fpr, "GPSTIME: %d\t%.3f\tSet_UP_B_fail\n", OBSTIME->Week, OBSTIME->SecOfWeek);
		delete blh;
		delete xyz;
		delete enu;
		return 0;
	default:
		break;
	}
	return 0;
}

void Result_DATA::KF_Print()
{
	cout << KF->getState().transpose() << endl;
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

int decodestream(Result_DATA* result, unsigned char Buff[], int& d)
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
		/*�ļ�Ԥ����*/
		for (; i < d - 2; i++) // ͬ��
		{
			if (Buff[i] == OEM4SYNC1 && Buff[i + 1] == OEM4SYNC2 && Buff[i + 2] == OEM4SYNC3)
			{
				break;
			}
		}
		key++;
		if (i + OEM4HLEN >= d) // ��Ϣ������������
			break;

		for (j = 0; j < OEM4HLEN; j++)
			TempBuff[j] = Buff[i + j]; // ������Ϣͷ�������뻺����

		len = U2(TempBuff + 8) + OEM4HLEN; // ��Ϣͷ+��Ϣ�峤��
		// if (OEM4HLEN != U1(TempBuff + 3))
		//{
		//	i += len + 4;
		//	continue;
		// }

		if ((len + 4 + i) > d || len > MAXRAWLEN) // ��Ϣ������������
			break;

		for (j = OEM4HLEN; j < len + 4; j++) // ������Ϣ�嵽�����뻺����
			TempBuff[j] = Buff[i + j];

		msgID = U2(TempBuff + 4);

		if (CRC32(TempBuff, len) != UCRC32(TempBuff + len, 4)) // ����CRC32
		{
			i += len + 4;
			continue;
		}

		/*¼���ļ�ͷ��Ϣ*/
		msgTYPE = (U1(TempBuff + 6) >> 4) & 0X3;
		gpstime = new GPSTIME(U2(TempBuff + 14), (double)(U4(TempBuff + 16) * 1e-3));
		if (msgTYPE != 0)
			continue;
		int prn = 0;
		/*����ͬ������Ϣ*/
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

		if (val == 1)  //�۲����ݽ���ɹ�
			break;
	}
	//---------------����󣬻���Ĵ���-------------------//
	for (j = 0; j < d - i; j++)
		Buff[j] = Buff[i + j];

	d = j; // ����󣬻�����ʣ�����δ������ֽ���
	//---------------����󣬻���Ĵ���-------------------//
	return val;
}