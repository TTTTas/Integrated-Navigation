#pragma once
#include <Eigen/Dense>
#include "transform.h"
#include "KF.h"
#include "LS.h"
#include "read.h"
#include "Configure.h"
#include <string>
#include <vector>
using namespace Eigen;
using namespace std;

#define GPS_SAT_QUAN 32
#define BDS_SAT_QUAN 63

/*������*/
#define UN_Solve 0
#define Success_Solve 1
#define Epoch_Loss -1
#define OBS_DATA_Loss -2
#define Set_UP_B_fail -3

/*�������洢��*/
class DATA_SET
{
public:
	GPSTIME* OBSTIME;	// �۲�ʱ��
	MatrixXd* Pos;		// λ��+�Ӳ�
	MatrixXd* Vel;		// �ٶ�+����
	MatrixXd* Q_Pos;	// Pos��Q����
	MatrixXd* Q_Vel;	// Vel��Q����
	double* thegma_Pos; // Pos�ĵ�λȨ�����
	double* thegma_Vel; // Vel�ĵ�λȨ�����
	double* PDOP;		// λ�õ�DOPֵ
	double* VDOP;		// �ٶȵ�DOPֵ
	int LS_GPS_num;		// GPS������
	int LS_BDS_num;		// BDS������
	int KF_GPS_num;		// GPS������
	int KF_BDS_num;		// BDS������
	string* LS_SATES;		// ʹ�õ����Ǽ���
	string* KF_SATES;
	int LS_result;	// ������
	int KF_result;
	XYZ* Real_Pos;		// �ο���ֵ

	/*���Ǵֲ�̽��洢����*/
	double GPS_GF[GPS_SAT_QUAN];
	double GPS_MW[GPS_SAT_QUAN];
	double GPS_PSE[6][GPS_SAT_QUAN];
	double GPS_PHA[6][GPS_SAT_QUAN];
	double GPS_DOP[6][GPS_SAT_QUAN];
	int GPS_COUNT[GPS_SAT_QUAN];

	double BDS_GF[BDS_SAT_QUAN];
	double BDS_MW[BDS_SAT_QUAN];
	double BDS_PSE[5][BDS_SAT_QUAN];
	double BDS_PHA[5][BDS_SAT_QUAN];
	double BDS_DOP[5][BDS_SAT_QUAN];
	int BDS_COUNT[BDS_SAT_QUAN];

	/*�۲�ֵ����*/
	OBS_DATA* range;

	/*����Ԫ�������ݴ洢*/
	EPHEMERIS* GPS_eph[GPS_SAT_QUAN];
	EPHEMERIS* BDS_eph[BDS_SAT_QUAN];

	/*�����̲���*/
	bool LS_first;
	bool KF_first;
	MatrixXd temp_ref;

	/*Least_Square*/
	Least_Squares* LS_Pos;
	Least_Squares* LS_Vel;

	/*KalmanFilter*/
	KalmanFilter* KF;

	DATA_SET(Configure cfg);

	void reset();

	int LS_print(Configure cfg);				// ����̨���
	int LS_Filewrite(FILE* fpr, Configure cfg); // �ļ����
	void KF_Print(FILE* fpr, Configure cfg);

	int Set_KF(Configure cfg);
	int Set_LS(Configure cfg);

	/*һ���Լ���*/
	int CheckOBSConsist(Satellate* sate, int sys, double t, int index, bool& PSE_flag, bool& PHA_flag);

	/*�ֲ�̽��*/
	int DetectOutlier(Satellate* sate, int sys, double t, int index1, int index2);

	void DetectOut(Configure cfg, double dt_e);
};

/*�����ļ���*/
int createDirectory(string path);
