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

/*解算结果*/
#define UN_Solve 0
#define Success_Solve 1
#define Epoch_Loss -1
#define OBS_DATA_Loss -2
#define Set_UP_B_fail -3

/*解算结果存储类*/
class DATA_SET
{
public:
	GPSTIME* OBSTIME;	// 观测时间
	MatrixXd* Pos;		// 位置+钟差
	MatrixXd* Vel;		// 速度+钟速
	MatrixXd* Q_Pos;	// Pos的Q矩阵
	MatrixXd* Q_Vel;	// Vel的Q矩阵
	double* thegma_Pos; // Pos的单位权中误差
	double* thegma_Vel; // Vel的单位权中误差
	double* PDOP;		// 位置的DOP值
	double* VDOP;		// 速度的DOP值
	int LS_GPS_num;		// GPS卫星数
	int LS_BDS_num;		// BDS卫星数
	int KF_GPS_num;		// GPS卫星数
	int KF_BDS_num;		// BDS卫星数
	string* LS_SATES;		// 使用的卫星集合
	string* KF_SATES;
	int LS_result;	// 解算结果
	int KF_result;
	XYZ* Real_Pos;		// 参考真值

	/*卫星粗差探测存储变量*/
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

	/*观测值数据*/
	OBS_DATA* range;

	/*单历元星历数据存储*/
	EPHEMERIS* GPS_eph[GPS_SAT_QUAN];
	EPHEMERIS* BDS_eph[BDS_SAT_QUAN];

	/*求解过程参量*/
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

	int LS_print(Configure cfg);				// 控制台输出
	int LS_Filewrite(FILE* fpr, Configure cfg); // 文件输出
	void KF_Print(FILE* fpr, Configure cfg);

	int Set_KF(Configure cfg);
	int Set_LS(Configure cfg);

	/*一致性检验*/
	int CheckOBSConsist(Satellate* sate, int sys, double t, int index, bool& PSE_flag, bool& PHA_flag);

	/*粗差探测*/
	int DetectOutlier(Satellate* sate, int sys, double t, int index1, int index2);

	void DetectOut(Configure cfg, double dt_e);
};

/*创建文件夹*/
int createDirectory(string path);
