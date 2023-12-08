#pragma once
#include <Eigen/Dense>
#include "transform.h"
#include "KF.h"
#include "read.h"
#include <string>
#include <vector>
using namespace Eigen;
using namespace std;

#define GPS_SAT_QUAN 32
#define BDS_SAT_QUAN 63

/*解算结果*/
#define UN_Solve 0
#define Success 1
#define Epoch_Loss -1
#define OBS_DATA_Loss -2
#define Set_UP_B_fail -3

/*解算结果存储类*/
class DATA_SET
{
public:
	GPSTIME *OBSTIME;	// 观测时间
	MatrixXd *Pos;		// 位置+钟差
	MatrixXd *Vel;		// 速度+钟速
	MatrixXd *Q_Pos;	// Pos的Q矩阵
	MatrixXd *Q_Vel;	// Vel的Q矩阵
	double *thegma_Pos; // Pos的单位权中误差
	double *thegma_Vel; // Vel的单位权中误差
	double *PDOP;		// 位置的DOP值
	double *VDOP;		// 速度的DOP值
	int GPS_num;		// GPS卫星数
	int BDS_num;		// BDS卫星数
	string *SATES;		// 使用的卫星集合
	int solve_result;	// 解算结果
	XYZ *Real_Pos;		// 参考真值

	/*观测值数据*/
	OBS_DATA *range;

	/*单历元星历数据存储*/
	EPHEMERIS *GPS_eph[GPS_SAT_QUAN];
	EPHEMERIS *BDS_eph[BDS_SAT_QUAN];

	/*KalmanFilter*/
	KalmanFilter *KF;

	DATA_SET();
	void reset();
	int OUTPUT();				// 控制台输出
	int WRITEOUTPUT(FILE *fpr); // 文件输出
	void KF_Print();
	MatrixXd Initial_KF();
};
/*网口和存储结构相关配置*/
class Configure
{
public:
	static const char *NetIP;			 // IP
	static const unsigned short NetPort; // 端口
	const char *ObsDatFile;				 // log文件路径
	const char *ResDatFile;				 // pos文件路径
};

/*创建文件夹*/
int createDirectory(string path);

/*网口下解码*/
int decodestream(DATA_SET *result, unsigned char Buff[], int &d);

// SPP单点定位KF
unsigned int KF_SPP(DATA_SET *data, double dt_e, bool first_flag);

unsigned int KF_1(DATA_SET *data, bool first_flag, double T);
