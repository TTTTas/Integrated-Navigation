#pragma once

/*FREQUNCY*/ // MHz
/*GPS*/
#define L1 1575.42
#define L2 1227.60
#define L5 1176.45
/*BDS*/
#define B1 1561.098
#define B1_C 1575.42
#define B2 1207.14
#define B2_a 1176.45
#define B3 1268.52

class Sate_Configure
{
public:
	bool used;

	double f1;
	double f2;

	Sate_Configure();
	Sate_Configure(double f1, double f2);
};

/*网口和存储结构相关配置*/
class Configure
{
public:
	const char* NetIP;			 // IP
	unsigned short NetPort; // 端口
	const char* ObsDatFile;				 // log文件路径
	const char* ResDatFile;				 // pos文件路径
	const char* KFDatFile;
	/*频点数、系统数配置*/
	int phase_num; // 单频or双频
	int SYS_num;   // 单星or双星
	int Hop_used;  // 是否启用对流层改正

	/*GPS系统配置*/
	Sate_Configure GPS_Cfg;

	/*BDS系统配置*/
	Sate_Configure BDS_Cfg;

	/*解算模式*/
	bool LS_used;
	bool KF_used;

	/*阈值*/
	double w_thresh;
	double Ele_Mask;
	Configure();
	void Load_cfg();
};
