#pragma once
#include "Common_value.h"

/*网口和存储结构相关配置*/
class Configure
{
public:
	const char* NetIP;			 // IP
	unsigned short NetPort; // 端口
	const char* ObsDatFile;				 // log文件路径
	const char* ResDatFile;				 // pos文件路径
	/*频点数、系统数配置*/
	int phase_num; // 单频or双频
	int SYS_num;   // 单星or双星
	int User_SYS;  // 单星下指定系统，默认GPS
	int Hop_used;  // 是否启用对流层改正

	/*GPS系统配置*/
	double GPS_f1;
	double GPS_f2;

	/*BDS系统配置*/
	double BDS_f1;
	double BDS_f2;

	void Load_cfg();
};