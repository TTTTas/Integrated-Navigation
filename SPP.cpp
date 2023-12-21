#include <iostream>
#include <fstream>
#include <vector>
#include <Eigen/Dense>
#include <iomanip>
#include <io.h>
#include <direct.h>

#include "transform.h"
#include "read.h"
#include "cal.h"
#include "sockets.h"
#include "KF.h"
#include "data.h"
#include "Configure.h"

using namespace std;
using namespace Eigen;

#define _CRT_SECURE_NO_WARNINGS

int main()
{
	/* 配置 */
	Configure CfgInfo;
	
    /* 结果存储变量 */
	DATA_SET* result = new DATA_SET(CfgInfo);
	
    /* 初始判断变量 */
	double dt_epoch = 1; // 文件流历元间时间差
	double temp_t = 0;
	FILE* DATA_Fobs;   // log文件指针
	FILE* Pos_Fobs;    // pos文件指针
	FILE* KF_Fobs;
	
    /* 网口输入数据相关变量 */
	int lenR, lenD;
	unsigned char curbuff[MAXRAWLEN];
	lenD = 0;
	unsigned char decBuff[2 * MAXRAWLEN];
	SOCKET NetGps;
	int ReadFlag;
	int choice = 0;

	printf("请选择输入方式\n1. 文件\t2. 网口\n");
	cin >> choice;
	
    /* 获取文件生成时间 */
	time_t nowtime;
	time(&nowtime); // 获取1970年1月1日0点0分0秒到现在经过的秒数
	tm p;
	localtime_s(&p, &nowtime); // 将秒数转换为本地时间,年从1900算起,需要+1900,月为0-11,所以要+1
	string filetime = to_string(p.tm_year + 1900) + "_" + to_string(p.tm_mon + 1) + "_" + to_string(p.tm_mday) + "_" + to_string(p.tm_hour) + "_" + to_string(p.tm_min) + "_" + to_string(p.tm_sec);
	string logpath = "D:\\GitHub\\Integrated-Navigation\\data\\logs\\";
	createDirectory(logpath);
	logpath += filetime + string(".log");
	string pospath = "D:\\GitHub\\Integrated-Navigation\\data\\Pos\\";
	string KFpath = "D:\\GitHub\\Integrated-Navigation\\data\\KF\\";
	createDirectory(pospath);
	createDirectory(KFpath);
	pospath += filetime + string(".pos");
	KFpath += filetime + string(".kf");
	CfgInfo.ObsDatFile = logpath.c_str();
	CfgInfo.ResDatFile = pospath.c_str();
	CfgInfo.KFDatFile = KFpath.c_str();

	FILE* file;
	char load_filepath[200];

	switch (choice)
	{
	case 1:
		cout << "请输入文件路径" << endl;
		cin >> load_filepath;
		decodefile(result, CfgInfo, load_filepath);
		break;

	case 2:
		if (OpenSocket(NetGps, CfgInfo.NetIP, CfgInfo.NetPort) == false)
		{
			printf("The ip %s was not opened\n", CfgInfo.NetIP);
			return 0;
		}

		if ((DATA_Fobs = fopen(CfgInfo.ObsDatFile, "wb")) == NULL)
		{
			printf("The obs file %s was not opened\n", CfgInfo.ObsDatFile);
			exit(0);
		}

		if ((Pos_Fobs = fopen(CfgInfo.ResDatFile, "w")) == NULL)
		{
			printf("The pos file %s was not opened\n", CfgInfo.ResDatFile);
			exit(0);
		}

		if ((KF_Fobs = fopen(CfgInfo.KFDatFile, "w")) == NULL)
		{
			printf("The kf file %s was not opened\n", "KF.pos");
			exit(0);
		}

		while (1)
		{
			Sleep(980);

			if ((lenR = recv(NetGps, (char*)curbuff, MAXRAWLEN, 0)) > 0) // 读取数据
			{
				printf("%5d\n", lenR);
				fwrite(curbuff, sizeof(unsigned char), lenR, DATA_Fobs); // 记录二进制数据流到文件中

				if ((lenD + lenR) > 2 * MAXRAWLEN)
					lenD = 0;

				memcpy(decBuff + lenD, curbuff, lenR); // 缓存拼接
				lenD += lenR;

				ReadFlag = decodestream(result, decBuff, lenD); // 解码

				if (ReadFlag != 1)
				{
					printf("Data acquisition and decode failed \n");
					continue;
				}
				else
				{
					if (result->LS_first)
						dt_epoch = 1;
					else
						dt_epoch = result->OBSTIME->SecOfWeek - temp_t;

					if (dt_epoch == 0)
						break;

					temp_t = result->OBSTIME->SecOfWeek;
					result->DetectOut(CfgInfo, dt_epoch);

					if (CfgInfo.LS_used)
					{
						if (LS_SPV(result, CfgInfo))
							result->LS_first = false;

						cout << "LS" << endl;
						result->LS_print(CfgInfo);              // 输出至控制台
						result->LS_Filewrite(Pos_Fobs, CfgInfo); // 输出至文件
					}

					if (CfgInfo.KF_used)
					{
						if (KF_SPV(result, dt_epoch, CfgInfo))
							result->KF_first = false;

						cout << "KF" << endl;
						result->KF_Print(KF_Fobs, CfgInfo);
					}

					result->reset();
				}
			}
			else
			{
				printf("NO MESSAGES IN!\n");
			}
		}

		break;
	default:
		break;
	}
	std::system("pause");
	return 0;
}