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
	/* ���� */
	Configure CfgInfo;
	
    /* ����洢���� */
	DATA_SET* result = new DATA_SET(CfgInfo);
	
    /* ��ʼ�жϱ��� */
	double dt_epoch = 1; // �ļ�����Ԫ��ʱ���
	double temp_t = 0;
	FILE* DATA_Fobs;   // log�ļ�ָ��
	FILE* Pos_Fobs;    // pos�ļ�ָ��
	FILE* KF_Fobs;
	
    /* ��������������ر��� */
	int lenR, lenD;
	unsigned char curbuff[MAXRAWLEN];
	lenD = 0;
	unsigned char decBuff[2 * MAXRAWLEN];
	SOCKET NetGps;
	int ReadFlag;
	int choice = 0;

	printf("��ѡ�����뷽ʽ\n1. �ļ�\t2. ����\n");
	cin >> choice;
	
    /* ��ȡ�ļ�����ʱ�� */
	time_t nowtime;
	time(&nowtime); // ��ȡ1970��1��1��0��0��0�뵽���ھ���������
	tm p;
	localtime_s(&p, &nowtime); // ������ת��Ϊ����ʱ��,���1900����,��Ҫ+1900,��Ϊ0-11,����Ҫ+1
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
		cout << "�������ļ�·��" << endl;
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

			if ((lenR = recv(NetGps, (char*)curbuff, MAXRAWLEN, 0)) > 0) // ��ȡ����
			{
				printf("%5d\n", lenR);
				fwrite(curbuff, sizeof(unsigned char), lenR, DATA_Fobs); // ��¼���������������ļ���

				if ((lenD + lenR) > 2 * MAXRAWLEN)
					lenD = 0;

				memcpy(decBuff + lenD, curbuff, lenR); // ����ƴ��
				lenD += lenR;

				ReadFlag = decodestream(result, decBuff, lenD); // ����

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
						result->LS_print(CfgInfo);              // ���������̨
						result->LS_Filewrite(Pos_Fobs, CfgInfo); // ������ļ�
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