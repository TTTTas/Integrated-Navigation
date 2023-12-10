#pragma once
#include "Common_value.h"

/*���ںʹ洢�ṹ�������*/
class Configure
{
public:
	const char* NetIP;			 // IP
	unsigned short NetPort; // �˿�
	const char* ObsDatFile;				 // log�ļ�·��
	const char* ResDatFile;				 // pos�ļ�·��
	/*Ƶ������ϵͳ������*/
	int phase_num; // ��Ƶor˫Ƶ
	int SYS_num;   // ����or˫��
	int User_SYS;  // ������ָ��ϵͳ��Ĭ��GPS
	int Hop_used;  // �Ƿ����ö��������

	/*GPSϵͳ����*/
	double GPS_f1;
	double GPS_f2;

	/*BDSϵͳ����*/
	double BDS_f1;
	double BDS_f2;

	void Load_cfg();
};