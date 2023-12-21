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

/*���ںʹ洢�ṹ�������*/
class Configure
{
public:
	const char* NetIP;			 // IP
	unsigned short NetPort; // �˿�
	const char* ObsDatFile;				 // log�ļ�·��
	const char* ResDatFile;				 // pos�ļ�·��
	const char* KFDatFile;
	/*Ƶ������ϵͳ������*/
	int phase_num; // ��Ƶor˫Ƶ
	int SYS_num;   // ����or˫��
	int Hop_used;  // �Ƿ����ö��������

	/*GPSϵͳ����*/
	Sate_Configure GPS_Cfg;

	/*BDSϵͳ����*/
	Sate_Configure BDS_Cfg;

	/*����ģʽ*/
	bool LS_used;
	bool KF_used;

	/*��ֵ*/
	double w_thresh;
	double Ele_Mask;
	Configure();
	void Load_cfg();
};
