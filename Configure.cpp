#include"Configure.h"
#include<iostream>

Configure::Configure()
{
	Load_cfg();
	int count = 0;
	if (GPS_Cfg.used)
		count++;
	if (BDS_Cfg.used)
		count++;
	SYS_num = count;
	std::cout << "系统数: " << SYS_num << "\n"
		<< "频点数: " << phase_num << "\n"
		<< std::endl;
}

void Configure::Load_cfg()
{
	NetIP = "8.140.46.126";
	NetPort = 5002;

	phase_num = 2;
	SYS_num = 0;
	Hop_used = 1;

	GPS_Cfg = Sate_Configure(L1, L2);
	BDS_Cfg = Sate_Configure(B3, B1);
	GPS_Cfg.used = false;
	BDS_Cfg.used = true;

	LS_used = true;
	KF_used = true;

	w_thresh = 10;
	Ele_Mask = 10;
}

Sate_Configure::Sate_Configure(double f1_, double f2_)
{
	used = false;
	f1 = f1_;
	f2 = f2_;
}

Sate_Configure::Sate_Configure()
{
	used = false;
	f1 = 0;
	f2 = 0;
}