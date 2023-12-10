#include"Configure.h"

void Configure::Load_cfg()
{
	NetIP = "8.140.46.126";
	NetPort = 5002;

	phase_num = 2;
	SYS_num = 1;
	User_SYS = SYS_GPS;
	Hop_used = 1;

	GPS_f1 = L1;
	GPS_f2 = L2;

	BDS_f1 = B1;
	BDS_f2 = B3;
}