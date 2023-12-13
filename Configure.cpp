#include"Configure.h"

Configure::Configure()
{
	Load_cfg();
}
void Configure::Load_cfg()
{
	NetIP = "8.140.46.126";
	NetPort = 5002;

	phase_num = 2;
	SYS_num = 1;
	Hop_used = 1;

	GPS_Cfg = Sate_Configure(L1, L2);
	BDS_Cfg = Sate_Configure(B1, B3);
	GPS_Cfg.used = true;
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