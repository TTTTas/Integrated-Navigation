#include <iostream>
#include <Eigen/Dense>
#include <vector>
using namespace std;
using namespace Eigen;

/*地球椭球相关参数*/
#define WGS84_e2 0.0066943799013
#define WGS84_a 6378137.0
#define CGCS2000_e2 0.00669438002290
#define CGCS2000_a 6378137.0
#define Pi 3.1415926
/*system IDs*/
#define SYS_GPS         0
#define SYS_BDS         4

#pragma region 时间域
/*世界时*/
typedef struct UTC
{
	unsigned short Year;		
	unsigned short Month;
	unsigned short Day;
	unsigned short Hour;
	unsigned short Min;
	unsigned short Sec;

	UTC(unsigned short y, unsigned short m, unsigned short d, unsigned short h, unsigned short min, unsigned short s)
	{
		Year = y;
		Month = m;
		Day = d;
		Hour = h;
		Min = min;
		Sec = s;
	}
	UTC()
	{
		Year = Month = Day = Hour = Min = Sec = 0;
	}
};
/*简化儒略日*/
typedef struct MJD
{
	int Days;
	double FracDay;

	MJD()
	{
		Days = 0;
		FracDay = 0.0;
	}
};
/*GPS时*/
typedef struct GPSTIME
{
	unsigned short Week;
	double SecOfWeek;

	GPSTIME()
	{
		Week = 0;
		SecOfWeek = 0.0;
	}
	GPSTIME(unsigned short w, double s)
	{
		Week = w;
		SecOfWeek = s;
	}
};
/*时间相互转化函数*/
MJD UTC2MJD(UTC utc);

UTC MJD2UTC(MJD mjd);

GPSTIME MJD2GPSTIME(MJD mjd);

MJD GPSTIME2MJD(GPSTIME gpst);

GPSTIME GPSTIME2BDSTIME(GPSTIME gpst);

GPSTIME BDSTIME2GPSTIME(GPSTIME bdst);

UTC GPSTIME2UTC(GPSTIME gpst);


#pragma endregion

#pragma region 坐标域
/*XYZ坐标系*/
typedef struct XYZ
{
	double X;		//m
	double Y;		//m
	double Z;		//m
	XYZ()
	{
		X = Y = Z = 0;
	}
	XYZ(double x, double y, double z)
	{
		X = x;
		Y = y;
		Z = z;
	}
};
/*大地坐标系*/
typedef struct BLH
{
	double Lon;				//经度,deg
	double Lat;				//纬度,deg
	double Height;			//高程,m

	BLH()
	{
		Lon = Lat = Height = 0;
	}
	BLH(double lon, double lat, double h)
	{
		Lon = lon;
		Lat = lat;
		Height = h;
	}
};
/*坐标系转换*/
XYZ BLH2XYZ(BLH blh, double e2, double a);

BLH XYZ2BLH(XYZ xyz, double e2, double a);

//转化为当地坐标系
/*****************************/
//	xyz1: 站心在xyz坐标系下坐标
//	xyz2: 目标在xyz坐标系下坐标
XYZ XYZ2ENU(XYZ xyz1, XYZ xyz2, int sys);

MatrixXd get_Rot(double B, double L);

MatrixXd get_XYZ(XYZ xyz);

MatrixXd get_BLH(BLH blh);

XYZ get_XYZ(MatrixXd xyz);

BLH get_BLH(MatrixXd blh);
#pragma endregion


#pragma once
