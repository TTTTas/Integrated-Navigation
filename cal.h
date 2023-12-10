#pragma once
#include "Common_value.h"
#include "read.h"
#include "transform.h"
#include "Configure.h"

using namespace std;

/*平方数*/
double SQR(double x);
/*模长*/
double Len(XYZ *pos);

/*角度单位转换*/
double degree2rad(double degree);

double rad2degree(double rad);
/*计算DOP值*/
double Cal_PDOP(MatrixXd Qxx);
/*最小二乘计算*/
unsigned int Cal_LEAST_SQR(MatrixXd B, MatrixXd l, MatrixXd P, MatrixXd *Qxx, MatrixXd &x, double *thegma, double *DOP);
/*解码频点*/
unsigned int decode_SYN(int sys, int signal);
/*解码频率*/
double CODE2FREQ(int code);

// 钟差改正
double CORRECT_CLK(double t, EPHEMERIS *eph);

// TGD计算
double TGD(EPHEMERIS *e, double f, int sys);

// 星历位置
unsigned int SAT_POS_CAL(double t, EPHEMERIS *eph, XYZ *xyz, double &clk, double dt, int SYS);

// 卫星高度角计算
double Ele_Angle(XYZ SatPos, XYZ RcvPos, int sys);

// Hopefiled对流层改正(m)
double Hopefield(double E, double H);

// Hopefiled对流层改正(m)
double Hopefield(XYZ SatPos, XYZ RcvPos, int sys);

/*Klobuchar模型改正电离层*/
double Klobuchar(XYZ RcvPos, double E, double A, double alpha[4], double beta[4], double UT, double code, int sys);

/*文件流下SPP解算*/
// 搭建位置解算构造矩阵(双频)
unsigned int setup_Pos(GPSTIME *OBS_TIME, MatrixXd Pos, vector<Satellate *> Sates, EPOCH *eph, bool first_flag, double f1, double f2, MatrixXd *B_Pos, MatrixXd *l_Pos, MatrixXd *P_Pos, string *sate_used);

// 搭建位置解算构造矩阵(单频)
unsigned int setup_Pos(GPSTIME *OBS_TIME, MatrixXd Pos, vector<Satellate *> Sates, EPOCH *eph, bool first_flag, double f, MatrixXd *B_Pos, MatrixXd *l_Pos, MatrixXd *P_Pos, string *sate_used);

// 搭建速度解算构造矩阵
unsigned int setup_Vel(GPSTIME *OBS_TIME, MatrixXd Pos, vector<Satellate *> Sates, EPOCH *eph, MatrixXd *B_Vel, MatrixXd *l_Vel, MatrixXd *P_Vel);

// 单星解算
unsigned int Cal_1(DATA_SET *data, OBS_DATA *obs, EPOCH *eph, bool first_flag);

// GPS、BDS双星解算
unsigned int Cal_2(DATA_SET *data, OBS_DATA *obs, EPOCH *gpseph, EPOCH *bdseph, bool first_flag);

// SPP单点定位
unsigned int Cal_SPP(DATA_SET *data, OBS_DATA *obs, EPOCH *gpseph, EPOCH *bdseph, double dt_e, bool first_flag);



/*网口下位置解算*/
// 搭建位置解算构造矩阵(双频)
unsigned int setup_Pos(GPSTIME *OBS_TIME, MatrixXd Pos, vector<Satellate *> Sates, EPHEMERIS **eph, bool first_flag, double f1, double f2, MatrixXd *B_Pos, MatrixXd *l_Pos, MatrixXd *P_Pos, string *sate_used);

unsigned int setup_Pos(DATA_SET* data, Configure cfg);

// 搭建位置解算构造矩阵(单频)
unsigned int setup_Pos(GPSTIME *OBS_TIME, MatrixXd Pos, vector<Satellate *> Sates, EPHEMERIS **eph, bool first_flag, double f, MatrixXd *B_Pos, MatrixXd *l_Pos, MatrixXd *P_Pos, string *sate_used);

// 搭建速度解算构造矩阵
unsigned int setup_Vel(GPSTIME *OBS_TIME, MatrixXd Pos, vector<Satellate *> Sates, EPHEMERIS **eph, MatrixXd *B_Vel, MatrixXd *l_Vel, MatrixXd *P_Vel);

// 单星解算
unsigned int Cal_1(DATA_SET *data, bool first_flag);

// GPS、BDS双星解算
unsigned int Cal_2(DATA_SET *data, bool first_flag);

// SPP单点定位
unsigned int Cal_SPP(DATA_SET *data, double dt_e, bool first_flag);