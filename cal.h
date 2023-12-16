#pragma once
#include "read.h"
#include "transform.h"
#include "Configure.h"

using namespace std;

/*地球椭球相关参数*/
#define WGS84_e2 0.0066943799013
#define WGS84_a 6378137.0
#define CGCS2000_e2 0.00669438002290
#define CGCS2000_a 6378137.0
#define Pi 3.1415926

/*code IDs*/
#define UNKOWN 0
/*GPS*/
#define CODE_L1C 1
#define CODE_L2P 2
#define CODE_L2W 3
#define CODE_L5Q 4
#define CODE_L1L 5
#define CODE_L2S 6
/*BDS*/
#define CODE_L2I 7
#define CODE_L7I 8
#define CODE_L6I 9
#define CODE_L1P 10
#define CODE_L5P 11

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

/*Hopefiled*/
#define H0 0
#define T0 288.16
#define P0 1013.25
#define RH0 0.5

/*粗差探测阈值*/
#define GF_THRESH 0.05
#define MW_THRESH 10

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

/*网口下位置解算*/
// 搭建位置解算构造矩阵(双频)
unsigned int setup_Pos(GPSTIME *OBS_TIME, MatrixXd Pos, vector<Satellate *> Sates, EPHEMERIS **eph, bool first_flag, double f1, double f2, MatrixXd *B_Pos, MatrixXd *l_Pos, MatrixXd *P_Pos, string *sate_used);

// 搭建位置解算构造矩阵(单频)
unsigned int setup_Pos(GPSTIME *OBS_TIME, MatrixXd Pos, vector<Satellate *> Sates, EPHEMERIS **eph, bool first_flag, double f, MatrixXd *B_Pos, MatrixXd *l_Pos, MatrixXd *P_Pos, string *sate_used);

// 搭建速度解算构造矩阵
unsigned int setup_Vel(GPSTIME *OBS_TIME, MatrixXd Pos, vector<Satellate *> Sates, EPHEMERIS **eph, MatrixXd *B_Vel, MatrixXd *l_Vel, MatrixXd *P_Vel);

unsigned int LS_SPV(DATA_SET* data, Configure cfg);

unsigned int setup_LS(DATA_SET* data, Configure cfg, int sys);

unsigned int setup_KF(DATA_SET* data, Configure cfg, int sys);

double get_measure(Satellate* sate, Configure cfg, EPHEMERIS* eph);

// SPP单点定位KF
unsigned int KF_SPV(DATA_SET* data, double dt_e, Configure cfg);