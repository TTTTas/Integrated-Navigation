#pragma once
#include "read.h"
#include "transform.h"
#include "Configure.h"

using namespace std;

/*����������ز���*/
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

/*�ֲ�̽����ֵ*/
#define GF_THRESH 0.05
#define MW_THRESH 10

/*ƽ����*/
double SQR(double x);
/*ģ��*/
double Len(XYZ *pos);

/*�Ƕȵ�λת��*/
double degree2rad(double degree);

double rad2degree(double rad);
/*����DOPֵ*/
double Cal_PDOP(MatrixXd Qxx);
/*��С���˼���*/
unsigned int Cal_LEAST_SQR(MatrixXd B, MatrixXd l, MatrixXd P, MatrixXd *Qxx, MatrixXd &x, double *thegma, double *DOP);
/*����Ƶ��*/
unsigned int decode_SYN(int sys, int signal);
/*����Ƶ��*/
double CODE2FREQ(int code);

// �Ӳ����
double CORRECT_CLK(double t, EPHEMERIS *eph);

// TGD����
double TGD(EPHEMERIS *e, double f, int sys);

// ����λ��
unsigned int SAT_POS_CAL(double t, EPHEMERIS *eph, XYZ *xyz, double &clk, double dt, int SYS);

// ���Ǹ߶ȽǼ���
double Ele_Angle(XYZ SatPos, XYZ RcvPos, int sys);

// Hopefiled���������(m)
double Hopefield(double E, double H);

// Hopefiled���������(m)
double Hopefield(XYZ SatPos, XYZ RcvPos, int sys);

/*Klobucharģ�͸��������*/
double Klobuchar(XYZ RcvPos, double E, double A, double alpha[4], double beta[4], double UT, double code, int sys);

/*�ļ�����SPP����*/
// �λ�ý��㹹�����(˫Ƶ)
unsigned int setup_Pos(GPSTIME *OBS_TIME, MatrixXd Pos, vector<Satellate *> Sates, EPOCH *eph, bool first_flag, double f1, double f2, MatrixXd *B_Pos, MatrixXd *l_Pos, MatrixXd *P_Pos, string *sate_used);

// �λ�ý��㹹�����(��Ƶ)
unsigned int setup_Pos(GPSTIME *OBS_TIME, MatrixXd Pos, vector<Satellate *> Sates, EPOCH *eph, bool first_flag, double f, MatrixXd *B_Pos, MatrixXd *l_Pos, MatrixXd *P_Pos, string *sate_used);

// ��ٶȽ��㹹�����
unsigned int setup_Vel(GPSTIME *OBS_TIME, MatrixXd Pos, vector<Satellate *> Sates, EPOCH *eph, MatrixXd *B_Vel, MatrixXd *l_Vel, MatrixXd *P_Vel);

/*������λ�ý���*/
// �λ�ý��㹹�����(˫Ƶ)
unsigned int setup_Pos(GPSTIME *OBS_TIME, MatrixXd Pos, vector<Satellate *> Sates, EPHEMERIS **eph, bool first_flag, double f1, double f2, MatrixXd *B_Pos, MatrixXd *l_Pos, MatrixXd *P_Pos, string *sate_used);

// �λ�ý��㹹�����(��Ƶ)
unsigned int setup_Pos(GPSTIME *OBS_TIME, MatrixXd Pos, vector<Satellate *> Sates, EPHEMERIS **eph, bool first_flag, double f, MatrixXd *B_Pos, MatrixXd *l_Pos, MatrixXd *P_Pos, string *sate_used);

// ��ٶȽ��㹹�����
unsigned int setup_Vel(GPSTIME *OBS_TIME, MatrixXd Pos, vector<Satellate *> Sates, EPHEMERIS **eph, MatrixXd *B_Vel, MatrixXd *l_Vel, MatrixXd *P_Vel);

unsigned int LS_SPV(DATA_SET* data, Configure cfg);

unsigned int setup_LS(DATA_SET* data, Configure cfg, int sys);

unsigned int setup_KF(DATA_SET* data, Configure cfg, int sys);

double get_measure(Satellate* sate, Configure cfg, EPHEMERIS* eph);

// SPP���㶨λKF
unsigned int KF_SPV(DATA_SET* data, double dt_e, Configure cfg);