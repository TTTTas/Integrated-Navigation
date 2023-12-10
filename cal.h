#pragma once
#include "Common_value.h"
#include "read.h"
#include "transform.h"
#include "Configure.h"

using namespace std;

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

// ���ǽ���
unsigned int Cal_1(DATA_SET *data, OBS_DATA *obs, EPOCH *eph, bool first_flag);

// GPS��BDS˫�ǽ���
unsigned int Cal_2(DATA_SET *data, OBS_DATA *obs, EPOCH *gpseph, EPOCH *bdseph, bool first_flag);

// SPP���㶨λ
unsigned int Cal_SPP(DATA_SET *data, OBS_DATA *obs, EPOCH *gpseph, EPOCH *bdseph, double dt_e, bool first_flag);



/*������λ�ý���*/
// �λ�ý��㹹�����(˫Ƶ)
unsigned int setup_Pos(GPSTIME *OBS_TIME, MatrixXd Pos, vector<Satellate *> Sates, EPHEMERIS **eph, bool first_flag, double f1, double f2, MatrixXd *B_Pos, MatrixXd *l_Pos, MatrixXd *P_Pos, string *sate_used);

unsigned int setup_Pos(DATA_SET* data, Configure cfg);

// �λ�ý��㹹�����(��Ƶ)
unsigned int setup_Pos(GPSTIME *OBS_TIME, MatrixXd Pos, vector<Satellate *> Sates, EPHEMERIS **eph, bool first_flag, double f, MatrixXd *B_Pos, MatrixXd *l_Pos, MatrixXd *P_Pos, string *sate_used);

// ��ٶȽ��㹹�����
unsigned int setup_Vel(GPSTIME *OBS_TIME, MatrixXd Pos, vector<Satellate *> Sates, EPHEMERIS **eph, MatrixXd *B_Vel, MatrixXd *l_Vel, MatrixXd *P_Vel);

// ���ǽ���
unsigned int Cal_1(DATA_SET *data, bool first_flag);

// GPS��BDS˫�ǽ���
unsigned int Cal_2(DATA_SET *data, bool first_flag);

// SPP���㶨λ
unsigned int Cal_SPP(DATA_SET *data, double dt_e, bool first_flag);