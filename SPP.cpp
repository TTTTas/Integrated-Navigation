#include <iostream>
#include <fstream>
#include <vector>
#include <Eigen/Dense>
#include <iomanip>
#include <io.h>
#include <direct.h>

#include "transform.h"
#include "read.h"
#include "cal.h"
#include "sockets.h"
#include "KF.h"
#include "data.h"

using namespace std;
using namespace Eigen;

#define _CRT_SECURE_NO_WARNINGS

int main()
{
    /*�ļ��������ݴ洢*/
    vector<OBS_DATA *> RANGE;
    EPOCH GPS_EPH[GPS_SAT_QUAN];
    EPOCH BDS_EPH[BDS_SAT_QUAN];
    /*�洢������ʼ��*/
    initial();
    /*����洢����*/
    DATA_SET *result = new DATA_SET();
    /*��ʼ�жϱ���*/
    double dt_epoch = 1; // �ļ�����Ԫ��ʱ���
    double temp_t = 0;
    bool first = true; // ��һ�ν����־
    bool KF_first = true;
    FILE *DATA_Fobs;   // log�ļ�ָ��
    FILE *Pos_Fobs;    // pos�ļ�ָ��
    /*��������������ر���*/
    int lenR, lenD;
    unsigned char curbuff[MAXRAWLEN];
    lenD = 0;
    unsigned char decBuff[2 * MAXRAWLEN];
    SOCKET NetGps;
    Configure CfgInfo;
    int ReadFlag;

    int choice = 0;
    printf("��ѡ�����뷽ʽ\n1. �ļ�\t2. ����\n");
    cin >> choice;
    /*��ȡ�ļ�����ʱ��*/
    time_t nowtime;
    time(&nowtime); // ��ȡ1970��1��1��0��0��0�뵽���ھ���������
    tm p;
    localtime_s(&p, &nowtime); // ������ת��Ϊ����ʱ��,���1900����,��Ҫ+1900,��Ϊ0-11,����Ҫ+1
    string filetime = to_string(p.tm_year + 1900) + "_" + to_string(p.tm_mon + 1) + "_" + to_string(p.tm_mday) + "_" + to_string(p.tm_hour) + "_" + to_string(p.tm_min) + "_" + to_string(p.tm_sec);
    string logpath = "C:\\Users\\Surface\\Desktop\\data\\logs\\";
    createDirectory(logpath);
    logpath += filetime + string(".log");
    string pospath = "C:\\Users\\Surface\\Desktop\\data\\Pos\\";
    createDirectory(pospath);
    pospath += filetime + string(".pos");
    CfgInfo.ObsDatFile = logpath.c_str();
    CfgInfo.ResDatFile = pospath.c_str();

    FILE *tempfile;
    string path = "D:/GitHub/SPP_Design/����/˫Ƶ˫ϵͳ.pos";

    FILE* KF_file;

    switch (choice)
    {
    case 1:
        readfile("NovatelOEM20211114-01.log", RANGE, GPS_EPH, BDS_EPH);
        if ((tempfile = fopen(path.c_str(), "w")) == NULL)
        {
            printf("The pos file % s was not opened\n", path.c_str());
            exit(0);
        }
        for (int i = 0; i < RANGE.size(); i++)
        {
            if (i > 0)
                dt_epoch = RANGE[i]->OBS_TIME->SecOfWeek - RANGE[i - 1]->OBS_TIME->SecOfWeek;
            if (Cal_SPP(result, RANGE[i], GPS_EPH, BDS_EPH, dt_epoch, first))
                first = false;
            result->LS_print();
            result->LS_Filewrite(tempfile);
        }
        fclose(tempfile);
        break;
    case 2:
        if (OpenSocket(NetGps, CfgInfo.NetIP, CfgInfo.NetPort) == false)
        {
            printf("The ip %s was not opened\n", CfgInfo.NetIP);
            return 0;
        }
        if ((DATA_Fobs = fopen(CfgInfo.ObsDatFile, "wb")) == NULL)
        {
            printf("The obs file %s was not opened\n", CfgInfo.ObsDatFile);
            exit(0);
        }
        if ((Pos_Fobs = fopen(CfgInfo.ResDatFile, "w")) == NULL)
        {
            printf("The pos file %s was not opened\n", CfgInfo.ResDatFile);
            exit(0);
        }
        if ((KF_file = fopen("KF.pos", "w")) == NULL)
        {
            printf("The pos file %s was not opened\n", "KF.pos");
            exit(0);
        }
        while (1)
        {
            Sleep(980);
            if ((lenR = recv(NetGps, (char *)curbuff, MAXRAWLEN, 0)) > 0) // ��ȡ����
            {
                printf("%5d\n", lenR);
                fwrite(curbuff, sizeof(unsigned char), lenR, DATA_Fobs); // ��¼���������������ļ���

                if ((lenD + lenR) > 2 * MAXRAWLEN)
                    lenD = 0;
                memcpy(decBuff + lenD, curbuff, lenR); // ����ƴ��
                lenD += lenR;

                ReadFlag = decodestream(result, decBuff, lenD); // ����
                if (ReadFlag != 1)
                {
                    printf("Data acquisition and decode failed \n");
                    continue;
                }
                else
                {
                    if (first)
                        dt_epoch = 1;
                    else
                        dt_epoch = result->OBSTIME->SecOfWeek - temp_t;
                    if (dt_epoch == 0)
                        break;
                    temp_t = result->OBSTIME->SecOfWeek;
                    if (Cal_SPP(result, dt_epoch, first))
                    {
                        first = false;
                        cout << "LS" << endl;
                        if (KF_SPP(result, dt_epoch, KF_first))
                            KF_first = false;
                     }
                     result->LS_print();              // ���������̨
                     result->LS_Filewrite(Pos_Fobs); // ������ļ�

 
                    result->KF_Print(KF_file);

                    result->reset();
                }
            }
            else
            {
                printf("NO MESSAGES IN!\n");
            }
        }

        break;
    default:
        break;
    }

    std::system("pause");
    return 0;
}