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
    /*文件输入数据存储*/
    vector<OBS_DATA *> RANGE;
    EPOCH GPS_EPH[GPS_SAT_QUAN];
    EPOCH BDS_EPH[BDS_SAT_QUAN];
    /*存储变量初始化*/
    initial();
    /*结果存储变量*/
    DATA_SET *result = new DATA_SET();
    /*初始判断变量*/
    double dt_epoch = 1; // 文件流历元间时间差
    double temp_t = 0;
    bool first = true; // 第一次解算标志
    bool KF_first = true;
    FILE *DATA_Fobs;   // log文件指针
    FILE *Pos_Fobs;    // pos文件指针
    /*网口输入数据相关变量*/
    int lenR, lenD;
    unsigned char curbuff[MAXRAWLEN];
    lenD = 0;
    unsigned char decBuff[2 * MAXRAWLEN];
    SOCKET NetGps;
    Configure CfgInfo;
    int ReadFlag;

    int choice = 0;
    printf("请选择输入方式\n1. 文件\t2. 网口\n");
    cin >> choice;
    /*获取文件生成时间*/
    time_t nowtime;
    time(&nowtime); // 获取1970年1月1日0点0分0秒到现在经过的秒数
    tm p;
    localtime_s(&p, &nowtime); // 将秒数转换为本地时间,年从1900算起,需要+1900,月为0-11,所以要+1
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
    string path = "D:/GitHub/SPP_Design/报告/双频双系统.pos";

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
            if ((lenR = recv(NetGps, (char *)curbuff, MAXRAWLEN, 0)) > 0) // 读取数据
            {
                printf("%5d\n", lenR);
                fwrite(curbuff, sizeof(unsigned char), lenR, DATA_Fobs); // 记录二进制数据流到文件中

                if ((lenD + lenR) > 2 * MAXRAWLEN)
                    lenD = 0;
                memcpy(decBuff + lenD, curbuff, lenR); // 缓存拼接
                lenD += lenR;

                ReadFlag = decodestream(result, decBuff, lenD); // 解码
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
                     result->LS_print();              // 输出至控制台
                     result->LS_Filewrite(Pos_Fobs); // 输出至文件

 
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