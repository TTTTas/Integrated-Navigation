#pragma once
#include <Eigen/Dense>
#include <iostream>
#include "Configure.h"
using namespace std;
using namespace Eigen;

class KalmanFilter
{
public:
    KalmanFilter(MatrixXd initial_state, MatrixXd initial_covariance);

    void predict();

    void update();

    void set_A(MatrixXd A);

    void set_H(MatrixXd H);

    void set_Q(MatrixXd Q);

    void set_R(MatrixXd R);

    void set_Z(MatrixXd Z);

    void setState(MatrixXd state);

    void reset();

    MatrixXd getState() const;

    MatrixXd getState_minus();

private:
    MatrixXd A_; // 状态转移矩阵
    MatrixXd Q_; // 过程噪声协方差
    MatrixXd H_; // 观测矩阵
    MatrixXd R_; // 测量噪声协方差

    MatrixXd x_hat_; // 状态估计
    MatrixXd P_;     // 估计误差协方差

    MatrixXd x_hat_minus_; // 预测状态估计
    MatrixXd P_minus_;     // 预测估计误差协方差

    MatrixXd K_; // 卡尔曼增益
    MatrixXd z_ ; //观测值
};

MatrixXd getA(double T, Configure cfg);

MatrixXd getH(MatrixXd B);

MatrixXd getQ(double T, Configure cfg);

MatrixXd getR(double ROW);

MatrixXd getz(MatrixXd l_P,MatrixXd l_V);