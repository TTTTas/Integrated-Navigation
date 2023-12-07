#pragma once
#include <Eigen/Dense>
#include <iostream>
using namespace std;
using namespace Eigen;

class KalmanFilter
{
public:
    KalmanFilter(MatrixXd initial_state, MatrixXd initial_covariance);

    void predict();

    void update(MatrixXd z);

    void set_A(MatrixXd A);

    void set_H(MatrixXd H);

    void set_Q(MatrixXd Q);

    void set_R(MatrixXd R);

    void setState(MatrixXd state);

    MatrixXd getState() const;

private:
    MatrixXd A_; // ״̬ת�ƾ���
    MatrixXd Q_; // ��������Э����
    MatrixXd H_; // �۲����
    MatrixXd R_; // ��������Э����

    MatrixXd x_hat_; // ״̬����
    MatrixXd P_;     // �������Э����

    MatrixXd x_hat_minus_; // Ԥ��״̬����
    MatrixXd P_minus_;     // Ԥ��������Э����

    MatrixXd K_; // ����������
};

MatrixXd getA(double T);

MatrixXd getH(MatrixXd B_pos, MatrixXd B_Vel);

MatrixXd getQ(double T);

MatrixXd getR(double ROW_P,double ROW_V);

MatrixXd getz(MatrixXd l_P,MatrixXd l_V);