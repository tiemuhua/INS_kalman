//
// Created by DELL on 2021/4/21.
//

#ifndef INS_INS_H
#define INS_INS_H

#include "KalmanFilter.hpp"
#include <eigen3/Eigen/Core>
#include "parameters.h"
#include <memory>
using namespace std;
using namespace Eigen;

class INS {
public:
    static constexpr int measure_dim = 3;
    static constexpr int state_dim = 18;
    static constexpr int sys_noise_dim = 9;
    typedef Matrix<double, state_dim, 1> Err;
    typedef Matrix<double, state_dim, state_dim> Cov;

private:
    Matrix<double, 9, 9> sys_noise_cov;
    Matrix3d measure_noise_cov;
    
    const double Tg = 3600;//陀螺仪Markov过程相关时间
    const double Ta = 1800;//加速度计Markov过程相关时间

private:
    typedef Matrix<double, INS::state_dim, INS::state_dim> ErrTrans;
    typedef Matrix<double, state_dim, sys_noise_dim> SysNoiseTrans;
    typedef Matrix<double, measure_dim, state_dim> State2Measure;
    ErrTrans err_trans;
    SysNoiseTrans sys_noise_transfer;
    State2Measure state2measure;

    Err err;
    Cov cov;

    double latitude;
    double sin_lati, sin_lati2, cos_lati;
    double longitude;
    double height;

    double yaw;
    double pitch;
    double row;

    double v_E;
    double v_N;
    double v_U;

    double Rm,Rn;

    Matrix3d Cnb;

    void inertial_update(const Vector3d &force_n, const Vector3d& omega_b_ib);

    typedef KalmanFilter<state_dim, sys_noise_dim, measure_dim> KF;
    unique_ptr<KF> kf_ptr;

    ErrTrans get_state_transition(Vector3d Fn) const;

    void update();
    void kalman_adjust(const Vector3d& gps);


public:
    INS(const Vector3d &posi_in, const Vector3d &atti_in, const Vector3d &v_in, const Err &err_in,
             const Cov &cov_in);

    typedef Vector3d Posi;
    typedef Vector3d Atti;
    typedef Vector3d Vecl;

    tuple<INS::Posi, INS::Atti, INS::Vecl> solve(const Vector3d &force_b, const Vector3d &omega_b_nb, const Vector3d &gps, bool is_use_gps);
};


#endif //INS_INS_H
