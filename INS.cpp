//
// Created by DELL on 2021/4/21.
//

#include "INS.h"
#include <eigen3/Eigen/Geometry>

/**
 * 规定，
 * 对于坐标系下标使用小写字母，
 * 对于方向下标使用大写字母，
 * 如Vn表示速度在导航系中的投影
 * vN表示北向速度*/

/**
 * 规定俯仰角theta介于[-pi/2,pi/2]，可以证明，对于任意的姿态，都有一个theta \in [-pi/2,pi/2]和对应的gamma,psi满足该姿态
 * 也就是说从(theta, psi, gamma)到姿态矩阵的映射是一个多对一的映射，在通常情况下是二对一，在theta=+-pi/2时是无穷对一
 * 那么我们可以直接通过反正弦函数来确定theta的值
 * 但是这么做会导致航向角psi的不连续，比如说规定psi=0,gamma=0，theta从0变到pi，theta=pi*t，t为时间参数
 * 那么实际上gamma会一直等于0，在前半段psi=0，当theta越过pi/2时，psi会跃变为pi，theta的实际取值为theta=pi/2-|pi/2-t|
 *
 * 确定俯仰角之后，如果T(2,2)<0或T(3,3)<0，那么便可以确定航向角/滚转角在后半球。
 * 而通过反正切函数求出来的航向角/滚转角一定在前半球，加减pi即可。
 * */

inline Matrix3d atti2matrix(double yaw, double row, double pitch) {
    return (Matrix3d)AngleAxisd(row, Eigen::Vector3d::UnitZ()) *
        AngleAxisd(pitch, Eigen::Vector3d::UnitY()) *
        AngleAxisd(yaw, Eigen::Vector3d::UnitX());
}

INS::INS(const Vector3d& posi_in, const Vector3d& atti_in, const Vector3d& v_in, const Err& err_in,
    const Cov& cov_in) {
    latitude = posi_in(0);
    longitude = posi_in(1);
    height = posi_in(2);
    sin_lati = sin(latitude);
    sin_lati2 = sin_lati * sin_lati;
    Rm = Re * (1 - 2 * f + 3 * f * sin_lati2) + height;
    Rm = Re * (1 + f * sin_lati2) + height;

    yaw = atti_in(0);
    pitch = atti_in(1);
    row = atti_in(2);

    v_E = v_in(0);
    v_N = v_in(0);
    v_U = v_in(0);

    err = err_in;
    cov = cov_in;

    Cnb = atti2matrix(yaw, row, pitch);

    kf_ptr=unique_ptr<KF>(new KF(sys_noise_cov,measure_noise_cov,ts));

    Vector3d measure_noise_std_cov_vec(1e-5, 1e-5, 10.3986);
    measure_noise_cov = (measure_noise_std_cov_vec.cwiseProduct(measure_noise_std_cov_vec)).asDiagonal();
    Matrix<double, 9, 1> sys_noise_std_cov_vec;
    sys_noise_std_cov_vec << 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.9780, 0.9780, 0.9780;
    sys_noise_cov
    = 1e-6 * (sys_noise_std_cov_vec.cwiseProduct(sys_noise_std_cov_vec)).asDiagonal();
}

tuple<INS::Posi, INS::Atti, INS::Vecl>
INS::solve(const Vector3d& force_b, const Vector3d& omega_b_ib, const Vector3d& gps, bool is_use_gps) {
    Vector3d force_n = Cnb * force_b;

    err_trans = get_state_transition(force_n);
    sys_noise_transfer = SysNoiseTrans::Zero();
    sys_noise_transfer.block(0, 0, 3, 3) = Cnb;
    sys_noise_transfer.block(12, 3, 3, 3) = Matrix3d::Identity();
    sys_noise_transfer.block(15, 6, 3, 3) = Matrix3d::Identity();
    state2measure = State2Measure::Zero();
    state2measure.block(0, 6, 3, 3) = Matrix3d::Identity();

    if (is_use_gps) {
        kalman_adjust(gps);
    } else {
        kf_ptr->update_without_measure(err_trans, sys_noise_transfer, cov, err);
    }

    inertial_update(force_n, omega_b_ib);
    Vector3d posi(latitude,longitude,height);
    Vector3d atti(pitch,row,yaw);
    Vector3d vec(v_E,v_N,v_U);
    return make_tuple(posi,atti,vec);
}

void INS::kalman_adjust(const Vector3d& gps) {
    Vector3d err_measure = gps;
    err_measure(0) = gps(0) - latitude;
    err_measure(1) = gps(1) - longitude;
    err_measure(2) = gps(2) - height;

    kf_ptr->update_with_measure(err_trans, sys_noise_transfer, state2measure, err_measure, cov, err);

    pitch += err(0);
    row += err(1);
    yaw += err(2);
    Cnb = atti2matrix(yaw, row, pitch);

    v_E += err(3);
    v_N += err(4);
    v_U += err(5);

    latitude += err(6);
    longitude += err(7);
    height += err(8);

    update();

    err.block(0, 0, 9, 1) = Matrix<double,9,1>::Zero();
}

inline Quaterniond q_add(Quaterniond q1, Quaterniond q2) {
    return Quaterniond(q1.w() + q2.w(), q1.x() + q2.x(), q1.y() + q2.y(), q1.z() + q2.z());
}

inline Quaterniond q_scalar_multiplication(Quaterniond q, double scalar) {
    return Quaterniond(q.w() * scalar, q.x() * scalar, q.y() * scalar, q.z() * scalar);
}

void INS::inertial_update(const Vector3d& force_b, const Vector3d& omega_b_ib) {
    Vector3d force_n = Cnb * force_b;

    Vector3d omega_n_ie(0,wie*cos_lati,wie*sin_lati);
    Vector3d omega_n_en(-v_N / Rm, v_E / Rn, v_E * tan(latitude) / Rn);

    Vector3d omega_b_nb = omega_b_ib - Cnb.transpose() * (omega_n_ie + omega_n_en);
    Quaterniond q(Cnb);
    Quaterniond delta_q(0, ts * omega_b_nb(0), ts * omega_b_nb(1), ts * omega_b_nb(2));
    q = q_add(q, q_scalar_multiplication(q * delta_q, 0.5));
    q.normalize();
    Cnb = q.toRotationMatrix();
    Vector3d atti = Cnb.eulerAngles(0, 1, 2);
    yaw = atti(0);
    row = atti(1);
    pitch = atti(2);

    double g = g0 + 0.051799 * sin_lati2 - 0.94114e-6 * height;
    Vector3d gn (0, 0, g);
    Vector3d v;
    v << v_E, v_N, v_U;
    v = v + ts * (force_n + gn - (2 * omega_n_ie + omega_n_en).cross(v));
    v_E = v(0);
    v_N = v(1);
    v_U = v(2);

    latitude += ts * v_N / Rm;
    longitude += ts * v_E / Rn / cos(latitude);
    height += ts * v_U;

    update();
}

void INS::update() {
    cos_lati = cos(latitude);
    sin_lati = sin(latitude);
    sin_lati2 = sin_lati * sin_lati;
    Rm = Re * (1 - 2 * f + 3 * f * sin_lati2) + height;
    Rm = Re * (1 + f * sin_lati2) + height;
}