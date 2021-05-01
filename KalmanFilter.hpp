//
// Created by DELL on 2021/4/20.
//

#ifndef INS_KALMANFILTER_H
#define INS_KALMANFILTER_H

#include <eigen3/Eigen/Core>

using namespace Eigen;

template<unsigned state_dim, unsigned sys_noise_dim, unsigned measure_dim>
class KalmanFilter {
private:

public:
    typedef Matrix<double, state_dim, state_dim> StateTrans;
    typedef Matrix<double, state_dim, sys_noise_dim> SysNoiseTrans;
    typedef Matrix<double, measure_dim, state_dim> State2Measure;

    typedef Matrix<double, measure_dim, 1> Measure;
    typedef Matrix<double, state_dim, state_dim> Cov;
    typedef Matrix<double, state_dim, 1> State;

    typedef Matrix<double, sys_noise_dim, sys_noise_dim> SysNoiseCov;
    typedef Matrix<double, measure_dim, measure_dim> MeasureNoiseCov;
    SysNoiseCov sys_noise_cov;
    MeasureNoiseCov measure_noise_cov;

    double ts;

    KalmanFilter(const SysNoiseCov& sys_noise_cov_, const MeasureNoiseCov& measure_noise_cov_, const double& ts_) {
        sys_noise_cov = sys_noise_cov_;
        measure_noise_cov = measure_noise_cov_;
        ts = ts_;
    }

    void update_with_measure(
        const StateTrans& state_trans,
        const SysNoiseTrans& sys_noise_trans,
        const State2Measure& state2measure,
        const Measure& measure,
        Cov& cov,
        State& x) {

        auto E = MatrixXd::Ones(state_dim, state_dim);
        StateTrans state_trans_discrete = E + ts * state_trans;
        SysNoiseTrans sys_noise_trans_discrete = (E + 0.5 * ts * state_trans) * sys_noise_trans * ts;
        State x_predict = state_trans_discrete * x;
        Cov cov_predict = state_trans_discrete * cov * state_trans_discrete.transpose()
            + sys_noise_trans_discrete * sys_noise_cov * sys_noise_trans_discrete.transpose();
        MatrixXd kalman_gain = cov_predict * state2measure.transpose() *
            (state2measure * cov_predict * state2measure.transpose() +
                measure_noise_cov).inverse();

        MatrixXd tmp = E - kalman_gain * state2measure;
        cov = tmp * cov_predict * tmp + kalman_gain * measure_noise_cov * kalman_gain.transpose();
        x = x + kalman_gain * (measure - state2measure * x_predict);
    }

    void update_without_measure(
        const StateTrans& state_trans,
        const SysNoiseTrans& sys_noise_trans,
        Cov& cov,
        State& x) {
        auto E = MatrixXd::Ones(state_dim, state_dim);
        StateTrans state_trans_discrete = E + ts * state_trans;
        SysNoiseTrans sys_noise_trans_discrete = (E + 0.5 * ts * state_trans) * sys_noise_trans * ts;
        x = state_trans_discrete * x;
        cov = state_trans_discrete * cov * state_trans_discrete.transpose()
            + sys_noise_trans_discrete * sys_noise_cov * sys_noise_trans_discrete.transpose();
    }

};

#endif //INS_KALMANFILTER_H
