//
// Created by DELL on 2021/4/23.
//
#include "INS.h"

INS::ErrTrans INS::get_state_transition(Vector3d Fn) const {
    double fe = Fn(1), fn = Fn(2), fu = Fn(3);
    ErrTrans state_transition_matrix = MatrixXd::Zero(18, 18);
    double ve2 = v_E * v_E, vn2 = v_N * v_N;
    double rn2 = Rn * Rn, rm2 = Rm * Rm;
    state_transition_matrix(0, 1) = wie * sin(longitude) + v_E * tan(longitude) / Rn;
    state_transition_matrix(0, 2) = -(wie * cos(longitude) + v_E / Rn);
    state_transition_matrix(0, 4) = -1 / Rm;
    state_transition_matrix(0, 8) = v_N / rm2;

    state_transition_matrix(1, 0) = -(wie * sin(longitude) + v_E * tan(longitude) / Rn);
    state_transition_matrix(1, 2) = -v_N / Rm;
    state_transition_matrix(1, 3) = 1 / Rn;
    state_transition_matrix(1, 6) = -wie * sin(longitude);
    state_transition_matrix(1, 8) = -v_E / rn2;

    state_transition_matrix(2, 0) = wie * cos(longitude) + v_E / Rn;
    state_transition_matrix(2, 1) = v_N / Rm;
    state_transition_matrix(2, 3) = tan(longitude) / Rn;
    state_transition_matrix(2, 6) = wie * cos(longitude) + v_E * (1 + tan(longitude) * tan(longitude)) / Rn;
    state_transition_matrix(2, 8) = -v_E * tan(longitude) / rn2;

    state_transition_matrix(3, 1) = -fu;
    state_transition_matrix(3, 2) = fn;
    state_transition_matrix(3, 3) = v_N * tan(longitude) / Rm - v_U / Rm;
    state_transition_matrix(3, 4) = 2 * wie * sin(longitude) + v_E * tan(longitude) / Rn;
    state_transition_matrix(3, 5) = -(2 * wie * cos(longitude) + v_E / Rn);
    state_transition_matrix(3, 6) = 2 * wie * cos(longitude) * v_N + v_E * v_N * (1 + tan(longitude) * tan(longitude)) / Rn + 2 * wie * sin(longitude) * v_U;
    state_transition_matrix(3, 8) = (v_E * v_U - v_E * v_N * tan(longitude)) / rn2;

    state_transition_matrix(4, 0) = fu;
    state_transition_matrix(4, 2) = -fe;
    state_transition_matrix(4, 3) = -2 * (wie * sin(longitude) + v_E * tan(longitude) / Rn);
    state_transition_matrix(4, 4) = -v_U / Rm;
    state_transition_matrix(4, 5) = -v_N / Rm;
    state_transition_matrix(4, 6) = -(2 * wie * cos(longitude) + v_E * ((1 + tan(longitude) * tan(longitude))) / Rn) * v_E;
    state_transition_matrix(4, 8) = (ve2 * tan(longitude) + v_N * v_U) / rn2;

    state_transition_matrix(5, 0) = -fn;
    state_transition_matrix(5, 1) = fe;
    state_transition_matrix(5, 3) = 2 * (wie * cos(longitude) + v_E / Rn);
    state_transition_matrix(5, 4) = 2 * v_N / Rm;
    state_transition_matrix(5, 6) = -2 * v_E * wie * sin(longitude);
    state_transition_matrix(5, 8) = -(vn2 + ve2) / rn2;

    state_transition_matrix(6, 4) = 1 / Rm;
    state_transition_matrix(7, 3) = 1 / (Rn * cos(longitude));
    state_transition_matrix(7, 6) = v_E * tan(longitude) / (Rn * cos(longitude));
    state_transition_matrix(8, 5) = 1;

    state_transition_matrix.block(0, 9, 3, 3) = Cnb;
    state_transition_matrix.block(0, 12, 3, 3) = Cnb;
    state_transition_matrix.block(3, 15, 3, 3) = Cnb;
    state_transition_matrix.block(12, 12, 3, 3) = -1 / Tg * Matrix3d::Zero();
    state_transition_matrix.block(15, 15, 3, 3) = -1 / Ta * Matrix3d::Zero();
    
    return state_transition_matrix;
}
