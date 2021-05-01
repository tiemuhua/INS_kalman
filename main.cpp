#include <iostream>
#include "INS.h"
#include <fstream>
#include "parameters.h"

int main() {
    Vector3d posi, atti, vec;
    posi << 0.5934119457, 1.8849555922, 2.0000000000;
    atti << 0, 0, 0.3491;
    vec << 0, 0, 0;
    INS::Err err = INS::Err::Zero();
    Matrix<double, 18, 1> tmp;
    double var = 0.1 * deg_per_hour;
    tmp << arg_min, arg_min, arg_min, 0.5, 0.5, 0.5, 30 / Re, 30 /
                                                              Re, 30, var, var, var, var, var, var, 1.e-3, 1.e-3, 1.e-3;
    INS::Cov cov = tmp.cwiseProduct(tmp).asDiagonal();
    INS ins(posi, atti, vec, err, cov);

    ifstream gps_fin, force_fin, omega_fin;
    gps_fin.open(gps_data);
    force_fin.open(force_data);
    omega_fin.open(omega_data);
    double tmp1, tmp2, tmp3;
    Vector3d gps, force, omega;
    int cnt = 0;
    while (!gps_fin.eof()) {
        force_fin >> tmp1 >> tmp2 >> tmp3;
        force << tmp1, tmp2, tmp3;
        omega_fin >> tmp1 >> tmp2 >> tmp3;
        omega << tmp1, tmp2, tmp3;
        if (cnt == 9) {
            gps_fin >> tmp1 >> tmp2 >> tmp3;
            gps << tmp1, tmp2, tmp3;
            ins.solve(force, omega, gps, true);
            cnt = 0;
        } else {
            ins.solve(force, omega, gps, false);
            cnt++;
        }
    }
}
