//
// Created by DELL on 2021/4/20.
//

#ifndef INS_PARAMETERS_H
#define INS_PARAMETERS_H

#include <eigen3/Eigen/Core>
#include <vector>
#include <string>

using namespace Eigen;
using namespace std;
const double Re = 6378160;//%地球半径(长半轴)
const double f = 1 / 298.3;//地球扁率
const double e = sqrt(2 * f - f * f);//偏心率
const double e2 = e * e;
const double Rp = (1 - f) * Re;//短半轴
const double ep = sqrt(Re * Re - Rp * Rp) / Rp;//第二偏心率,  此处有改动，加号改为减号
const double ep2 = ep * ep;
const double wie = 7.2921151467e-5;//地球自转角速率
const double g0 = 9.7803267714;//重力加速度
const double ts = 0.1;//sample time

const double pi = 3.1415926;
const double arg_deg = pi / 180;
const double arg_min = arg_deg / 60;
const double arg_sec = arg_min / 60;
const double deg_per_hour = arg_deg / 3600;

const string gps_data = "gps.txt";
const string force_data = "force.txt";
const string omega_data = "omega.txt";


constexpr int loop_num = 10;

#endif //INS_PARAMETERS_H
