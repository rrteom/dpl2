#ifndef VMESH_HPP
#define VMESH_HPP

#include "structures.hpp"

double initDistrFunction(double v_x, double v_y, double T_1, double T_2);
double expTemp1(double v_x, double v_y);
double expTemp2(double v_x, double v_y, double temp_1, double temp_2);

class VMesh {
    OneDArr v_x_mesh, v_y_mesh, init_v_distr;
    TwoDArr index_to_p;
    double v_cut;
    int max_p;
    int n_v_x, n_v_y;
public:
    OneDArr exp_t1, exp_t2, abs_v_x, sign_v_x, abs_v_y, sign_v_y;
    VMesh(double v_cut, int n_v_x, int n_v_y, double temp_1, double temp_2);
    int indexToP(Index2d i); 
    int safeIndexToP(Index2d i);
    int maxP();
    OneDArr getInitVDistr();
    double getVCut();
    double gridLockX(double v_x);
    double gridLockY(double v_y);

    Index2d getClosestVelocityIndex(Vector2d v);
    Index2d getFloorVelocityIndex(Vector2d v);
    double getNodeEnergy(Index2d node, Vector2d v_c);
};

#endif