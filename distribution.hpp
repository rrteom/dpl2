#ifndef DISTRIBUTION_HPP
#define DISTRIBUTION_HPP

#include "VMesh.hpp"
#include <iostream>

double sgn (double x);

double mc_limiter_func (double f_i_minus, double f_i, double f_i_plus);

class Distribution {
    int n_h_x, n_h_y, n_v_x, n_v_y, ix_chip_start, iy_chip_start, ix_chip_end, iy_chip_end, max_p;
    double tau, L_x, L_y, h_x, h_y, temp_1, temp_2, v_cut;
    VMesh v_mesh;
    OneDArr x_mesh, y_mesh;
    TwoDArr chip_mask;
    ThreeDArr distr;
public:
    Distribution(double L_x, double L_y, int n_h_x, int n_h_y,
                 double v_cut, int n_v_x, int n_v_y,
                 double tau, double temp_1, double temp_2,
                 int ix_chip_start, int iy_chip_start, int ix_chip_end, int iy_chip_end);
    TwoDArr getChipMask();
    void stepX(double tau_x);
    void stepY(double tau_y);
    void solveOneDTask(TwoDArr& distr_1d, double t_1, double t_2, double tau, double h,
                                 OneDArr& v_proj_abs, OneDArr& v_proj_sign, OneDArr& v_exp_t_1, OneDArr& v_exp_t_2);
};

#endif