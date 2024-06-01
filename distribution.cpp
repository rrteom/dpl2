#include "distribution.hpp"

double sgn (double x) {
    return (x > 0) ? 1 : ((x < 0) ? -1 : 0);
}

double mc_limiter_func (double f_i_minus, double f_i, double f_i_plus) {
    if ((f_i_plus - f_i) * (f_i - f_i_minus) <= 0)
        return 0;
    double pre_min = std::min(std::abs(f_i_plus - f_i_minus) / 4,std::abs(f_i_plus - f_i));
    return std::min(pre_min, std::abs(f_i - f_i_minus)) * sgn(f_i_plus - f_i);
}

Distribution::Distribution(double L_x, double L_y, int n_h_x, int n_h_y,
                 double v_cut, int n_v_x, int n_v_y,
                 double tau, double temp_1, double temp_2,
                 int ix_chip_start, int iy_chip_start, int ix_chip_end, int iy_chip_end) : 
            L_x(L_x), L_y(L_y), n_h_x(n_h_x), n_h_y(n_h_y), h_x(L_x / n_h_x), h_y(L_y / n_h_y),
            v_cut(v_cut), n_v_x(n_v_x), n_v_y(n_v_y),
            tau(tau), temp_1(temp_1), temp_2(temp_2),
            ix_chip_start(ix_chip_start), iy_chip_start(iy_chip_start), ix_chip_end(ix_chip_end), iy_chip_end(iy_chip_end),
            x_mesh(0, L_x, n_h_x), y_mesh(0, L_y, n_h_y), v_mesh(v_cut, n_v_x, n_v_y, temp_1, temp_2),
            chip_mask(n_h_x, n_h_y) 
{
    // 1 if no chip, 0 for tiles filled with chip
    chip_mask.fill(1);
    for (int ix = ix_chip_start; ix <= ix_chip_end; ix++) {
        for (int iy = iy_chip_start; iy <= iy_chip_end; iy++) {
            chip_mask.at(ix, iy) = 0;
        }
    }
    max_p = v_mesh.maxP();
    OneDArr init_v_distr = v_mesh.getInitVDistr();
    distr.resize(n_h_x, n_h_y, max_p);
    for (int ix = 1; ix <= n_h_x; ix ++) {
        for (int iy = 1; iy <= n_h_y; iy++) {
            if (chip_mask.at(ix, iy)) {
                for (int p = 1; p <= max_p; p++) {
                    distr.at(ix, iy, p) = init_v_distr.at(p);
                }
            }
        }
    }
}

TwoDArr Distribution::getChipMask() {return chip_mask;}

void Distribution::stepX(double tau_x) {
    for (int iy = 1; iy <= n_h_y; iy++) {
        if ((iy >= iy_chip_start) and (iy <= iy_chip_end)) {
            TwoDArr f_one_d_left(ix_chip_start - 1, max_p), f_one_d_right(n_h_x - ix_chip_end, max_p);
            for (int ix = 1; ix < ix_chip_start; ix++) {
                for (int p = 1; p <= max_p; p++) {
                    f_one_d_left.at(ix, p) = distr.at(ix, iy, p);
                }
            }
            for (int ix = ix_chip_end + 1; ix <= n_h_x; ix++) {
                for (int p = 1; p <= max_p; p++) {
                    f_one_d_right.at(ix - ix_chip_end, p) = distr.at(ix, iy, p);
                }
            }
            solveOneDTask(f_one_d_left, temp_2, temp_1, tau_x, h_x, v_mesh.abs_v_x, v_mesh.sign_v_x, v_mesh.exp_t2, v_mesh.exp_t1);    
            solveOneDTask(f_one_d_right, temp_1, temp_2, tau_x, h_x, v_mesh.abs_v_x, v_mesh.sign_v_x, v_mesh.exp_t1, v_mesh.exp_t2);
            //update distr
            for (int ix = 1; ix < ix_chip_start; ix++) {
                for (int p = 1; p <= max_p; p++) {
                    distr.at(ix, iy, p) = f_one_d_left.at(ix, p);
                }
            }
            for (int ix = ix_chip_end + 1; ix <= n_h_x; ix++) {
                for (int p = 1; p <= max_p; p++) {
                    distr.at(ix, iy, p) = f_one_d_right.at(ix - ix_chip_end, p);
                }
            }
        }
        else {
            TwoDArr f_one_d(n_h_x, max_p);
            for (int ix = 1; ix <= n_h_x; ix++) {
                for (int p = 1; p <= max_p; p++) {
                    f_one_d.at(ix, p) = distr.at(ix, iy, p);
                }
            }
            solveOneDTask(f_one_d, temp_2, temp_2, tau_x, h_x, v_mesh.abs_v_x, v_mesh.sign_v_x, v_mesh.exp_t2, v_mesh.exp_t2);
            // update distr
            for (int ix = 1; ix <= n_h_x; ix++) {
                for (int p = 1; p <= max_p; p++) {
                    distr.at(ix, iy, p) = f_one_d.at(ix, p);
                }
            }
        }
        
    }
}
void Distribution::stepY(double tau_y) {
    for (int ix = 1; ix <= n_h_x; ix++) {
        if ((ix >= ix_chip_start) and (ix <= ix_chip_end)) {
            TwoDArr f_one_d_upper(iy_chip_start - 1, max_p), f_one_d_lower(n_h_y - iy_chip_end, max_p);
            for (int iy = 1; iy < iy_chip_start; iy++) {
                for (int p = 1; p <= max_p; p++) {
                    f_one_d_upper.at(iy, p) = distr.at(ix, iy, p);
                }
            }
            for (int iy = iy_chip_end + 1; iy <= n_h_y; iy++) {
                for (int p = 1; p <= max_p; p++) {
                    f_one_d_lower.at(iy - iy_chip_end, p) = distr.at(ix, iy, p);
                }
            }
            solveOneDTask(f_one_d_upper, temp_2, temp_1, tau_y, h_y, v_mesh.abs_v_y, v_mesh.sign_v_y, v_mesh.exp_t2, v_mesh.exp_t1);    
            solveOneDTask(f_one_d_lower, temp_1, temp_2, tau_y, h_y, v_mesh.abs_v_y, v_mesh.sign_v_y, v_mesh.exp_t1, v_mesh.exp_t2);
            //uodate distr
            for (int iy = 1; iy < iy_chip_start; iy++) {
                for (int p = 1; p <= max_p; p++) {
                    distr.at(ix, iy, p) = f_one_d_upper.at(iy, p);
                }
            }
            for (int iy = iy_chip_end + 1; iy <= n_h_y; iy++) {
                for (int p = 1; p <= max_p; p++) {
                    distr.at(ix, iy, p) = f_one_d_lower.at(iy - iy_chip_end, p);
                }
            }
        }
        else {
            TwoDArr f_one_d(n_h_y, max_p);
            for (int iy = 1; iy <= n_h_y; iy++) 
                for (int p = 1; p <= max_p; p++) 
                    f_one_d.at(iy, p) = distr.at(ix, iy, p);
            solveOneDTask(f_one_d, temp_2, temp_2, tau_y, h_y, v_mesh.abs_v_y, v_mesh.sign_v_y, v_mesh.exp_t2, v_mesh.exp_t2);
            //update distr
            for (int iy = 1; iy <= n_h_y; iy++) 
                for (int p = 1; p <= max_p; p++) 
                    distr.at(ix, iy, p) = f_one_d.at(iy, p);
        }
    }
}

void Distribution::step1DY(double tau_y) {
    TwoDArr f_one_d(n_h_y, max_p);
    int ix = 1;
    for (int iy = 1; iy <= n_h_y; iy++) 
        for (int p = 1; p <= max_p; p++) 
            f_one_d.at(iy, p) = distr.at(ix, iy, p);
    solveOneDTask(f_one_d, temp_2, temp_1, tau_y, h_y, v_mesh.abs_v_y, v_mesh.sign_v_y, v_mesh.exp_t2, v_mesh.exp_t1);
    //update distr
    for (int iy = 1; iy <= n_h_y; iy++) 
        for (int p = 1; p <= max_p; p++) 
            distr.at(ix, iy, p) = f_one_d.at(iy, p);
}


void Distribution::solveOneDTask(TwoDArr& f, double t_1, double t_2, double tau, double h_step,
                                 OneDArr& v_proj_abs, OneDArr& v_proj_sign, OneDArr& v_exp_t_1, OneDArr& v_exp_t_2) 
{
    int n_h = f.dimx_, p_max = v_mesh.maxP();
    TwoDArr f_half(n_h + 1, p_max); // plus one in h space
    OneDArr f_zero(p_max), f_plus(p_max);

    // far from edges
    for (int p = 1; p <= p_max; p++) {
        double current_gamma = v_proj_abs.at(p) * tau / h_step;
        if (v_proj_sign.at(p) < 0) {
            for (int ih = 2; ih <= n_h - 1; ih++) {
                f_half.at(ih, p) = f.at(ih, p) - (1 - current_gamma) * mc_limiter_func(f.at(ih - 1, p), f.at(ih, p), f.at(ih + 1, p));
            }
        }
        else {
            for (int ih = 3; ih <= n_h; ih++) {
                f_half.at(ih, p) = f.at(ih - 1, p) + (1 - current_gamma) * mc_limiter_func(f.at(ih - 2, p), f.at(ih - 1, p), f.at(ih, p));
            }
        }
    }

    // "coming to" edge
    for (int p = 1; p <= p_max; p++) {
        double current_gamma = v_proj_abs.at(p) * tau / h_step;
        if (v_proj_sign.at(p) < 0) {
            f_zero.at(p) = std::max(0.0d, 2 * f.at(1, p) - f.at(2, p));
            f_half.at(1, p) = f.at(1, p) - (1 - current_gamma) * mc_limiter_func(f_zero.at(p), f.at(1, p), f.at(2, p));
        }
        else {
            f_plus.at(p) = std::max(0.0d, 2 * f.at(n_h, p) - f.at(n_h - 1, p));
            f_half.at(n_h + 1, p) = f.at(n_h, p) + (1 - current_gamma) * mc_limiter_func(f.at(n_h - 1, p), f.at(n_h, p), f_plus.at(p));
        }
    }

    // diffuse reflexion norms
    double wall_1_out = 0, wall_2_out = 0, wall_1_in = 0, wall_2_in = 0, border_1_in = 0, border_2_in = 0;
    for (int p = 1; p <= max_p; p++) {
        if (v_proj_sign.at(p) < 0) {
            wall_2_out += v_proj_abs.at(p) * v_exp_t_2.at(p);
            wall_1_in += v_proj_abs.at(p) * f_half.at(1, p);
            border_1_in += v_proj_abs.at(p) * (f_zero.at(p) + f.at(1, p)) / 2;
        }
        else {
            wall_1_out += v_proj_abs.at(p) * v_exp_t_1.at(p);
            wall_2_in += v_proj_abs.at(p) * f_half.at(n_h + 1, p);
            border_2_in += v_proj_abs.at(p) * (f_plus.at(p) + f.at(n_h, p)) / 2;
        }
    }

    //diffuse reflexion
    for (int p = 1; p <= max_p; p++) {
        double current_gamma = v_proj_abs.at(p) * tau / h_step;
        if (v_proj_sign.at(p) < 0) {
            f_half.at(n_h + 1, p) = v_exp_t_2.at(p) * wall_2_in / wall_2_out;
            f_plus.at(p) = std::max(0.0d, 2 * v_exp_t_2.at(p) * border_2_in / wall_2_out - f.at(n_h, p));
            f_half.at(n_h, p) = f.at(n_h, p) - (1 - current_gamma) * mc_limiter_func(f.at(n_h - 1, p), f.at(n_h, p), f_plus.at(p));
        }
        else {
            f_half.at(1, p) = v_exp_t_1.at(p) * wall_1_in / wall_1_out;
            f_zero.at(p) = std::max(0.0d, 2 * v_exp_t_1.at(p) * border_1_in / wall_1_out - f.at(1, p));
            f_half.at(2, p) = f.at(1, p) + (1 - current_gamma) * mc_limiter_func(f_zero.at(p), f.at(1, p), f.at(2, p));
        }
    }

    // updating distr
    for (int p = 1; p <= max_p; p++) {
        double current_gamma = v_proj_abs.at(p) * tau / h_step * v_proj_sign.at(p);
        for (int ih = 1; ih <= n_h; ih++) {
            double new_f = f.at(ih, p) - current_gamma * (f_half.at(ih + 1, p) - f_half.at(ih, p));
            f.at(ih, p) = new_f;
        }
    }
    return;
}


TwoDArr Distribution::getConcentration() {
    TwoDArr n(distr.dimx_, distr.dimy_);
    for (int i = 1; i <= distr.dimx_; i++) {
        for (int j = 1; j <= distr.dimy_; j++) {
            double n_current = 0;
            for (int p = 1; p <= distr.dimp_; p++) {
                n_current += distr.at(i, j, p);
            }
            n.at(i, j) = n_current / v_mesh.c_norm;
        }
    }
    return n;
}

TwoDArr Distribution::getTemperature() {
    TwoDArr n = getConcentration();
    TwoDArr temperature(distr.dimx_, distr.dimy_);
    for (int i = 1; i <= distr.dimx_; i++) {
        for (int j = 1; j <= distr.dimy_; j++) {
            double t_current = 0;
            for (int p = 1; p <= distr.dimp_; p++) {
                t_current += distr.at(i, j, p) * v_mesh.v_squared.at(p);
            }
            if (n.at(i, j) != 0)
                temperature.at(i, j) = t_current  / 2 / n.at(i, j) * temp_1 / v_mesh.c_norm;
        }
    }
    return temperature;
}

double Distribution::getTau0() {
    return std::min(h_x, h_y) / v_cut;
}