#include "VMesh.hpp"

double initDistrFunction(double v_x, double v_y, double temp_1, double temp_2) {
    // return exp(-0.5 * temp_1 / temp_2 * (v_x * v_x + v_y * v_y)) + exp(-0.5 * temp_1 / temp_2 / 0.95 * (v_x * v_x + v_y * v_y))/ pow(0.95, 1.5);
     return exp(-0.5 * temp_1 / temp_2 * (v_x * v_x + v_y * v_y));
}

double expTemp1(double v_x, double v_y) {
    return exp(-0.5 * (v_x * v_x + v_y * v_y));
}

double expTemp2(double v_x, double v_y, double temp_1, double temp_2) {
    return exp(-0.5 * temp_1 / temp_2 * (v_x * v_x + v_y * v_y));
}


VMesh::VMesh(double v_cut, int n_v_x, int n_v_y, double temp_1, double temp_2) : v_cut(v_cut), v_x_mesh(-v_cut, v_cut, n_v_x), 
                        v_y_mesh(-v_cut, v_cut, n_v_y), n_v_x(n_v_x), n_v_y(n_v_y) {
    double n_0 = 1;
    int p = 0;
    c_norm = 0;
    index_to_p.resize(n_v_x, n_v_y);
    for (int a_x = 1; a_x <= n_v_x; a_x++) {
        for (int a_y = 1; a_y <= n_v_y; a_y++) {
            double v_2 = pow(v_x_mesh.at(a_x), 2) + pow(v_y_mesh.at(a_y), 2);
            if (v_2 <= pow(v_cut, 2)) {
                p++;
                index_to_p.at(a_x, a_y) = p;
                p_to_index_x.push_back(static_cast<double>(a_x));
                p_to_index_y.push_back(static_cast<double>(a_y));
                double init_distr_value = n_0 * initDistrFunction(v_x_mesh.at(a_x), v_y_mesh.at(a_y), temp_1, temp_2);
                c_norm += init_distr_value;
                v_squared.push_back(v_2);
                init_v_distr.push_back(init_distr_value);
                exp_t1.push_back(expTemp1(v_x_mesh.at(a_x), v_y_mesh.at(a_y)));
                exp_t2.push_back(expTemp2(v_x_mesh.at(a_x), v_y_mesh.at(a_y), temp_1, temp_2));
                if (a_x <= n_v_x / 2) {
                    abs_v_x.push_back(v_x_mesh.at(n_v_x - a_x + 1));
                    sign_v_x.push_back(-1);
                }
                else {
                    abs_v_x.push_back(v_x_mesh.at(a_x));
                    sign_v_x.push_back(1);
                }
                if (a_y <= n_v_y / 2) {
                    abs_v_y.push_back(v_y_mesh.at(n_v_y - a_y + 1));
                    sign_v_y.push_back(-1);
                }
                else {
                    abs_v_y.push_back(v_y_mesh.at(a_y));
                    sign_v_y.push_back(1);
                }
            }
        }
    }
    max_p = p;
    v_ph_vol = 4 * v_cut * v_cut / n_v_x / n_v_y;
    for (int p = 1; p <= max_p; p++) {
        init_v_distr.at(p) = init_v_distr.at(p) / c_norm / v_ph_vol;
    }
}

int VMesh::indexToP(Index2d i) {
    return index_to_p.at(i.ix, i.iy);
}

int VMesh::safeIndexToP(Index2d i) {
    if ((i.ix < 1) or (i.iy < 1) or (i.ix > n_v_x) or (i.iy > n_v_y))
        return 0;
    return index_to_p.at(i.ix, i.iy);
}

int VMesh::maxP() {return max_p;}

OneDArr VMesh::getInitVDistr() {return init_v_distr;}

double VMesh::getVCut() {return v_cut;}


double VMesh::gridLockX(double v_x) {
    int index = static_cast<int>(floor(n_v_x / 2 * (1 + v_x / v_cut))) + 1;
    if (index < 1) 
        return v_x_mesh.at(1);
    else if (index > n_v_x) 
        return v_x_mesh.at(n_v_x);
    else 
        return v_x_mesh.at(index);
}

double VMesh::gridLockY(double v_y) {
    int index = static_cast<int>(floor(n_v_y / 2 * (1 + v_y / v_cut))) + 1;
    if (index < 1) 
        return v_y_mesh.at(1);
    else if (index > n_v_y) 
        return v_y_mesh.at(n_v_y);
    else 
        return v_y_mesh.at(index);
}

Index2d VMesh::getClosestVelocityIndex(Vector2d v) {
    int ix = static_cast<int>(std::round((v.x + v_cut) / (2 * v_cut) * n_v_x + 0.5));
    int iy = static_cast<int>(std::round((v.y + v_cut) / (2 * v_cut) * n_v_y + 0.5));
    return Index2d(ix, iy);
}

Index2d VMesh::getFloorVelocityIndex (Vector2d v) {
    int ix = static_cast<int>(std::floor((v.x + v_cut) / (2 * v_cut) * n_v_x + 0.5));
    int iy = static_cast<int>(std::floor((v.y + v_cut) / (2 * v_cut) * n_v_y + 0.5));
    return Index2d(ix, iy);
}

double VMesh::getNodeEnergy(Index2d node, Vector2d v_c) {
    double v_node_x = v_x_mesh.at(node.ix);
    double v_node_y = v_y_mesh.at(node.iy);
    Vector2d delta_v(v_node_x - v_c.x, v_node_y - v_c.y);
    return delta_v.pow2();
}