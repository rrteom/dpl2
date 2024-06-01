#include <iostream>
#include "structures.hpp"
#include "VMesh.hpp"
#include "distribution.hpp"
#include "CollisionNodes.hpp"
#include <fstream>

using namespace std;

int main() {
    int n_v_x = 20, n_v_y = 20;
    std::vector<unsigned int> koefs = {1, 11281, 7537, 39218, 32534, 11977};
    VMesh v_mesh(4.8, n_v_x, n_v_y, 1, 2);


    // ofstream itop_file;
    // itop_file.open("col_test/index_to_p");
    // for (int ix = 1; ix <= n_v_x; ix++) {
    //     for (int iy = 1; iy <= n_v_y; iy++) {
    //         itop_file << v_mesh.safeIndexToP(Index2d(ix, iy)) << '\t';
    //     }
    //     itop_file << endl;
    // }
    // itop_file.close();
    
    CollisionNodes nodes(50021, &v_mesh, 1);

    nodes.randomizeNodes(koefs);
    nodes.scaleValues();
    nodes.checkOutOfSphere();
    nodes.calculateRelVelocitiesAfterCollision();
    
    nodes.findInterpNodes();

    vector<double> collisions = nodes.collisions;
    ofstream col_file;
    col_file.open("col_test/cols");
    col_file << "v_a_x\tv_a_y\tv_b_x\tv_b_y\tb\teps\ttrash" << endl; 
    for (int i = 0; i < collisions.size(); i++) {
        col_file << collisions[i] << '\t';
        if (i % 6 == 5)
            col_file << endl;
    }
    col_file.close();

    vector<double> rel_velocities = nodes.rel_velocities;
    vector<double> thetas = nodes.thetas;
    ofstream relv_file;
    col_file.open("col_test/relv_theta");
    col_file << "rel_v_x_new\trel_v_y_new\ttheta" << endl; 
    for (int i = 0; i < thetas.size(); i++) {
        col_file << rel_velocities[2 * i] << '\t' << rel_velocities[2 * i + 1] << '\t' << thetas[i] << endl;
    }
    col_file.close();

    std::vector<int> int_p_alpha = nodes.int_p_alpha, int_p_beta = nodes.int_p_beta, 
                        int_p_l = nodes.int_p_l, int_p_ls = nodes.int_p_ls, int_p_m = nodes.int_p_m, int_p_ms = nodes.int_p_ms;
    std::vector<bool> is_active = nodes.is_active;
    std::vector<double> interp_r = nodes.interp_r;
    ofstream interp_file;  
    interp_file.open("col_test/interp_active");
    interp_file << "is_active\tp_alpha\tp_beta\tp_l\tp_ls\tp_m\tp_ms\tinterp_r" << endl; 
    for (int i = 0; i < is_active.size(); i++) {
        interp_file << is_active[i] << '\t';
        interp_file << int_p_alpha[i] << '\t';
        interp_file << int_p_beta[i] << '\t';
        interp_file << int_p_l[i] << '\t';
        interp_file << int_p_ls[i] << '\t';
        interp_file << int_p_m[i] << '\t';
        interp_file << int_p_ms[i] << '\t';
        interp_file << interp_r[i];
        interp_file << endl;
    }
    interp_file.close();  




    return 0;
}