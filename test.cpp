#include <iostream>
#include "structures.hpp"
#include "VMesh.hpp"
#include "distribution.hpp"
#include "CollisionNodes.hpp"
#include <fstream>

using namespace std;

void logging(Distribution);
void logging_append(Distribution);
void logging_v11(Distribution d, VMesh v_m, int n_v_x, int n_v_y, ThreeDArr);

int main() {
    int n_v_x = 20, n_v_y = 20;
    std::vector<unsigned int> koefs = {1, 832685, 1269182, 1228431, 532894, 174792};
    VMesh v_mesh(4.8, n_v_x, n_v_y, 2, 1);
    Distribution distribution(50, 50, 30, 30,
                4.8, n_v_x, n_v_y,
                2, 1,
                11, 11, 20, 20);
    CollisionNodes collision_nodes(2000003, &v_mesh, 1.0d);
    // std::cout << collision_nodes.s_max;
    collision_nodes.prepareNodes(koefs);

    double tau_0 = distribution.getTau0();
    double tau_x = tau_0 / 4, tau_y = tau_0 / 2; 
    cout << "tau_0: " << tau_0 << endl;
    ThreeDArr init_distr = distribution.distr;
    logging(distribution);
    logging_v11(distribution, v_mesh, n_v_x, n_v_y, init_distr);
    for (int time_step = 0; time_step < 100; time_step++) {
        cout << "time_step: " << time_step << endl;
        //space step / 2
        distribution.stepX(tau_x);
        distribution.stepY(tau_y);
        distribution.stepX(tau_x);

        //collision step
        distribution.collisionStep(collision_nodes, tau_0);

        //space step / 2
        distribution.stepX(tau_x);
        distribution.stepY(tau_y);
        distribution.stepX(tau_x);  

        //save      
        logging_append(distribution);
    }

    return 0;
}

void logging(Distribution distribution) {
    TwoDArr n = distribution.getConcentration();
    ofstream n_file;
    n_file.open("intgr_test/n.log");
    for (int i = 1; i <= n.dimx_; i++) {
        for (int j = 1; j <= n.dimy_; j++) {
            n_file << n.at(i, j) << '\t';
        }
        n_file << endl;
    }
    n_file << endl;
    n_file.close();

    TwoDArr tt = distribution.getTemperature();
    ofstream tt_file;
    tt_file.open("intgr_test/tt.log");
    for (int i = 1; i <= n.dimx_; i++) {
        for (int j = 1; j <= n.dimy_; j++) {
            tt_file << tt.at(i, j) << '\t';
        }
        tt_file << endl;
    }
    tt_file << endl;
    tt_file.close();    
}

void logging_append(Distribution distribution) {
    TwoDArr n = distribution.getConcentration();
    ofstream n_file;
    n_file.open("intgr_test/n.log", ios_base::app);
    for (int i = 1; i <= n.dimx_; i++) {
        for (int j = 1; j <= n.dimy_; j++) {
            n_file << n.at(i, j) << '\t';
        }
        n_file << endl;
    }
    n_file << endl;
    n_file.close();

    TwoDArr tt = distribution.getTemperature();
    ofstream tt_file;
    tt_file.open("intgr_test/tt.log", ios_base::app);
    for (int i = 1; i <= n.dimx_; i++) {
        for (int j = 1; j <= n.dimy_; j++) {
            tt_file << tt.at(i, j) << '\t';
        }
        tt_file << endl;
    }
    tt_file << endl;
    tt_file.close();    
}

void logging_v11(Distribution d, VMesh v_m, int n_v_x, int n_v_y, ThreeDArr init_distr) {
    double curr;
    ofstream file;
    file.open("intgr_test/v11.log", ios_base::app);
    for (int a_x = 1; a_x <= n_v_x; a_x++) {
        for (int a_y = 1; a_y <= n_v_y; a_y++) {
            int p_index = v_m.indexToP(Index2d(a_x, a_y));
            if (p_index == 0) {
                file << 0 << ' ';
            }
            else {
                file << d.distr.at(1, 1, p_index) - init_distr.at(1, 1, p_index) << ' ';
            }
        }
        file << endl;
    }
    file << endl;
    file.close();
}