#include <iostream>
#include "structures.hpp"
#include "VMesh.hpp"
#include "distribution.hpp"
#include "CollisionNodes.hpp"
#include <fstream>

using namespace std;

void logging(Distribution);
void logging_append(Distribution);

int main() {
    int n_v_x = 20, n_v_y = 20;
    // std::vector<unsigned int> koefs = {1, 11281, 7537, 39218, 32534, 11977};
    VMesh v_mesh(4.8, n_v_x, n_v_y, 1, 2);
    
    // CollisionNodes nodes(50021, &v_mesh, 1);
    // nodes.prepareNodes(koefs);
    Distribution distribution(50, 50, 20, 20,
                 4.8, n_v_x, n_v_y,
                 1, 2, 1,
                 6, 6, 15, 15);

    double tau_0 = distribution.getTau0();
    double tau_x = tau_0 / 2, tau_y = tau_0; 
    cout << "tau_0: " << tau_0 << endl;

    logging(distribution);
    for (int time_step = 0; time_step < 300; time_step++) {
        distribution.stepX(tau_x);
        // logging_append(distribution);
        distribution.stepY(tau_y);
        // logging_append(distribution);
        distribution.stepX(tau_x);
        // distribution.step1DY(tau_y);
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