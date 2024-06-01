#include <iostream>
#include "structures.hpp"
#include "VMesh.hpp"
#include "distribution.hpp"
#include "CollisionNodes.hpp"
#include <fstream>

using namespace std;

int main() {
    int n_v_x = 20, n_v_y = 20;
    // std::vector<unsigned int> koefs = {1, 11281, 7537, 39218, 32534, 11977};
    VMesh v_mesh(4.8, n_v_x, n_v_y, 1, 2);
    
    // CollisionNodes nodes(50021, &v_mesh, 1);
    // nodes.prepareNodes(koefs);
    Distribution distribution(50, 50, 20, 20,
                 4.8, n_v_x, n_v_y,
                 1, 2, 1,
                 8, 8, 13, 13);


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

    return 0;
}