#include <iostream>
#include "structures.hpp"
#include "VMesh.hpp"
#include "distribution.hpp"
#include "CollisionNodes.hpp"

using namespace std;

int main() {
    int n_v_x = 6, n_v_y = 8;
    // std::vector<unsigned int> koefs = {1, 11281, 7537, 39218, 32534, 11977};
    // std::vector<double> collisions;
    VMesh v_mesh(4.8, n_v_x, n_v_y, 1, 2);

    for (int ix = 0; ix <= n_v_x; ix++) {
        for (int iy = 0; iy <= n_v_y; iy++) {
            cout << v_mesh.safeIndexToP(Index2d(ix, iy)) << ' ';
        }
        cout << endl;
    }
    
    // CollisionNodes nodes(50021, &v_mesh, 1);

    // collisions = nodes.getCollisions();
    // for (int i = 0; i < collisions.size(); i++) {
    //     cout << collisions[i] << ' ';
    //     if (i % 6 == 5)
    //         cout << endl;
    // }

    // nodes.randomizeNodes(koefs);

    // collisions = nodes.getCollisions();
    // for (int i = 0; i < collisions.size(); i++) {
    //     cout << collisions[i] << ' ';
    //     if (i % 6 == 5)
    //         cout << endl;
    //     if (i == 59) 
    //         break;
    // }

    // nodes.scaleValues();
    // cout << "scaled collisions" << endl;
    // collisions = nodes.getCollisions();
    // for (int i = 0; i < collisions.size(); i++) {
    //     cout << collisions[i] << ' ';
    //     if (i % 6 == 5)
    //         cout << endl;
    //     if (i == 59) 
    //         break;
    // }


    return 0;
}