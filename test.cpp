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
    
    CollisionNodes nodes(50021, &v_mesh, 1);

    nodes.randomizeNodes(koefs);
    nodes.scaleValues();
    nodes.checkOutOfSphere();
    nodes.calculateRelVelocitiesAfterCollision();
    nodes.findInterpNodes();

    return 0;
}