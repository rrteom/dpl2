#ifndef COLLISION_NODES_HPP
#define COLLISION_NODES_HPP

#include "VMesh.hpp"
#include "structures.hpp"

class CollisionNodes {
private:
    VMesh* p_v_mesh;

public:
    double s_max;
    int n_nodes;
    std::vector<int> int_p_alpha, int_p_beta, int_p_l, int_p_ls, int_p_m, int_p_ms;
    std::vector<double> collisions, rel_velocities, thetas, interp_r;
    std::vector<bool> is_active;

    CollisionNodes(int n_nodes, VMesh* p_v_mesh, double s_max);
    void randomizeNodes(std::vector<unsigned int> korobov_coefs);
    void scaleValues();
    void checkOutOfSphere();
    void calculateRelVelocitiesAfterCollision();
    void findInterpNodes();
    bool isActive(int coll_no);
    std::vector<double> getCollisions();
    void prepareNodes(std::vector<unsigned int> korobov_coefs);
    std::vector<int> generatePermutation();
};

#endif