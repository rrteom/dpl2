#include "CollisionNodes.hpp"

CollisionNodes::CollisionNodes(int n_nodes, VMesh* p_v_mesh, double s_max_) : n_nodes(n_nodes), p_v_mesh(p_v_mesh), s_max(s_max_) {
    collisions = std::vector<double>(6 * n_nodes);
    rel_velocities = std::vector<double>(2 * n_nodes);
    thetas = std::vector<double>(n_nodes);
    is_active = std::vector<bool>(n_nodes);
    interp_r = std::vector<double>(n_nodes);
    int_p_alpha = std::vector<int>(n_nodes);
    int_p_beta = std::vector<int>(n_nodes);
    int_p_l = std::vector<int>(n_nodes);
    int_p_ls = std::vector<int>(n_nodes);
    int_p_m = std::vector<int>(n_nodes);
    int_p_ms = std::vector<int>(n_nodes);
}

std::vector<double> CollisionNodes::getCollisions() {return collisions;}

void CollisionNodes::randomizeNodes(std::vector<unsigned int> korobov_coefs) {
    double trashed_int_part;
    for (int curr_node = 0; curr_node < n_nodes; curr_node++) {
        for (int i = 0; i < 6; i++) {
            collisions[curr_node * 6 + i] = modf(korobov_coefs[i] * (curr_node + 1) / static_cast<double>(n_nodes), &trashed_int_part);
            // collisions[curr_node * 8 + i] = std::rand() / static_cast<double>(RAND_MAX); // uniform [0, 1]
        }
    }
}

void CollisionNodes::scaleValues() {
    double v_cut = p_v_mesh->getVCut();
    for (int curr_node = 0; curr_node < n_nodes; curr_node++) {
        for (int i = 0; i < 6; i++) {
            if ((i == 0) or (i == 2)) { // v_x scaling
                collisions[curr_node * 6 + i] = p_v_mesh->gridLockX(2 * v_cut * collisions[curr_node * 6 + i] - v_cut);
            }
            else if ((i == 1) or (i == 3)) { // v_y scaling
                collisions[curr_node * 6 + i] = p_v_mesh->gridLockY(2 * v_cut * collisions[curr_node * 6 + i] - v_cut);
            }
            else if (i == 4) { // S scaling
                collisions[curr_node * 6 + i] *= s_max;
            }
            else { // sign of rotation
                double sign_val = (collisions[curr_node * 6 + i] > 0.5) ? 1: -1;
                collisions[curr_node * 6 + i] = sign_val;
            }
        }
    }
}

void CollisionNodes::checkOutOfSphere() {
    double v_cut_2 = pow(p_v_mesh->getVCut(), 2);
    for (int coll_no = 0; coll_no < n_nodes; coll_no++) {
        double v_2, v1_2;
        v_2 = pow(collisions[6 * coll_no], 2) + pow(collisions[6 * coll_no + 1], 2);
        v1_2 = pow(collisions[6 * coll_no + 2], 2) + pow(collisions[6 * coll_no + 3], 2);
        if ((v_2 <= v_cut_2) and (v1_2 <= v_cut_2)) {
            is_active[coll_no] = true;
        }
        else {
            is_active[coll_no] = false;
        }
    }
}

bool CollisionNodes::isActive(int coll_no) {return is_active[coll_no];}

void CollisionNodes::calculateRelVelocitiesAfterCollision() {
    double g_x, g_y;
    for (int coll_no = 0; coll_no < n_nodes; coll_no++) {
        g_x = collisions[6 * coll_no + 2] - collisions[6 * coll_no];
        g_y = collisions[6 * coll_no + 3] - collisions[6 * coll_no + 1];
        if ((g_x == 0) and (g_y == 0)) {
            is_active[coll_no] = false;
        }
        else if (is_active[coll_no]) {
            thetas[coll_no] = 2 * acos(sqrt(collisions[coll_no * 6 + 4]));
            // g`_x = g_x * cos(theta) - sign * g_y * sin(theta)
            rel_velocities[coll_no * 2] = g_x * cos(thetas[coll_no]) - collisions[6 * coll_no + 5] * g_y * sin(thetas[coll_no]);
            // g`_x = g_y * cos(theta) + sign * g_x * sin(theta)
            rel_velocities[coll_no * 2 + 1] = g_y * cos(thetas[coll_no]) + collisions[6 * coll_no + 5] * g_x * sin(thetas[coll_no]);
        }
    }
    return;
}

void CollisionNodes::findInterpNodes() {
    for(int coll_no = 0; coll_no < n_nodes; coll_no++) {
        bool nodes_found = false;
        if (!is_active[coll_no])
            continue;
        Vector2d v_alpha = Vector2d(collisions[6 * coll_no], collisions[6 * coll_no + 1]);
        Vector2d v_beta = Vector2d(collisions[6 * coll_no + 2], collisions[6 * coll_no + 3]);
        Index2d alpha = p_v_mesh->getClosestVelocityIndex(v_alpha);
        Index2d beta = p_v_mesh->getClosestVelocityIndex(v_beta);

        Vector2d g_rel = Vector2d(rel_velocities[coll_no * 2], rel_velocities[coll_no * 2 + 1]);
        Vector2d v_c = Vector2d((v_alpha.x + v_beta.x) / 2, (v_alpha.y + v_beta.y) / 2);
        Vector2d v_alpha_new = Vector2d((v_alpha.x + v_beta.x - g_rel.x) / 2, (v_alpha.y + v_beta.y - g_rel.y) / 2);
        Index2d eta_near = p_v_mesh->getClosestVelocityIndex(v_alpha_new);
        Index2d eta_near_1(alpha.ix + beta.ix - eta_near.ix, alpha.iy + beta.iy - eta_near.iy);
        if ((p_v_mesh->safeIndexToP(eta_near) == 0) or (p_v_mesh->safeIndexToP(eta_near_1) == 0)) {
            is_active[coll_no] = false;
            continue;
        }
        double energy_near = p_v_mesh->getNodeEnergy(eta_near, v_c);

        Index2d eta_floor = p_v_mesh->getFloorVelocityIndex(v_alpha_new);
        double energy_true = Vector2d(v_alpha_new.x - v_c.x, v_alpha_new.y - v_c.y).pow2();

        double r_interp;
        Index2d l, ls, m, ms;
        if (energy_near == energy_true) {
            l = eta_near;
            ls = eta_near;
            m = eta_near_1;
            ms = eta_near_1;
            r_interp = 1;
            nodes_found = true;
        }
        else if (energy_near > energy_true) {
            ls = eta_near;
            ms = eta_near_1;
            double curr_min_distance, energy_eta_res;
            for (int i = 0; i <= 1; i++) {
                for (int j = 0; j <= 1; j++) {
                    Index2d current_eta = Index2d(eta_floor.ix + i, eta_floor.iy + j);
                    Index2d current_eta_1 = Index2d(alpha.ix + beta.ix - current_eta.ix, alpha.iy + beta.iy - current_eta.iy);
                    if (p_v_mesh->safeIndexToP(current_eta) == 0)
                        continue;
                    double current_energy = p_v_mesh->getNodeEnergy(current_eta, v_c);
                    if ((p_v_mesh->safeIndexToP(current_eta) != 0) and (p_v_mesh->safeIndexToP(current_eta_1) != 0) and (current_energy < energy_true)) {
                        double curr_eta_distance = p_v_mesh->getNodeEnergy(current_eta, v_alpha_new);
                        if (!nodes_found) {
                            nodes_found = true;
                            curr_min_distance = curr_eta_distance;
                            energy_eta_res = current_energy;
                            l = current_eta;
                            m = current_eta_1;
                        }
                        else if (curr_eta_distance < curr_min_distance){
                            curr_min_distance = curr_eta_distance;
                            energy_eta_res = current_energy;
                            l = current_eta;
                            m = current_eta_1;
                        }
                    }
                }
            }
            if (!nodes_found) {
                is_active[coll_no] = false;
                continue;
            }
            r_interp = (energy_true - energy_eta_res) / (energy_near - energy_eta_res);
        }
        else {
            // energy_near < energy_true
            l = eta_near;
            m = eta_near_1;
            double curr_min_distance, energy_eta_res;
            for (int i = 0; i <= 1; i++) {
                for (int j = 0; j <= 1; j++) {
                    Index2d current_eta = Index2d(eta_floor.ix + i, eta_floor.iy + j);
                    Index2d current_eta_1 = Index2d(alpha.ix + beta.ix - current_eta.ix, alpha.iy + beta.iy - current_eta.iy);
                    if (p_v_mesh->safeIndexToP(current_eta) == 0)
                        continue;
                    double current_energy = p_v_mesh->getNodeEnergy(current_eta, v_c);
                    if ((p_v_mesh->safeIndexToP(current_eta) != 0) and (p_v_mesh->safeIndexToP(current_eta_1) != 0) and (current_energy > energy_true)) {
                        double curr_eta_distance = p_v_mesh->getNodeEnergy(current_eta, v_alpha_new);
                        if (!nodes_found) {
                            nodes_found = true;
                            curr_min_distance = curr_eta_distance;
                            energy_eta_res = current_energy;
                            ls = current_eta;
                            ms = current_eta_1;
                        }
                        else if (curr_eta_distance < curr_min_distance){
                            curr_min_distance = curr_eta_distance;
                            energy_eta_res = current_energy;
                            ls = current_eta;
                            ms = current_eta_1;
                        }
                    }
                }
            }
            if (!nodes_found) {
                is_active[coll_no] = false;
                continue;
            }
            r_interp = (energy_true - energy_near) / (energy_eta_res - energy_near);
        }
        if (!nodes_found) {
            is_active[coll_no] = false;
            continue;
        }
        interp_r[coll_no] = r_interp;
        int_p_alpha[coll_no] = p_v_mesh->indexToP(alpha);
        int_p_beta[coll_no] = p_v_mesh->indexToP(beta);
        int_p_l[coll_no] = p_v_mesh->indexToP(l);
        int_p_ls[coll_no] = p_v_mesh->indexToP(ls);
        int_p_m[coll_no] = p_v_mesh->indexToP(m);
        int_p_ms[coll_no] = p_v_mesh->indexToP(ms);
    }
    return;
}

void CollisionNodes::prepareNodes(std::vector<unsigned int> korobov_coefs) {
    randomizeNodes(korobov_coefs);
    scaleValues();
    checkOutOfSphere();
    calculateRelVelocitiesAfterCollision();
    findInterpNodes();
    return;
}

std::vector<int> CollisionNodes::generatePermutation() {
    std::vector<int> initial, result;
    for (int i = 0; i < n_nodes; i++)
        initial.push_back(i);
    for (int i = 0; i < n_nodes; i++) {
        int rand_index = rand() % initial.size();
        result.push_back(initial[rand_index]);
        std::swap(initial[rand_index], initial[initial.size() - 1]);
        initial.pop_back();
    }
    return result;
}