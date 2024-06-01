#ifndef STRUCTURES_HPP
#define STRUCTURES_HPP

#include <cmath>
#include <vector>
#include <stdexcept>

class Vector2d {
public:
    Vector2d();
    Vector2d(double x, double y);
    double pow2();
    double x, y;
};

class Index2d {
public:
    Index2d();
    Index2d(int ix, int iy);
    int ix, iy;
};

class OneDArr {
    // idx starts from 1
    std::vector<double> data_;
    int dim_;
public:
    OneDArr ();
    OneDArr (int n);
    OneDArr (double start_val, double end_val, int n_points);
    void resize(int dim);
    void push_back(double val);
    double& at(int i);
    int size();
};

class TwoDArr {
    // idx starts from 1
    std::vector<double> data_;
    
public:
    int dimx_, dimy_;
    TwoDArr();
    TwoDArr (int dimx, int dimy);
    void resize(int dimx, int dimy);
    double& at(int x,  int y);
    void fill(double val);
};

class ThreeDArr {
    // idx starts from 1
    std::vector<double> data_;
public:
    int dimx_, dimy_, dimp_;
    ThreeDArr();
    ThreeDArr (int dimx, int dimy, int dimp);
    void resize(int dimx, int dimy, int dimp);
    double& at(int x,  int y, int p);
    void fill(double val);
};

#endif