#include <cmath>
#include <vector>
#include "structures.hpp"

using namespace std;

Vector2d::Vector2d() : x(0), y(0) {}
Vector2d::Vector2d(double x, double y) : x(x), y(y) {}
double Vector2d::pow2() {return x*x + y*y;}


Index2d::Index2d() : ix(0), iy(0) {}
Index2d::Index2d(int ix, int iy) : ix(ix), iy(iy) {}

OneDArr::OneDArr () : dim_(0) {}
OneDArr::OneDArr (int n) : dim_ (n) {
    data_.resize(n);
}
OneDArr::OneDArr (double start_val, double end_val, int n_points) : dim_ (n_points){
    // mesh constructor
    double step = (end_val - start_val) / n_points;
    for (int i = 0; i < n_points; i++) {
        data_.push_back(start_val + step * (i + 0.5));
    }
}
void OneDArr::resize(int dim) {
    data_.resize(dim);
    dim_ = dim;
    return;
}
void OneDArr::push_back(double val) {
    dim_++;
    data_.push_back(val);
    return;
}
double& OneDArr::at(int i) {
    if (i > dim_ || i < 1)
        throw std::out_of_range("matrix indices out of range");
    return data_[i - 1];
}
int OneDArr::size() {return dim_;}


TwoDArr::TwoDArr () : dimx_ (0), dimy_ (0) {}
TwoDArr::TwoDArr (int dimx, int dimy) : dimx_ (dimx), dimy_ (dimy) {
    data_.resize(dimx_ * dimy_);
}
void TwoDArr::resize(int dimx, int dimy) {
    dimx_ = dimx;
    dimy_ = dimy;
    data_.resize(dimx * dimy);
}
double& TwoDArr::at(int x,  int y) {
    if (x > dimx_ || y > dimy_ || x < 1 || y < 1)
        throw std::out_of_range("matrix indices out of range");
    return data_[dimx_ * (y - 1) + x - 1];
}
void TwoDArr::fill(double val) {
    std::fill(data_.begin(), data_.end(), val);
    return;
}

ThreeDArr::ThreeDArr () : dimx_ (0), dimy_ (0), dimp_ (0) {}
ThreeDArr::ThreeDArr (int dimx, int dimy, int dimp_) : dimx_ (dimx), dimy_ (dimy), dimp_(dimp_) {
    data_.resize(dimx_ * dimy_ * dimp_);
}
void ThreeDArr::resize(int dimx, int dimy, int dimp) {
    dimx_ = dimx;
    dimy_ = dimy;
    dimp_ = dimp;
    data_.resize(dimx * dimy * dimp);
}
double& ThreeDArr::at(int x,  int y, int p) {
    if (x > dimx_ || y > dimy_ || p > dimp_ || x < 1 || y < 1 || p < 1)
        throw std::out_of_range("matrix indices out of range");
    return data_[ dimy_ * dimx_ * (p - 1) + dimx_ * (y - 1) + x - 1];
}
void ThreeDArr::fill(double val) {
    std::fill(data_.begin(), data_.end(), val);
    return;
}
