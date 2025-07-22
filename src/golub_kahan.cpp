#include "golub_kahan.hpp"
#include <Eigen/src/Core/Matrix.h>
#include <Eigen/src/Core/util/Constants.h>
#include <algorithm>
#include <vector>
#include <iostream>
#include <Eigen/Dense>

bool isDiagonal(Eigen::MatrixXd B){

    int numRows = B.rows();
    int numCols = B.cols();

    bool diagonal = false;
    double epsilon = 1e-12;
    for (int i = 0; i < numRows; i++){
            
        if (i + 1 > numCols){
            continue;
        }

        if (B(i, i + 1) > epsilon){
            return diagonal;
        }
    }

    return true;
}

Eigen::MatrixXd chooseSubmatrix(Eigen::MatrixXd B){

    return B.block(B.rows() - 2, B.cols() - 2, 2, 2);
}

double wilkinsonShift(Eigen::MatrixXd trail){

    double a = trail(0, 0);
    double b = trail(0, 1);
    double c = trail(1, 1);

    double delta = (a - c) / 2;    
    double denominator = std::sqrt(std::pow(delta, 2) + std::pow(b, 2)) + std::abs(delta);

    double sign = (delta >= 0) ? 1.0 : 0.0;
    return  c - (sign * std::pow(b, 2)) / denominator;
}

Eigen::MatrixXd shiftHelper(Eigen::MatrixXd A, double miu){

    double first = std::pow(A(0, 0), 2) - miu;
    double second = A(0, 0) *A(1, 0);

    double radius = std::sqrt(std::pow(first, 2) + std::pow(second, 2));
    double cos = first / radius;
    double sin = second / radius;

    Eigen::MatrixXd G = Eigen::MatrixXd(2, 2);
    G(0,0) = cos;
    G(0, 1) = sin;
    G(1, 0) = -sin;
    G(1, 1) = cos;

    return G;
}

Eigen::MatrixXd applyShift(Eigen::MatrixXd trail, double miu, Eigen::MatrixXd &U, Eigen::MatrixXd &V_transposed){

    Eigen::MatrixXd G = shiftHelper(trail, miu);
    U = G * U;
    Eigen::MatrixXd newTrail = trail * G.transpose();
    
    G = shiftHelper(newTrail, miu);
    V_transposed = V_transposed * U;
    return G.transpose() * newTrail;
}

Eigen::MatrixXd golubKahan(Eigen::MatrixXd U, Eigen::MatrixXd &B, Eigen::MatrixXd V_transposed){

    B = B * B.transpose();
    while (!isDiagonal(B)){

        Eigen::MatrixXd trail = chooseSubmatrix(B);
        double miu = wilkinsonShift(trail);

        // the trail doesn't need to be another variable
        // and the result needs to be updated into B;
        applyShift(trail, miu, U, V_transposed);
    }
}