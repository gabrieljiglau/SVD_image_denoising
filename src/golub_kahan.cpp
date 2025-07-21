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

double wilkinsonShift(Eigen::MatrixXd B){

    double a = B(0, 0);
    double b = B(0, 1);
    double c = B(1, 1);

    double delta = (a - c) / 2;
    double sign = (delta > 0) ? 1.0 : -1.0;
    
    double denominator = std::sqrt(std::pow(delta, 2) + std::pow(b, 2) + std::abs(delta));
    return c - (sign * std::pow(b, 2)) / denominator;
}

Eigen::MatrixXd golubKahan(Eigen::MatrixXd U, Eigen::MatrixXd &B, Eigen::MatrixXd V_transposed){

    while (!isDiagonal(B)){

        B = B * B.transpose();

        Eigen::MatrixXd trail = chooseSubmatrix(B);
        double miu = wilkinsonShift(B);

        applyShift(B, miu);
    }
}