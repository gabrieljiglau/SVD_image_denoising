#include "golub_kahan.hpp"
#include <Eigen/src/Core/Matrix.h>
#include <Eigen/src/Core/util/Constants.h>
#include <algorithm>
#include <vector>
#include <iostream>
#include <Eigen/Dense>


bool isDiagonal(Eigen::MatrixXd B, double epsilon){

    int numRows = B.rows();
    int numCols = B.cols();

    bool diagonal = false;
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

void deflateValues(Eigen::MatrixXd &B, double epsilon){

    int numRows = B.rows();
    int numCols = B.cols();

    for (int i = 0; i < numRows; i++){

        if (B(i, i) < epsilon){
            B(i, i) = 0;
        }

        if (i + 1 > numCols){
            continue;
        }

        if (B(i, i + 1) < epsilon){
            B(i, i + 1) = 0;
        }
    }

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

Eigen::MatrixXd shiftHelper(Eigen::MatrixXd &A, double miu, int i, double epsilon){

    double first = std::pow(A(i, i), 2) - miu;

    if (first < epsilon){
        first = 0;
        return Eigen::MatrixXd::Identity(A.rows(), A.cols());
    }

    double second = A(i, i) *A(i + 1, i);

    double radius = std::sqrt(std::pow(first, 2) + std::pow(second, 2));
    double cos = first / radius;
    double sin = second / radius;  

    int k = A.rows();
    Eigen::MatrixXd G = Eigen::MatrixXd::Identity(k, k);
    G(k - 2,k - 2) = cos;
    G(k - 2, k - 1) = sin;
    G(k - 1, k - 2) = -sin;
    G(k - 1, 1) = cos;

    return G;
}

// shiftul e altfel
Eigen::MatrixXd applyShift(Eigen::MatrixXd &B, Eigen::MatrixXd &U, Eigen::MatrixXd &V_transposed, double miu, double epsilon){

    Eigen::MatrixXd G;
    Eigen::MatrixXd newB;
    for (int i = 0; i < B.rows() - 2; i += 2){
        G = shiftHelper(B, miu, i, epsilon);
        U = G * U;
        newB = B * G.transpose();
        
        G = shiftHelper(newB, miu, i, epsilon);
        V_transposed = V_transposed * G;
    }

    return G.transpose() * newB;
}

Eigen::MatrixXd golubKahan(Eigen::MatrixXd U, Eigen::MatrixXd &B, Eigen::MatrixXd V_transposed){

    double epsilon = 1e-12;
    Eigen::MatrixXd newB;

    while (!isDiagonal(B, epsilon)){

        Eigen::MatrixXd trail = chooseSubmatrix(B);
        double miu = wilkinsonShift(trail);

        newB = applyShift(B, U, V_transposed, miu, epsilon);
        B = newB;
        
        deflateValues(B, epsilon);
    }

    return B;
}