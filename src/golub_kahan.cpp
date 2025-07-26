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

        if (i + 1 >= numCols){
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
    double denominatorSquared = std::pow(delta, 2) + std::pow(b, 2);

    if (denominatorSquared < 1e-20){ // avoid division by 0
        return c;
    }
  
    double sign = (delta >= 0) ? 1.0 : 0.0;
    return  c - (sign * std::pow(b, 2)) / (std::abs(delta) + std::sqrt(denominatorSquared));
}

Eigen::MatrixXd givensCommon(int numRows, int numCols, int i, double a, double b){

    Eigen::MatrixXd G = Eigen::MatrixXd::Identity(numRows, numCols);

    double radius = std::hypot(a, b);
    if (radius == 0){
        return G;
    }

    double cos = a / radius;
    double sin = -b / radius;

    G(i, i) = cos;
    G(i, i + 1) = sin;
    G(i + 1, i) = -sin;
    G(i + 1, i + 1) = cos;

    return G;
}

Eigen::MatrixXd applyShift(Eigen::MatrixXd B, Eigen::MatrixXd &U, Eigen::MatrixXd &V_transposed, double miu, double epsilon){

    double a = 0.0;
    double b = 0.0;
    Eigen::MatrixXd G = Eigen::MatrixXd::Identity(B.rows(), B.cols());

    assert(B.rows() == B.cols() && "The input matrix must be square");
    for (int i = 0; i < B.rows() - 1; i++){

        std::cout << "Now applying shift for row " << i << std::endl;

        if (i != 0){ 
            a = B(i, i);

            if (i + 1 > B.cols() - 1){
                break;
            }
            b = B(i, i+1);
        } else { // smart miu
            a = std::pow(B(0, 0), 2) - miu;
            b = B(0, 0) * B(0, 1);
        }
        
        G = givensCommon(B.rows(), B.cols(), i, a, b);
        V_transposed = V_transposed * G;
        B = B * G.transpose();
        
        a = B(i, i);

        if (i + 1 >= B.rows() - 1){
            break;
        }
        b = B(i + 1, i);
        G = givensCommon(B.rows(), B.cols(), i, a, b);

        U = G * U;
        B = G * B;

        std::cout << B << std::endl;
    }

    return B;
}

std::vector<Eigen::MatrixXd> golubKahan(Eigen::MatrixXd B, Eigen::MatrixXd &U, Eigen::MatrixXd &V_transposed){

    double epsilon = 1e-12;
    int iterationNumber = 0;

    while (!isDiagonal(B, epsilon)){

        std::cout<< "Now at iteration " << iterationNumber << std::endl;

        Eigen::MatrixXd trail = chooseSubmatrix(B);
        double miu = wilkinsonShift(trail);

        // miu are mereu aceeasi valoare ???
        std::cout << "miu = " << miu << std::endl;

        B = applyShift(B, U, V_transposed, miu, epsilon);

        deflateValues(B, epsilon);
        iterationNumber += 1;

        std::cout << "B : " << std::endl;
        std::cout << B << std::endl;

        std::cout << "U : " << std::endl;
        std::cout << U << std::endl;

        std::cout << "V_transposed : " << std::endl;
        std::cout << V_transposed << std::endl;
        
        if (iterationNumber > 25){
            break;
        }
    }

    return std::vector{B, U, V_transposed};
}