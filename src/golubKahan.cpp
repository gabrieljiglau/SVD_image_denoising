#include "golubKahan.hpp"
#include <Eigen/src/Core/Matrix.h>
#include <Eigen/src/Core/util/Constants.h>
#include <vector>
#include <iostream>
#include <Eigen/Dense>


bool isDiagonal(Eigen::MatrixXd B){

    int numRows = B.rows();
    int numCols = B.cols();

    for (int i = 0; i < numRows; i++){
        for (int j = 0; j < numCols; j++){
            if (i != j && B(i, j) != 0){
                return false;
            }

            if (i == j && B(i, j) == 0){
                return false;
            }

            if (i == j && B(i, j) != 0){
                continue;
            }
        }
    }

    return true;
}

void deflateValues(Eigen::MatrixXd &B, double epsilon){

    int numRows = B.rows();
    int numCols = B.cols();

    for (int i = 0; i < numRows; i++){
        for (int j = 0; j < numCols; j++){

            if (std::abs(B(i, j)) < epsilon){
                B(i, j) = 0;
            }
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

    if (denominatorSquared < 1e-12){ // avoid division by 0
        return c;
    }

    double sign = (delta >= 0) ? 1.0 : -1.0;
    return  c - (sign * std::pow(b, 2)) / (std::abs(delta) + std::sqrt(denominatorSquared));
}

std::tuple<double, double> givensCommon(double a, double b, double epsilon){

    double radius = std::hypot(a, b);

    if (radius < epsilon || std::abs(b) < epsilon){
        return std::tuple{0.0, 1.0};
    }

    double sin = b / radius;
    double cos = a / radius;

    return std::tuple{sin, cos};
}

void matMul(Eigen::MatrixXd &full, double sin, double cos, int i, bool right){

    if (!right){        
        Eigen::RowVectorXd temp1 = full.row(i) * cos + full.row(i + 1) * sin;
        Eigen::RowVectorXd temp2 = -full.row(i) * sin + full.row(i + 1) * cos;

        full.row(i) = temp1;
        full.row(i + 1) = temp2;
    } else {
        Eigen::VectorXd temp1 = full.col(i) * cos + full.col(i + 1) * sin;
        Eigen::VectorXd temp2 = -full.col(i) * sin + full.col(i + 1) * cos;

        full.col(i) = temp1;
        full.col(i + 1) = temp2;
    }
}

void matMulTranspose(Eigen::MatrixXd &full, double sin, double cos, int i){
    
    Eigen::VectorXd temp1 = full.col(i) * cos + full.col(i + 1) * sin;
    Eigen::VectorXd temp2 = -full.col(i) * sin + full.col(i + 1) * cos;

    full.col(i) = temp1;
    full.col(i + 1) = temp2;
}


void applyShift(Eigen::MatrixXd &B, Eigen::MatrixXd &U, Eigen::MatrixXd &V_transposed, double epsilon){
    assert(B.rows() == B.cols() && "The input matrix must be square");

    std::tuple<double, double> tuple;
    double a = 0.0;
    double b = 0.0;
    bool right = true;
    double miu = wilkinsonShift(chooseSubmatrix(B));

    // numRows() - 1 since the last non-zero element from the last row is the one from the main diagonal 
    for (int i = 0; i < B.rows() - 1; i++){ 
        assert(i >= 0 && i + 1 < B.rows() && i + 1 < B.cols());

        if (i != 0){ 
            a = B(i, i);
            b = B(i, i+1);
        } else { // smart miu
            a = std::pow(B(0, 0), 2) - miu;
            b = B(0, 0) * B(0, 1);
        }
        
        auto [sin, cos] = givensCommon(a, b, epsilon);
        matMul(V_transposed, sin, cos, i, right);
        matMulTranspose(B, sin, cos, i);

        a = B(i, i);
        b = B(i + 1, i);
        
        auto [sin1, cos1] = givensCommon(a, b, epsilon);
        matMul(U, sin1, cos1, i, !right);
        matMul(B, sin1, cos1, i, !right);
    }
}


std::vector<Eigen::MatrixXd> svdGolubKahan(Eigen::MatrixXd B, Eigen::MatrixXd &U, Eigen::MatrixXd &V_transposed){

    double epsilon = 1e-8;
    int iterationNumber = 0;

    while (!isDiagonal(B)){

        double offDiagonalNorm = (B - B.diagonal().asDiagonal().toDenseMatrix()).norm();
        if (iterationNumber % 20 == 0){
            std::cout << " now at iteration " << iterationNumber << std::endl;
            std::cout << "norm = " << offDiagonalNorm << std::endl;        
        }

        applyShift(B, U, V_transposed, epsilon);
        deflateValues(B, epsilon);

        if (iterationNumber >= 1500) {
            break;
        }

        iterationNumber += 1;
    }

    std::cout << "B : " << std::endl;
    std::cout << B.rows() << " x " << B.cols() << std::endl;

    std::cout << "U : " << std::endl;
    std::cout << U.rows() << " x " << U.cols() << std::endl;

    std::cout << "V_transposed : " << std::endl;
    std::cout << V_transposed.rows() << " x " << V_transposed.cols() << std::endl;

    return std::vector{B, U, V_transposed};
}