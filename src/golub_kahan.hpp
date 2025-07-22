#pragma once
#include <vector>
#include <Eigen/Dense>


bool isDiagonal(Eigen::MatrixXd B);

Eigen::MatrixXd chooseSubmatrix(Eigen::MatrixXd B);

double wilkinsonShift(Eigen::MatrixXd B);

Eigen::MatrixXd applyShift(Eigen::MatrixXd trail, double miu);

Eigen::MatrixXd golubKahan(Eigen::MatrixXd U, Eigen::MatrixXd &B, Eigen::MatrixXd V_transposed);