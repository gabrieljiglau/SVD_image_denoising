#pragma once
#include <vector>
#include <Eigen/Dense>


Eigen::MatrixXd chooseSubmatrix(Eigen::MatrixXd B);

double wilkinsonShift(Eigen::MatrixXd B);

void applyShift(Eigen::MatrixXd &B, Eigen::MatrixXd &U, Eigen::MatrixXd &V_transposed, double epsilon);

std::vector<Eigen::MatrixXd> svdGolubKahan(Eigen::MatrixXd B, Eigen::MatrixXd &U, Eigen::MatrixXd &V_transposed);