#pragma once
#include <vector>
#include <Eigen/Dense>


bool isDiagonal(Eigen::MatrixXd B, double epsilon);

Eigen::MatrixXd chooseSubmatrix(Eigen::MatrixXd B);

double wilkinsonShift(Eigen::MatrixXd B);

void deflateValues(Eigen::MatrixXd &B, double epsilon);

Eigen::MatrixXd applyShift(Eigen::MatrixXd B, Eigen::MatrixXd &U, Eigen::MatrixXd &V_transposed, double miu, double epsilon);

std::vector<Eigen::MatrixXd> golubKahan(Eigen::MatrixXd B, Eigen::MatrixXd &U, Eigen::MatrixXd &V_transposed);