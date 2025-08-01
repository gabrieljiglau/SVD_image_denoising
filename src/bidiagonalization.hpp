#pragma once
#include <Eigen/Dense>


Eigen::MatrixXd leftReflection(Eigen::MatrixXd &inputArr, int numReflection, int maxReflection);

Eigen::MatrixXd rightReflection(Eigen::MatrixXd &inputArr, int numReflection, int maxReflection);

Eigen::MatrixXd resizeH(Eigen::MatrixXd H, const int maxSize, const int currentSize);