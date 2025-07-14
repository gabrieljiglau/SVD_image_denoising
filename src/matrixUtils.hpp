#pragma once
#include <vector>
#include <optional>
#include <Eigen/Dense>

struct ChannelData{
    std::vector<int> matrixChannel;
    int channelCols;
    int channelRows;
};

std::optional<ChannelData> rgbChannel(std::vector<unsigned char> image, int height, int width, char channel);

// dynamic matrix of ints
Eigen::MatrixXd bidiagonalize(std::vector<int> rgbMatrix);

void leftReflection(Eigen::MatrixXd &inputArr, const Eigen::MatrixXi I, const Eigen::VectorXi e1, int numReflection);

void rightReflection(Eigen::MatrixXd inputArr, const Eigen::MatrixXi I, const Eigen::VectorXi e1, int numReflection);