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
std::vector<Eigen::MatrixXd> bidiagonalize(std::vector<int> rgbMatrix, int numRows, int numCols);

Eigen::MatrixXd leftReflection(Eigen::MatrixXd &inputArr, int numReflection, int maxReflection);

Eigen::MatrixXd rightReflection(Eigen::MatrixXd &inputArr, int numReflection, int maxReflection);

Eigen::MatrixXd resizeH(Eigen::MatrixXd H, const int maxSize, const int currentSize);

std::string checkBidiagonality(Eigen::MatrixXd B, int numRows, int numCols);