# pragma once
#include <filesystem>
#include <fmt/core.h>
#include <vector>
#include <string>
#include <Eigen/Dense>
namespace fs = std::filesystem;

struct ChannelData{
    std::vector<int> matrixChannel;
    int channelCols;
    int channelRows;
};

std::string checkBidiagonality(Eigen::MatrixXd B, int numRows, int numCols);

std::vector<ChannelData> rgbChannel(std::vector<unsigned char> image, int height, int width);

// dynamic matrix of ints
std::vector<Eigen::MatrixXd> bidiagonalize(Eigen::MatrixXd A, int numRows, int numCols);

bool isDiagonal(Eigen::MatrixXd B);

Eigen::MatrixXd deflateValues(Eigen::MatrixXd B, double epsilon);

void reconstructImage(std::vector<Eigen::MatrixXd> truncatedChannels, std::vector<unsigned char> &newImage,
    std::string newPath, int height, int width, int idx);

void addNoise(std::vector<unsigned char>image, int height, int width);

void printImage(std::vector<unsigned char> imagePixels, int height, int width);

void matrixToCsv(Eigen::MatrixXd &matrix, fs::path csvOutput, bool overwrite);

void cleanDiagonalMatrix(Eigen::MatrixXd &B);

Eigen::MatrixXd truncateU(Eigen::MatrixXd &U, int k);

Eigen::MatrixXd truncateB(Eigen::MatrixXd &B, int k);

Eigen::MatrixXd truncateV(Eigen::MatrixXd &V, int k);

Eigen::MatrixXd csvToMatrix(const fs::path &filename);