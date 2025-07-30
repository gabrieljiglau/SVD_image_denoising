# pragma once
#include <filesystem>
#include <fmt/core.h>
#include <vector>
#include <string>
#include <Eigen/Dense>
namespace fs = std::filesystem;


int addPixelNoise(int originalPixel, int mean, int stdDev);

void generateNoisyImage(std::vector<unsigned char> originalImage, std::vector<unsigned char> &newImage, 
    std::string newPath, int height, int width, int mean, int stdDev);

void printImage(std::vector<unsigned char> imagePixels, int height, int width);

void matrixToCsv(Eigen::MatrixXd &matrix, fs::path csvOutput, bool overwrite);

void cleanDiagonalMatrix(Eigen::MatrixXd &B);

Eigen::MatrixXd csvToMatrix(const fs::path &filename);