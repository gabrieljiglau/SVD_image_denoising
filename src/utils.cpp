#include "utils.hpp"
#include "lodepng.h"
#include <Eigen/src/Core/Matrix.h>
#include <random>
#include <cmath>
#include <iostream>
#include <filesystem>
#include <fmt/core.h>
#include <vector>
#include <fstream>
#include <string>
#include <Eigen/Dense>
namespace fs = std::filesystem;

float gaussianGenerator(int mean, int stdDev){

    constexpr float PI = 3.14159265358979323846f;

    std::random_device rd;
    std::mt19937 generator(rd());

    // uniform distribution
    std::uniform_real_distribution<float> uniformDistribution(0.0f, 1.0f);

    float u1 = uniformDistribution(generator);
    float u2 = uniformDistribution(generator);

    // Box-Muller generator fod Gaussian pdf
    float z0 = std::sqrt(-2 * std::log(u1)) * std::cos(2 * PI * u2);

    return mean + stdDev * z0;
}

int addPixelNoise(int originalPixel, int mean, int stdDev){

    int noise = (int) gaussianGenerator(mean, stdDev);
    return (originalPixel + noise) % 256;

}

Eigen::MatrixXd truncateU(Eigen::MatrixXd &U, int k){

    int rowNum = U.rows();
    return U.block(0, 0, rowNum, k);
}

Eigen::MatrixXd truncateB(Eigen::MatrixXd &B, int k){

    return B.block(0, 0, k, k);
}

Eigen::MatrixXd truncateV(Eigen::MatrixXd &V, int k){

    int colNum = V.cols();
    return V.block(0, 0, k, colNum);
}



void reconstructImage(std::vector<Eigen::MatrixXd> truncatedChannels, std::vector<unsigned char> &newImage,
    std::string newPath, int height, int width, int id){
    
    std::vector<Eigen::RowVectorXd> rowVectors;

    for (Eigen::MatrixXd channel: truncatedChannels){
        Eigen::RowVectorXd innerVector = Eigen::Map<Eigen::RowVectorXd>(channel.data(), channel.size());
        rowVectors.push_back(innerVector);
    }

    int size = rowVectors[0].size();
    for (int i = 0; i < size; i++){

        int r = (int) rowVectors[0](i);
        int g = (int) rowVectors[1](i);
        int b = (int) rowVectors[2](i);
        int a = 255; // I already know a channel is 255 (the image isn't transparent)

        newImage.push_back(r);
        newImage.push_back(g);
        newImage.push_back(b);
        newImage.push_back(a);
    }

    int error = lodepng::encode(newPath, newImage, width, height);
    if (error) {
        std::cerr << "Encoder error: " << error << lodepng_error_text(error) << std::endl;
    }
}

void generateNoisyImage(std::vector<unsigned char> originalImage, std::vector<unsigned char> &newImage, 
    std::string newPath, int height, int width, int mean, int stdDev){

    for (unsigned int col = 0; col < height; col++){
        for (unsigned int row = 0; row < width; row++){
            int idx = 4 * (col * width + row);

            int r = (int) originalImage[idx];
            int g = (int) originalImage[idx + 1];
            int b = (int) originalImage[idx + 2];
            int a = (int) originalImage[idx + 3];

            r += addPixelNoise(r, mean, stdDev);
            g += addPixelNoise(g, mean, stdDev);
            b += addPixelNoise(b, mean, stdDev);

            newImage.push_back(r);
            newImage.push_back(g);
            newImage.push_back(b);
            newImage.push_back(a);
        }
    }    

    int error = lodepng::encode(newPath, newImage, width, height);
    if (error) {
        std::cerr << "Encoder error: " << error << lodepng_error_text(error) << std::endl;
    }

}

void addNoise(std::vector<unsigned char>image, int height, int width){

    int mean = 0; //mu
    std::vector<int> stdDevs = {10, 20, 30};
    std::string newPath;
    std::vector<unsigned char> newImage;

    for (int stdDev : stdDevs){
        newPath = fmt::format("noisy_mu{mu:.2f}_std{:.2f}.png", mean, stdDev);
        generateNoisyImage(image, newImage, newPath, height, width, mean, stdDev);
    }
}


void cleanDiagonalMatrix(Eigen::MatrixXd &B){

    for (int i = 0; i < B.rows(); i++){
        for (int j = 0; j < B.cols(); j++){
            if (i != j){
                B(i, j) = 0;
            }
        }
    }
}

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

Eigen::MatrixXd deflateValues(Eigen::MatrixXd B, double epsilon){

    int numRows = B.rows();
    int numCols = B.cols();

    for (int i = 0; i < numRows; i++){
        for (int j = 0; j < numCols; j++){

            if (std::abs(B(i, j)) < epsilon){
                B(i, j) = 0;
            }
        }
    }

    return B;
}

std::string checkBidiagonality(Eigen::MatrixXd B, int numRows, int numCols){
    bool is_bidiagonal = true;

    for (int i = 0; i < numRows; i++){
        for (int j = 0; j < numCols; j++){
            if ( (i == j) || (j == i + 1)){
                continue;
            } else {
                if ((std::abs(B(i, j)) > 1e-10)){
                    is_bidiagonal = false;
                    break;
                }
            }
        }
    }

    if (is_bidiagonal == 1) {
        return " It is bidiagonal !!";
    }
    
    return" It is NOT bidiagonal !!";
}

void printImage(std::vector<unsigned char> image, int height, int width){
    
    for (unsigned int col = 0; col < height; col++){
        for (unsigned int row = 0; row < width; row++){
            int idx = 4 * (col * width + row);

            std::cout << "(" << (int) image[idx] << ", "
                      << (int) image[idx + 1] << ", "
                      << (int) image[idx + 2] << ", "
                      << (int) image[idx + 3] << ")" << std::endl;
        }
    }    

    printf("\n");
}

std::vector<ChannelData> rgbChannel(std::vector<unsigned char> image, int height, int width){

    std::vector<std::vector<int>> channelMatrices(3);
    int channelCols = 0;
    int channelRows = 0;
    bool firstTime = true;

    for (unsigned int col = 0; col < height; col++){ // height
        channelCols += 1;
        for (unsigned int row = 0; row < width; row++){ //width
            int idx = 4 * (col * width + row);
            
            channelMatrices[0].push_back(image[idx]);
            channelMatrices[1].push_back(image[idx + 1]);
            channelMatrices[2].push_back(image[idx + 2]);

            if (firstTime) {
                channelRows += 1;
            }
        }
        firstTime = false;
    }

    ChannelData R_channel = ChannelData{channelMatrices[0], channelCols, channelRows};
    ChannelData G_channel = ChannelData{channelMatrices[1], channelCols, channelRows};
    ChannelData B_channel = ChannelData{channelMatrices[2], channelCols, channelRows};
    
    return std::vector{R_channel, G_channel, B_channel};
}

void matrixToCsv(Eigen::MatrixXd &matrix, fs::path csvOutput, bool overwrite){

    if (!overwrite and fs::exists(csvOutput)){
        return;
    }

    std::ofstream csv(csvOutput);
    if (!csv) {
        std::cerr << "Could not open csv " << csvOutput << " for writing" << std::endl;
    }

    for (Eigen::Index i = 0; i < matrix.rows(); i++){
        for (Eigen::Index j = 0; j < matrix.cols(); j++){
            csv << matrix(i, j);
            csv << ",";
            if (j + 1 > matrix.cols()){
                csv << ";";
            }
        }
        csv << "\n";
    }
    csv.close();
    // std::cout << "Successfully written matrix to " << csvOutput << std::endl;
}

Eigen::MatrixXd csvToMatrix(const fs::path &filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file: " + filename.string());
    }

    std::string line;
    std::vector<std::vector<double>> values;
    size_t cols = 0;

    while (std::getline(file, line)) {
        std::stringstream linestream(line);
        std::string element;
        std::vector<double> row;

        while (std::getline(linestream, element, ',')) {
            if (!element.empty()) {
                try {
                    row.push_back(std::stod(element) * 255.0);
                } catch (...) {
                    row.push_back(0.0); // if stod fails
                }
            } else {
                row.push_back(0.0); // empty cell -> 0
            }
        }

        cols = std::max(cols, row.size());
        values.push_back(row);
    }

    size_t rows = values.size();
    Eigen::MatrixXd matrix = Eigen::MatrixXd::Zero(rows, cols);

    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < values[i].size(); j++) {
            matrix(i, j) = values[i][j];
        }
    }

    return matrix;
}
