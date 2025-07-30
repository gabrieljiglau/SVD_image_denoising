#include "imageUtils.hpp"
#include "lodepng.h"
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

void cleanDiagonalMatrix(Eigen::MatrixXd &B){

    for (int i = 0; i < B.rows(); i++){
        for (int j = 0; j < B.cols(); j++){
            if (i != j){
                B(i, j) = 0;
            }
        }
    }
}

void matrixToCsv(Eigen::MatrixXd &matrix, fs::path csvOutput, bool overwrite){

    if (!overwrite and !fs::exists(csvOutput)){
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
                    row.push_back(std::stod(element));
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

void printImage(std::vector<unsigned char> image, int height, int width){
    
    for (unsigned int col = 0; col < height; col++){
        for (unsigned int row = 0; row < width; row++){
            int idx = 4 * (col * width + row);

            std::cout << "(" << (int) image[idx] << ", "
                      << (int) image[idx + 1] << ", "
                      << (int) image[idx + 2] << ", "
                      << ")" << std::endl;
        }
    }    

    printf("\n");
}
