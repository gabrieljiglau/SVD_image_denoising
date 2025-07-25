#include "lodepng.h"
#include "bidiagonalization.hpp"
#include "imageUtils.hpp"
#include <Eigen/src/Core/util/Constants.h>
#include <Eigen/src/SVD/JacobiSVD.h>
#include <filesystem>
#include <fstream>
#include <vector>
#include <optional>
#include <string>
#include <iostream>
#include <fmt/core.h>
namespace fs = std::filesystem;

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

void matrixToCsv(Eigen::MatrixXd &matrix, fs::path csvOutput, bool overwrite){

    if (!overwrite){
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

Eigen::MatrixXd csvToMatrix(const fs::path &filename){

    std::ifstream file(filename);
    std::string line;
    std::vector<std::vector<double>> values;

    size_t cols = 0;
    while (std::getline(file, line)){

        std::stringstream linestream(line);
        std::string element;
        std::vector<double> row;

        while (std::getline(linestream, element, ',')){
            
            if (!element.empty()){
                row.push_back(std::stod(element));
            }
        }

        if (!row.empty()){
            cols = std::max(cols, row.size());
            values.push_back(row);
        }
    }

    size_t rows = values.size();
    Eigen::MatrixXd matrix(rows, cols);
    for (size_t i = 0; i < rows; i++){
        for (size_t j = 0; j < values[i].size(); j++){
            matrix(i,j) = values[i][j];
        }
    }

    return matrix;
}

int main(){

    std::vector<unsigned char> image;
    unsigned width; // img_row
    unsigned height; // img_col

    std::string myImage = "/home/gabriel/Documents/HolyC/SVD_image_denoising/images/noisy_mu0_std10.png";
    unsigned error = lodepng::decode(image, width, height, myImage);

    if (error) {
        std::cerr << "Error when decoding" << error << " : " 
                  << lodepng_error_text(error) << std::endl;

        return 1;
    }

    std::cout << "Image is " << width << " x " << height << "\n";

    // printImage(image, height, width);
    std::optional<ChannelData> channelR = rgbChannel(image, height, width, 'r');
    
    if (!channelR){
        std::cerr << "Wrong input as the RGB channel" << std::endl;
    }

    // adding noise
    // addNoise(image, height, width);

    //std::vector<Eigen::MatrixXd> matrices = bidiagonalize(channelR->matrixChannel, channelR->channelRows, channelR->channelCols);
    
    //Eigen::MatrixXd U = matrices[0];
    //Eigen::MatrixXd B_mine = matrices[1];
    //Eigen::MatrixXd V_transpose = matrices[2];

    /*
    fs::path imgPath = myImage;
    fs::path csvPathB = "/home/gabriel/Documents/HolyC/SVD_image_denoising/matrices/channelR/B.csv";
    fs::path csvPathU = "/home/gabriel/Documents/HolyC/SVD_image_denoising/matrices/channelR/U.csv";
    fs::path csvPathV = "/home/gabriel/Documents/HolyC/SVD_image_denoising/matrices/channelR/V.csv";
    */

    //matrixToCsv(B_mine, csvPathB, true); 
    //matrixToCsv(U, csvPathU, true); 
    //matrixToCsv(V_transpose, csvPathV, true); 

    Eigen::Map<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> mappedMatrix(channelR->matrixChannel.data(), 
        channelR->channelRows, channelR->channelCols);
    
    Eigen::MatrixXd A = mappedMatrix.cast<double>();
    /*
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(A);
    double cond = svd.singularValues()(0) / svd.singularValues().tail(1)(0);
    std::cout << "Condition number: " << cond << std::endl;
    */

    /*=
    Eigen::JacobiSVD<Eigen::MatrixXd> svdDecomp(A, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Eigen::MatrixXd U = svdDecomp.matrixU();
    Eigen::MatrixXd V = svdDecomp.matrixV();
    Eigen::MatrixXd B = U.transpose() * A * V;

    std::cout << B.rows() << " X " << B.cols() << std::endl;

    Eigen::MatrixXd B_mine = csvToMatrix(csvPathB);

    std::cout << "(B_mine - B).norm = " << (B_mine - B).norm() << std::endl;
    std::cout << B_mine.rows() << " X " << B_mine.cols() << std::endl;
    */
    // the bidiagonalization step must be applied for each channel
    //functie care sa returneze toate canalele bidiagonalizate

    return 0;
}