#include "lodepng.h"
#include "bidiagonalization.hpp"
#include "imageUtils.hpp"
#include "golub_kahan.hpp"
#include <Eigen/src/Core/Matrix.h>
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

    Eigen::Map<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> mappedMatrix(channelR->matrixChannel.data(), 
    channelR->channelRows, channelR->channelCols);

    Eigen::MatrixXd A = mappedMatrix.cast<double>();
    //std::vector<Eigen::MatrixXd> matrices = bidiagonalize(A, A.rows(), A.cols);

    Eigen::MatrixXd B_test(5, 5);
    B_test << 0.3, 0.7, 0.8, 0.1, 0.3,
        0.25, 0.8, 0.8, 0.5, 0.6,
        0.35 ,0.9 , 0.738, 0.87, 0.33,
        0.45 ,0.6, 0, 0.3634, 0.9,
        0.55, 0.5, 0.19, 0.22, 0.9;

    std::vector<Eigen::MatrixXd> matrices = bidiagonalize(B_test, 5, 5);

    Eigen::MatrixXd U = matrices[0];
    Eigen::MatrixXd B_mine = matrices[1];
    Eigen::MatrixXd V_transpose = matrices[2];
    std::cout << "U = " << std::endl << U << std::endl;
    std::cout << "B = " << std::endl << B_mine << std::endl;
    std::cout << "V_transpose = " << std::endl << V_transpose << std::endl;

    /*
    fs::path imgPath = myImage;
    fs::path csvPathB = "/home/gabriel/Documents/HolyC/SVD_image_denoising/matrices/channelR/B.csv";
    fs::path csvPathU = "/home/gabriel/Documents/HolyC/SVD_image_denoising/matrices/channelR/U.csv";
    fs::path csvPathV = "/home/gabriel/Documents/HolyC/SVD_image_denoising/matrices/channelR/V.csv";

    fs::path csvFinalB = "/home/gabriel/Documents/HolyC/SVD_image_denoising/matrices/channelR/B_final.csv";
    fs::path csvFinalU = "/home/gabriel/Documents/HolyC/SVD_image_denoising/matrices/channelR/U_final.csv";
    fs::path csvFinalV = "/home/gabriel/Documents/HolyC/SVD_image_denoising/matrices/channelR/V_final.csv";
    */
    //matrixToCsv(B_mine, csvPathB, true); 
    //matrixToCsv(U, csvPathU, true); 
    //matrixToCsv(V_transpose, csvPathV, true); 

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

    /*
    Eigen::MatrixXd B = csvToMatrix(csvPathB);
    Eigen::MatrixXd U = csvToMatrix(csvPathU);
    Eigen::MatrixXd V_transposed = csvToMatrix(csvPathV);

    Eigen::MatrixXd B_test(5, 5);
    B_test << -8.02055, 129.93, 0, 0, 0,
        0, -83.5067, 8.72759, 0, 0,
        0 ,0 ,7.01538, 8.72961, 0,
        0 ,0, 0, -12.3634, 9.37061,
        0, 0, 0, 0, -9.53761;

    Eigen::MatrixXd U_test (5, 5);
    U_test << -0.089476, 0.0258819, 0.075377, -0.0699437, 0.0618868,
    -0.0860535,0.0208654,0.078938,-0.0730571,0.0582035,
    -0.0958323,0.0333659,0.0655085,-0.0646888,0.0517716,
    -0.0968101,0.0348304,0.0577696,-0.0538016,0.0463264,
    -0.0904539,0.0252435,0.0793154,-0.0668693,0.0581317;

    Eigen::MatrixXd V_test(5, 5);
    V_test << 0.0556736, 0.0534652, 0.0549524, 0.053282, 0.0464544,
    -0.374102, -0.25997, -0.262998, -0.177746, -0.0951319,
    0.0906982, -0.0372941, -0.00531366, -0.0473925, 0.041103,
    -0.103568, 0.103578, 0.0771154, 0.11051, 0.0401614,
    0.0611464, -0.0704854, -0.0227363, -0.0229629, -0.0318511;

    
    std::vector<Eigen::MatrixXd> newValues = golubKahan(B_test, U_test, V_test);
    matrixToCsv(newValues[0], csvFinalB, true);

    /*
    std::vector<Eigen::MatrixXd> finalMatrices = golubKahan(B, U, V_transposed);
    matrixToCsv(finalMatrices[0], csvFinalB, true);
    matrixToCsv(finalMatrices[1], csvFinalU, true);
    matrixToCsv(finalMatrices[2], csvFinalV, true);

    */

    return 0;
}