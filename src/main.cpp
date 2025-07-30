#include "lodepng.h"
#include "bidiagonalization.hpp"
#include "imageUtils.hpp"
#include "golubKahan.hpp"
#include "tests.hpp"
#include <Eigen/src/Core/Matrix.h>
#include <Eigen/src/Core/util/Constants.h>
#include <Eigen/src/SVD/JacobiSVD.h>
#include <filesystem>
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

void svdChannels(){

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

    // tests
    //bidiagonalizeTest();
    
    //golubKahanTest();

    Eigen::Map<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> mappedMatrix(channelR->matrixChannel.data(), 
    channelR->channelRows, channelR->channelCols);

    Eigen::MatrixXd R_channel = mappedMatrix.cast<double>();
    //std::vector<Eigen::MatrixXd> matrices = bidiagonalize(A, A.rows(), A.cols);

    fs::path imgPath = myImage;
    fs::path csvPathB = "/home/gabriel/Documents/HolyC/SVD_image_denoising/matrices/channelR/B.csv";
    fs::path csvPathU = "/home/gabriel/Documents/HolyC/SVD_image_denoising/matrices/channelR/U.csv";
    fs::path csvPathV = "/home/gabriel/Documents/HolyC/SVD_image_denoising/matrices/channelR/V.csv";

    fs::path csvFinalB = "/home/gabriel/Documents/HolyC/SVD_image_denoising/matrices/channelR/B_final.csv";
    fs::path csvFinalU = "/home/gabriel/Documents/HolyC/SVD_image_denoising/matrices/channelR/U_final.csv";
    fs::path csvFinalV = "/home/gabriel/Documents/HolyC/SVD_image_denoising/matrices/channelR/V_final.csv";
    //matrixToCsv(B_mine, csvPathB, true); 
    //matrixToCsv(U, csvPathU, true); 
    //matrixToCsv(V_transpose, csvPathV, true); 
    
    Eigen::JacobiSVD<Eigen::MatrixXd> svdDecomp(R_channel, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Eigen::MatrixXd R_svd = svdDecomp.matrixU()
                          * svdDecomp.singularValues().asDiagonal()
                          * svdDecomp.matrixV().transpose();
    
    if (fs::exists(csvFinalU)){
        Eigen::MatrixXd U_final = csvToMatrix(csvFinalU);
        U_final *= 255.0;
    }
    Eigen::MatrixXd U_final = csvToMatrix(csvFinalU);
    U_final *= 255.0;
    Eigen::MatrixXd B_final = csvToMatrix(csvFinalB);
    cleanDiagonalMatrix(B_final);
    B_final *= 255.0;
    Eigen::MatrixXd V_final = csvToMatrix(csvFinalV);
    V_final *= 255.0;

    Eigen::MatrixXd R_recon = U_final * B_final * V_final.transpose();

    double d2NormDiff = (R_svd - R_recon).norm();
    double d2Norm = R_channel.norm();

    std::cout << "Relative construction error " << d2NormDiff / (d2Norm * 65025) << std::endl;
    
    auto mySingularValues = B_final.diagonal();
    auto eigenSingularValues = svdDecomp.singularValues();

    std::cout << mySingularValues.transpose() << std::endl;

    std::cout << std::endl << " S E P A R A T O R " << std::endl;

    std::cout << eigenSingularValues.transpose() << std::endl;

    // the bidiagonalization step must be applied for each channel
    //functie care sa returneze toate canalele bidiagonalizate
    


    return 0;
}