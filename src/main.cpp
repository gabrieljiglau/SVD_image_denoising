#include "lodepng.h"
#include "matrixUtils.hpp"
#include "imageUtils.hpp"
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

void matrixToCsv(Eigen::MatrixXd &matrix, fs::path imagePath, fs::path csvOutput){
    
    if (!fs::exists(imagePath)) {
        std::cerr << "Could not open image: " << imagePath << std::endl;
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
    std::cout << "Successfully written matrix to " << csvOutput << std::endl;
}

// ar mai fi nevoie de o functie csvToMatrix

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

    Eigen::MatrixXd B = bidiagonalize(channelR->matrixChannel, channelR->channelRows, channelR->channelCols);
    fs::path imgPath = myImage;
    fs::path csvPath = "/home/gabriel/Documents/HolyC/SVD_image_denoising/matrices/channel_r.csv";
    matrixToCsv(B, imgPath, csvPath); 

    //functie care sa returneze toate canalele bidiagonalizate

    return 0;
}