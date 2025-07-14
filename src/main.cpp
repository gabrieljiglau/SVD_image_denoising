#include "lodepng.h"
#include "matrixUtils.hpp"
#include "imageUtils.hpp"
#include <vector>
#include <optional>
#include <string>
#include <iostream>
#include <fmt/core.h>



int main(){

    std::vector<unsigned char> image;
    unsigned width; // img_row
    unsigned height; // img_col

    std::string myImage = "/home/gabriel/Documents/HolyC/SVD_image_denoising/images/original.png";
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

    if (channelR->matrixChannel.empty()){
        std::cerr << "The vector is empty" << std::endl;
    } else {
        std::cout << "inside main" << std::endl;
        for (const int pixel : channelR->matrixChannel){
            // std::cout << pixel << std::endl;
        }
        
        std::cout << "Cols = " << channelR->channelCols << std::endl;
        std::cout << "Rows = " << channelR->channelRows << std::endl;
    }



    // adding noise

    /*
    int mean = 0; //mu
    std::vector<int> stdDevs = {10, 20, 30};
    std::string newPath;
    std::vector<unsigned char> newImage;

    for (int stdDev : stdDevs){
        newPath = fmt::format("noisy_mu{mu:.2f}_std{:.2f}.png", 
                             fmt::arg("mu", mean),
                             fmt::arg("std", stdDev));
        generateNoisyImage(image, newImage, newPath, height, width, mean, stdDev);
    }
    */

    return 0;
}