#include "lodepng.h"
#include "matrixUtils.hpp"
#include <vector>
#include <string>
#include <iostream>



int main(){

    std::vector<unsigned char> image;
    unsigned width; // img_row
    unsigned height; // img_col

    std::string myImage = "../images/original.png";
    unsigned error = lodepng::decode(image, width, height, myImage);

    if (error) {
        std::cerr << "Error when decoding" << error << " : " 
                  << lodepng_error_text(error) << std::endl;

        return 1;
    }

    std::cout << "Image is " << width << " x " << height << "\n";

    std::vector<unsigned char> newImage;


    // adding noise
    /*
    for (unsigned int col = 0; col < height; col++){
        for (unsigned int row = 0; row < width; row++){
            int idx = 4 * (col * width + row);

            int r = (int) image[idx];
            int g = (int) image[idx + 1];
            int b = (int) image[idx + 2];
            int a = (int) image[idx + 3];

            r += addNoise(r, 0, 30);
            g += addNoise(g, 0, 30);
            b += addNoise(b, 0, 30);

            newImage.push_back(r);
            newImage.push_back(g);
            newImage.push_back(b);
            newImage.push_back(a);
        }
    }    

    error = lodepng::encode("../images/noisy_mu0_std30.png", newImage, width, height);
    if (error) {
        std::cerr << "Encoder error: " << error << lodepng_error_text(error) << std::endl;
    }
    */

    return 0;
}