#include "lodepng.h"
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

    for (unsigned int col = 0; col < height; col++){  // images are stored row-by-row
        for (unsigned int row = 0; row < width; row++){
            
            unsigned int idx = 4 * (col * width + row); // RGBA format

            unsigned char r = image[idx];
            unsigned char g = image[idx + 1];
            unsigned char b = image[idx + 2];
            unsigned char a = image[idx + 3];
            
            std::cout << "Pixel (" << row << ", " << col << "): "
                      << "R = " << (int)r << ", G = " << (int)g 
                      << ", B = "<< (int)b << ", A = " << (int)a << std::endl;
        }
    }

    return 0;
}