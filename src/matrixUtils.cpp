#include "matrixUtils.hpp"
#include <vector>
#include <iostream>
#include <optional>
#include <Eigen/Dense>


std::optional<std::vector<int>> rgbChannel(std::vector<unsigned char> image, int height, int width, char channel){

    std::vector<int> matrix;

    for (unsigned int col = 0; col < height; col++){
        for (unsigned int row = 0; row < width; row++){
            int idx = 4 * (col * width + row);
            
            int pixel = 0;
            switch (channel) {
                case 'r':
                    pixel = (int) image[idx];
                    break;
                case 'g':
                    pixel = (int) image[idx + 1];
                    break;
                case 'b':
                    pixel = (int) image[idx + 2];
                    break;
                default:
                    std::cout << "The only correct channels are 'r', 'g' and 'b'" << std::endl;
                    return std::nullopt;
            }

            matrix.push_back(pixel);
        }

    }

    return matrix;
}

std::vector<int> bidiagonalize();

std::vector<int> leftReflection();

std::vector<int> rightReflection()l