#include "matrixUtils.hpp"
#include <random>
#include <cmath>
#include <iostream>

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

int addNoise(int originalPixel, int mean, int stdDev){

    int noise = (int) gaussianGenerator(mean, stdDev);
    return (originalPixel + noise) % 256;

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