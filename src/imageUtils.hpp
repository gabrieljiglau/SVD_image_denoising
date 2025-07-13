#pragma once
#include <vector>
#include <string>


int addPixelNoise(int originalPixel, int mean, int stdDev);

void generateNoisyImage(std::vector<unsigned char> originalImage, std::vector<unsigned char> &newImage, 
    std::string newPath, int height, int width, int mean, int stdDev);

void printImage(std::vector<unsigned char> imagePixels, int height, int width);