#pragma once
#include <vector>


float gaussianPdf(int mean, int stdDev);

int addNoise(int originalPixel, int mean, int stdDev);

void printImage(std::vector<unsigned char> imagePixels, int height, int width);
