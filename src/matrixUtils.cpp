#include "matrixUtils.hpp"
#include <Eigen/src/Core/util/Constants.h>
#include <algorithm>
#include <vector>
#include <iostream>
#include <optional>
#include <Eigen/Dense>


std::optional<ChannelData> rgbChannel(std::vector<unsigned char> image, int height, int width, char channel){

    std::vector<int> matrix;
    int channelCols = 0;
    int channelRows = 0;
    bool firstTime = true;

    for (unsigned int col = 0; col < height; col++){
        channelCols += 1;
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
            if (firstTime) {
                channelRows += 1;
            }
        }
        firstTime = false;
    }

    return ChannelData{matrix, channelCols, channelRows};
}


// dynamic matrix of ints
Eigen::MatrixXd bidiagonalize(std::vector<int> rgbMatrix, int numRows, int numCols){
    // the bidiagonalization step must be applied for each channel
    
    Eigen::Map<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> mappedMatrix(rgbMatrix.data(), numRows, numCols);
    Eigen::MatrixXd A = mappedMatrix.cast<double>();


    // se aplica alternativ leftReflection si rightReflection
    // aplici leftReflection de min(rows, cols)
    // aplici rightReflection de min(rows, cols) - 1

    Eigen::VectorXi e1 = Eigen::VectorXi::Unit(numCols, 0); // first column from the identity matrix
    Eigen::MatrixXi I = Eigen::MatrixXi::Identity(numRows, numCols);
    int numReflections = std::min(numRows, numCols);
    int sign = 1;

    for (int reflection = 0; reflection < numReflections; reflection++){

        leftReflection(A, I, e1, reflection);
        rightReflection(A, I, e1, reflection);
    }

    
}

void leftReflection(Eigen::MatrixXd &inputArr, const Eigen::MatrixXi I, const Eigen::VectorXi e1, int numReflection){
    
    Eigen::VectorXd H;
    Eigen::VectorXd col = inputArr.col(numReflection);
    auto slice = col.segment(numReflection, col.size() - numReflection); // the arguments are (start_idx, length)

    double norm = col.norm();
    norm = (norm > 0) ? norm : (-1) * norm;
    Eigen::VectorXd v = col + norm * e1.cast<double>();

    H = I.cast<double>() - 2 * v * v.transpose();
    inputArr *= H;
}

void rightReflection(Eigen::MatrixXd &inputArr, const Eigen::MatrixXi I, const Eigen::VectorXi e1, int numReflection){

    Eigen::VectorXd H;
    Eigen::VectorXd row = inputArr.row(numReflection);
    // numReflection +1 as the starting index, since the element from the main diagonal was already diagonalised 
    auto slice = row.segment(numReflection + 1, row.size() - numReflection);

    double norm = row.norm();
    norm = (norm > 0) ? norm: (-1) * norm;
    Eigen::VectorXd w = row + norm * e1.cast<double>();
    H = I.cast<double>() - 2 * w * w.transpose();
    inputArr *= H;
}