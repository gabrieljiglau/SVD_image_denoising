#include "matrixUtils.hpp"
#include <Eigen/src/Core/Matrix.h>
#include <Eigen/src/Core/util/Constants.h>
#include <algorithm>
#include <fmt/base.h>
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

Eigen::MatrixXd clampSmallValues(Eigen::MatrixXd &A){

    for (Eigen::Index i = 0; i < A.rows(); i++){
        for(Eigen::Index j = 0; j < A.cols(); j++){
            if (!(std::abs(A(i, j)) > 1e-10)) {
                A(i,j) = 0.0;
            }
        }
    }

    return A;
}

Eigen::MatrixXd resizeH(Eigen::MatrixXd H, const int fullSize, const int reflectionIndex){

    Eigen::MatrixXd I = Eigen::MatrixXd::Identity(fullSize, fullSize);

    int k = H.rows();
    I.block(reflectionIndex, reflectionIndex, k, k) = H;

    return I;
}

// dynamic matrix of ints
std::vector<Eigen::MatrixXd> bidiagonalize(std::vector<int> rgbMatrix, int numRows, int numCols){
    
    Eigen::Map<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> mappedMatrix(rgbMatrix.data(), numRows, numCols);
    Eigen::MatrixXd A = mappedMatrix.cast<double>();

    int numReflections = std::min(numRows, numCols) - 1; // 1 less rotations possible than the min of rows, cols
    Eigen::MatrixXd U = Eigen::MatrixXd::Identity(numRows, numRows);
    Eigen::MatrixXd V = Eigen::MatrixXd::Identity(numCols, numCols);

    // apply alternately leftReflection and rightReflection
    for (int reflection = 0; reflection < numReflections; reflection++){

        std::cout << "Now at iteration " << reflection << std::endl;
        U *= leftReflection(A, reflection, numRows);
        // leftReflection(A, reflection, numReflections + 1);
        V *= rightReflection(A, reflection, numCols);
    }

    std::vector<Eigen::MatrixXd> matricesList;
    matricesList.push_back(clampSmallValues(U));
    matricesList.push_back(clampSmallValues(A));
    std::cout << "V^T (last row):\n" << V.transpose().row(V.cols() - 1) << "\n";
    matricesList.push_back(V.transpose());

    Eigen::MatrixXd V_transpose = V.transpose();
    Eigen::MatrixXd V_mult = V * V_transpose.transpose();
    std::cout << "(v_mult - v).norm = " << (V_mult - Eigen::MatrixXd::Identity(V_transpose.rows(), V_transpose.cols())).norm() << std::endl; 


    return matricesList;
}

Eigen::MatrixXd leftReflection(Eigen::MatrixXd &inputArr, int numReflection, int maxReflection){

    int rowNum = inputArr.rows();
    int colNum = inputArr.cols();
    
    Eigen::MatrixXd H;
    Eigen::VectorXd col = inputArr.col(numReflection);

    auto slice = col.segment(numReflection, col.size() - numReflection); // the arguments are (start_idx, length)
    int k = slice.size();
    Eigen::MatrixXi I = Eigen::MatrixXi::Identity(k, k);
    Eigen::MatrixXd U = I.cast<double>();
    Eigen::VectorXi e1 = I.col(0); // first column from the identity matrix

    double norm = slice.norm();
    norm = (norm > 0) ? norm : (-1) * norm;
    
    Eigen::VectorXd v = slice + norm * e1.cast<double>();
    double normalizationConstant = 2.0 / v.squaredNorm();
    H = I.cast<double>() - normalizationConstant * v * v.transpose();

    int rowsLeft = rowNum - numReflection;
    int colsLeft = colNum - numReflection;
    inputArr.block(numReflection, numReflection, rowsLeft, colsLeft) = 
        H * inputArr.block(numReflection, numReflection, rowsLeft, colsLeft);

    return resizeH(H, maxReflection, numReflection);
}

Eigen::MatrixXd rightReflection(Eigen::MatrixXd &inputArr, int numReflection, int maxReflection){

    int rowNum = inputArr.rows();
    int colNum = inputArr.cols();

    Eigen::MatrixXd H;
    Eigen::VectorXd row = inputArr.row(numReflection);

    int remaining = colNum - (numReflection + 1);
    if (remaining == 1) {
        return Eigen::MatrixXd::Identity(maxReflection, maxReflection);
    }

    // numReflection + 1 as the starting index, 
    // since the element from the main diagonal was already diagonalised by the left reflection
    auto slice = row.segment(numReflection + 1, row.size() - (numReflection + 1));
    int k = slice.size();
    Eigen::MatrixXi I = Eigen::MatrixXi::Identity(k, k);
    Eigen::VectorXi e1 = I.col(0);

    double norm = slice.norm();
    norm = (norm > 0) ? norm: (-1) * norm;

    Eigen::VectorXd w = slice + norm * e1.cast<double>(); //e1 should be the same size as the slice
    double normalizationConstant = 2.0 / w.squaredNorm();
    H = I.cast<double>() - normalizationConstant * w * w.transpose();

    int rowsLeft = rowNum - numReflection;
    int colsLeft = colNum - (numReflection + 1);

    // numReflection + 1 as the starting index, since the element from the main diagonal was already diagonalised 
    inputArr.block(numReflection, numReflection + 1, rowsLeft, colsLeft) = 
        inputArr.block(numReflection, numReflection + 1, rowsLeft, colsLeft) * H;

    return resizeH(H, maxReflection, numReflection);
}