#include "bidiagonalization.hpp"
#include <Eigen/src/Core/Matrix.h>
#include <Eigen/src/Core/util/Constants.h>
#include <algorithm>
#include <fmt/base.h>
#include <vector>
#include <iostream>
#include <optional>
#include <Eigen/Dense>


std::string checkBidiagonality(Eigen::MatrixXd B, int numRows, int numCols){
    bool is_bidiagonal = true;

    for (int i = 0; i < numRows; i++){
        for (int j = 0; j < numCols; j++){
            if ( (i == j) || (j == i + 1)){
                continue;
            } else {
                if ((std::abs(B(i, j)) > 1e-10)){
                    is_bidiagonal = false;
                    break;
                }
            }
        }
    }

    if (is_bidiagonal == 1) {
        return " It is bidiagonal !!";
    }
    
    return" It is NOT bidiagonal !!";
}

std::optional<ChannelData> rgbChannel(std::vector<unsigned char> image, int height, int width, char channel){

    std::vector<int> matrix;
    int channelCols = 0;
    int channelRows = 0;
    bool firstTime = true;

    for (unsigned int col = 0; col < height; col++){ // height
        channelCols += 1;
        for (unsigned int row = 0; row < width; row++){ //width
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
std::vector<Eigen::MatrixXd> bidiagonalize(Eigen::MatrixXd A, int numRows, int numCols){
    
    Eigen::MatrixXd A_balanced = A / 255.0;
    Eigen::MatrixXd A_copy = A_balanced;


    int numReflections = std::min(numRows, numCols) - 1 ;
    Eigen::MatrixXd U = Eigen::MatrixXd::Identity(numRows, numRows);
    Eigen::MatrixXd V = Eigen::MatrixXd::Identity(numCols, numCols);

    // apply alternately leftReflection and rightReflection
    for (int reflection = 0; reflection < numReflections; ++reflection){

        std::cout << "Now at iteration " << reflection << std::endl;
        
        U = leftReflection(A_balanced, reflection, numRows) * U ;

        if (reflection + 1 != numReflections){ // 1 less rotations possible, since the left rotations reduces the problem by 1
            V = V * rightReflection(A_balanced, reflection, numCols);
        }
        
    }

    std::vector<Eigen::MatrixXd> matricesList;

    // sanity checks
    std::cout << "V_T * V - I= " << (((V.transpose() * V) - Eigen::MatrixXd::Identity(V.rows(), V.cols())).norm())<< std::endl;
    std::cout << "U_T * U -I = " << (((U.transpose() * U) - Eigen::MatrixXd::Identity(U.rows(), U.cols())).norm())<< std::endl;

    // more sanity checks, it's impossible to perfectly reconstruct the original matrix ??? 
    Eigen::MatrixXd B_explicit = U.transpose() * A_copy * V;

    matricesList.push_back(clampSmallValues(U));
    matricesList.push_back(clampSmallValues(A_balanced));
    matricesList.push_back(clampSmallValues(V).transpose());
    
    auto recon1 = U * A_balanced * V.transpose();    
    auto recon2 = U * B_explicit * V.transpose();
    
    std::cout << "||A_copy - U * B_bidiagonal * V^T|| = " << (A_copy - recon1).norm() << std::endl; 
    std::cout << "||A_copy - U * B_explicit * V^T|| = " << (A_copy - recon2).norm() << std::endl;   

    std::cout << "IS B Bidiagonal ?? " << checkBidiagonality(A_balanced, numRows, numCols) << std::endl;
    std::cout << "IS B_explicit Bidiagonal ?? " << checkBidiagonality(B_explicit, numRows, numCols) << std::endl;

    return matricesList;
}

Eigen::MatrixXd leftReflection(Eigen::MatrixXd &inputArr, int numReflection, int maxReflection){

    int rowNum = inputArr.rows();
    int colNum = inputArr.cols();
    
    Eigen::VectorXd col = inputArr.col(numReflection);

    auto slice = col.segment(numReflection, col.size() - numReflection); // the arguments are (start_idx, length)
    int k = slice.size();
    
    Eigen::MatrixXd I = Eigen::MatrixXd::Identity(k, k);
    Eigen::VectorXd e1 = I.col(0); // first column from the identity matrix

    double alpha = slice[0];
    double norm = slice.norm();
    double sign = (alpha >= 0) ? 1.0 : -1.0;
    
    Eigen::VectorXd v = slice;
    v[0] += sign * norm;

    double normalizationConstant = 2.0 / v.squaredNorm();
    Eigen::MatrixXd H = I - normalizationConstant * v * v.transpose();

    Eigen::Block<Eigen::MatrixXd> trail = inputArr.block(numReflection, numReflection, rowNum - numReflection, colNum - numReflection);
    trail = H * trail;

    return resizeH(H, maxReflection, numReflection);
}

Eigen::MatrixXd rightReflection(Eigen::MatrixXd &inputArr, int numReflection, int maxReflection){

    int rowNum = inputArr.rows();
    int colNum = inputArr.cols();

    Eigen::MatrixXd H;
    Eigen::VectorXd row = inputArr.row(numReflection);

    int remaining = colNum - (numReflection + 1);
    if (remaining <= 0) {
        return Eigen::MatrixXd::Identity(maxReflection, maxReflection);
    }

    // numReflection + 1 as the starting index, 
    // since the element from the main diagonal was already diagonalised by the left reflection
    auto slice = row.segment(numReflection + 1, row.size() - (numReflection + 1));
    int k = slice.size();
    Eigen::MatrixXd I = Eigen::MatrixXd::Identity(k, k);

    double alpha = slice[0];
    double norm = slice.norm();
    double sign = (alpha >= 0) ? 1.0 : -1.0;
    
    if (norm < 1e-10){
        return I;
    } 

    Eigen::VectorXd w = slice;
    w[0] += sign * norm;
    double normalizationConstant = 2.0 / w.squaredNorm();
    H = I - normalizationConstant * w * w.transpose();

    // do not update the matrix in place !
    Eigen::Block<Eigen::MatrixXd> trail = inputArr.block(numReflection, numReflection + 1, rowNum - numReflection, colNum - (numReflection + 1));
    trail = trail * H;

    return resizeH(H, maxReflection, numReflection);
}