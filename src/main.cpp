#include "lodepng.h"
#include "bidiagonalization.hpp"
#include "utils.hpp"
#include "golubKahan.hpp"
#include "tests.hpp"
#include <Eigen/src/Core/Matrix.h>
#include <Eigen/src/Core/util/Constants.h>
#include <Eigen/src/SVD/JacobiSVD.h>
#include <filesystem>
#include <vector>
#include <optional>
#include <string>
#include <iostream>
#include <fmt/core.h>
namespace fs = std::filesystem;

std::vector<Eigen::MatrixXd> svdChannel(Eigen::MatrixXd &channelMatrix, std::vector<fs::path> rgbChannel,int index){
    

    // channelMatrices -> holds the input RGB matrices
    // outChannel -> holds the final matrices

    std::vector<Eigen::MatrixXd> outChannel;
    if (!fs::exists(rgbChannel[0]) || !fs::exists(rgbChannel[1]) || !fs::exists(rgbChannel[2])){

        for (fs::path &currentPath : rgbChannel){
            fs::create_directories(currentPath.parent_path());
        }
        
        std::cout << "currentPath[0]: " << rgbChannel[0] << " Bidiagonalizing" << std::endl;
        std::vector<Eigen::MatrixXd> bidiagonalOut = bidiagonalize(channelMatrix, channelMatrix.rows(), channelMatrix.cols());
        std::vector<Eigen::MatrixXd> bidiagonalMatrices = bidiagonalOut;

        std::cout << "Now applying golub kahan" << std::endl;
        outChannel = svdGolubKahan(bidiagonalOut[0], bidiagonalOut[1], bidiagonalOut[2]);

        cleanDiagonalMatrix(outChannel[0]);
        // writing the final matrices
        for (int k = 0; k < outChannel.size(); k++){
            outChannel[k] *= 255.0;
            matrixToCsv(outChannel[k], rgbChannel[k], true);
        }
        
        std::cout << "Now saving the decomposition" << std::endl;
    } else {
        for (int k = 0; k < rgbChannel.size(); k++){
            outChannel.push_back(csvToMatrix(rgbChannel[k]));
        }

        std::cout << "currentPath[0]: " << rgbChannel[0] <<" Loading the decomposition" << std::endl;
    }

    return outChannel;
}

std::optional<std::vector<std::vector<Eigen::MatrixXd>>> svdChannels(std::string inputImage, unsigned &width, unsigned &height, std::vector<Eigen::MatrixXd> &channelMatrices, std::vector<std::vector<Eigen::MatrixXd>> &outMatrices,
                std::vector<std::vector<fs::path>> rgbPaths){

    std::vector<unsigned char> image;
    unsigned error = lodepng::decode(image, width, height, inputImage);
    std::vector<std::vector<Eigen::MatrixXd>> svdDecomp(3, std::vector<Eigen::MatrixXd>(3));
                
    // printImage(image, height, width);

    if (error) {
        std::cerr << "Error when decoding" << error << " : " << lodepng_error_text(error) << std::endl;
        return std::nullopt;
    }

    std::vector<ChannelData> channels = rgbChannel(image, height, width);
    
    for (int i = 0; i < channels.size(); i++){
        ChannelData channel = channels[i];
        Eigen::Map<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> mappedMatrix(channel.matrixChannel.data(), channel.channelRows, channel.channelCols);
        channelMatrices[i] = mappedMatrix.cast<double>();
    }

    for (int i = 0; i < channelMatrices.size(); i++){
        svdDecomp[i] = svdChannel(channelMatrices[i], rgbPaths[i], i);
    } 

    return svdDecomp;
}

// B, U, V (the order in the 'decomp' vector)
void svdDenoising(std::string inputImage, unsigned width, unsigned height, std::vector<Eigen::MatrixXd> &channelMatrices, std::vector<std::vector<Eigen::MatrixXd>> &outMatrices,
    std::vector<std::vector<fs::path>> rgbPaths, int k, int id){
    
    /* reconstruction step:
    in: original matrix M (of size m X n), kSVD 
    
    U_k : first m X kSVD blocks
    Σ_k: kSVD X kSVD diagonal matrix
    V_k: kSVD X n
    --------------------------------
    
    out reconstructed final matrix M' (m X n)
    */

    std::optional<std::vector<std::vector<Eigen::MatrixXd>>> svd = svdChannels(inputImage, width, height, channelMatrices, outMatrices, rgbPaths);
    if (!svd){
        std::cerr << "Error before denoising" << std::endl;
    }

    std::vector<Eigen::MatrixXd> truncatedChannels;
    for (int i = 0; i < svd->size(); i++){

        std::vector<Eigen::MatrixXd> currentDecomp = svd->at(i);
        // B, U, V (the order in the 'decomp' vector)
        Eigen::MatrixXd finalB = truncateB(currentDecomp[0], k);
        Eigen::MatrixXd finalU = truncateU(currentDecomp[1], k);
        Eigen::MatrixXd finalV = truncateV(currentDecomp[2], k);
        
        truncatedChannels.push_back((finalU * finalB) * finalV);
        assert(truncatedChannels[i].rows() == height && truncatedChannels[i].cols() == width);
    }

    std::vector<unsigned char> newImage;
    std::string path = fmt::format("/home/gabriel/Documents/HolyC/SVD_image_denoising/images/denoised{}_{}.png", id, k);
    reconstructImage(truncatedChannels, newImage, path, height, width, id);
}


int main(){

    std::vector<unsigned char> image;
    unsigned width; // img_row
    unsigned height; // img_col

    std::string myImage = "/home/gabriel/Documents/HolyC/SVD_image_denoising/images/noisy_mu0_std10.png";


    // printImage(image, height, width);

    // adding noise
    // addNoise(image, height, width);
    
    // In
    Eigen::MatrixXd R;
    Eigen::MatrixXd G;
    Eigen::MatrixXd B;
    std::vector<Eigen::MatrixXd> channelMatrices = {R, G, B};

    // OUT

    std::vector<Eigen::MatrixXd> outR;
    std::vector<Eigen::MatrixXd> outG;
    std::vector<Eigen::MatrixXd> outB;
    std::vector<std::vector<Eigen::MatrixXd>> outMatrices = {outR, outG, outB};

    // saving paths for SVD decomposition
    std::vector<fs::path> channelR = {"matrices/channelR/B_final.csv", "matrices/channelR/U_final.csv", "matrices/channelR/V_final.csv"};
    std::vector<fs::path> channelG = {"matrices/channelG/B_final.csv", "matrices/channelG/U_final.csv", "matrices/channelG/V_final.csv"};
    std::vector<fs::path> channelB = {"matrices/channelB/B_final.csv", "matrices/channelB/U_final.csv", "matrices/channelB/V_final.csv"};

    std::vector<std::vector<fs::path>> rgbPaths = {channelR, channelG, channelB};


    // tests
    //bidiagonalizeTest();
    
    //golubKahanTest();

    //std::vector<Eigen::MatrixXd> matrices = bidiagonalize(A, A.rows(), A.cols);

    fs::path imgPath = myImage;

    int k = 50; // here, I should try more k's to see which one is the best
    int idx = 1; // the first image
    svdDenoising(myImage, width, height, channelMatrices, outMatrices, rgbPaths, k, idx);

    //golubKahanTest();
    
    //matrixToCsv(B_mine, csvPathB, true); 
    //matrixToCsv(U, csvPathU, true); 
    //matrixToCsv(V_transpose, csvPathV, true); 
    
    /*
    Eigen::JacobiSVD<Eigen::MatrixXd> svdDecomp(R_channel, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Eigen::MatrixXd R_svd = svdDecomp.matrixU()
                          * svdDecomp.singularValues().asDiagonal()
                          * svdDecomp.matrixV().transpose();
    
    Eigen::MatrixXd R_recon = U_final * B_final * V_final.transpose();

    double d2NormDiff = (R_svd - R_recon).norm();
    double d2Norm = R_channel.norm();

    std::cout << "Relative construction error " << d2NormDiff / (d2Norm * 65025) << std::endl;
    
    auto mySingularValues = B_final.diagonal();
    auto eigenSingularValues = svdDecomp.singularValues();

    std::cout << mySingularValues.transpose() << std::endl;

    std::cout << std::endl << " S E P A R A T O R " << std::endl;

    std::cout << eigenSingularValues.transpose() << std::endl;
    */
    


    return 0;
}