#include "tests.hpp"
#include "imageUtils.hpp"
#include "bidiagonalization.hpp"
#include "golubKahan.hpp"
#include <iostream>

void bidiagonalizeTest(){
    Eigen::MatrixXd B_test(5, 5);
    B_test << 0.3, 0.7, 0.8, 0.1, 0.3,
        0.25, 0.8, 0.8, 0.5, 0.6,
        0.35 ,0.9 , 0.738, 0.87, 0.33,
        0.45 ,0.6, 0, 0.3634, 0.9,
        0.55, 0.5, 0.19, 0.22, 0.9;

    std::vector<Eigen::MatrixXd> matrices = bidiagonalize(B_test, 5, 5);

    Eigen::MatrixXd U = matrices[0];
    Eigen::MatrixXd B_mine = matrices[1];
    Eigen::MatrixXd V_transpose = matrices[2];
    std::cout << "U = " << std::endl << U << std::endl;
    std::cout << "B = " << std::endl << B_mine << std::endl;
    std::cout << "V_transpose = " << std::endl << V_transpose << std::endl;
}

void golubKahanTest(){

    fs::path csvFinalB = "/home/gabriel/Documents/HolyC/SVD_image_denoising/matrices/channelR/B_final_test.csv";

    Eigen::MatrixXd B_test(5, 5);
    Eigen::MatrixXd U_test (5, 5);
    Eigen::MatrixXd V_test(5, 5);

    B_test << -8.02055, 129.93, 0, 0, 0,
        0, -83.5067, 8.72759, 0, 0,
        0 ,0 ,7.01538, 8.72961, 0,
        0 ,0, 0, -12.3634, 9.37061,
        0, 0, 0, 0, -9.53761;

    U_test << -0.089476, 0.0258819, 0.075377, -0.0699437, 0.0618868,
    -0.0860535,0.0208654,0.078938,-0.0730571,0.0582035,
    -0.0958323,0.0333659,0.0655085,-0.0646888,0.0517716,
    -0.0968101,0.0348304,0.0577696,-0.0538016,0.0463264,
    -0.0904539,0.0252435,0.0793154,-0.0668693,0.0581317;

    V_test << 0.0556736, 0.0534652, 0.0549524, 0.053282, 0.0464544,
    -0.374102, -0.25997, -0.262998, -0.177746, -0.0951319,
    0.0906982, -0.0372941, -0.00531366, -0.0473925, 0.041103,
    -0.103568, 0.103578, 0.0771154, 0.11051, 0.0401614,
    0.0611464, -0.0704854, -0.0227363, -0.0229629, -0.0318511;

    std::vector<Eigen::MatrixXd> newValues = svdGolubKahan(B_test, U_test, V_test);
    matrixToCsv(newValues[0], csvFinalB, true);
}