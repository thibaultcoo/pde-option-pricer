#include <iostream>
#include <cstddef>
#include <string>
#include <vector>
#include <random>
#include <chrono>
#include "matrix.h"

// sandbox cpp file for testing matrix class (working well, still corner cases to look into)
int main()
{
    size_t nRows = 6;
    size_t nCols = 6;
    matrix mat(nRows, nCols);
    matrix invert(nRows, nCols);
    matrix testingInvert(nRows, nCols);

    auto currentTime = std::chrono::system_clock::now().time_since_epoch().count();
    std::mt19937 rng(currentTime);
    std::uniform_real_distribution<double> dist(-0.001, 0.001);

    // randomly filling the initial matrix
    for (int i = 0; i < nRows; i++) {
        for (int j = 0; j < nCols; j++) {
            mat.m_M[i][j] = dist(rng);
            std::cout << " " << mat.m_M[i][j];
        }
        std::cout << std::endl;
    }
    
    invert = mat.inversion();
    testingInvert = mat * invert;
    
    double sum = 0;
    for (int i = 0; i < nRows; i++) {
        for (int j = 0; j < nCols; j++) {
            sum += testingInvert.m_M[i][j];
            std::cout << " " << testingInvert.m_M[i][j];
        }
        std::cout << std::endl;
    }
    
    if (sum = nCols) {
        std::cout << "Inversion was successful" << std::endl;
    } else {
        std::cout << "Inversion failed !" << std::endl;
    }

}

double random_double(double min, double max) {
    // Random number generator using the Mersenne Twister engine
    static std::mt19937 rng(std::random_device{}());

    // Uniform distribution in range [min, max]
    std::uniform_real_distribution<double> dist(min, max);

    // Generate and return the random double
    return dist(rng);
}