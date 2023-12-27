#include <iostream>
#include <cstddef>
#include <vector>
#include <random>
#include <chrono>
#include <functional>
#include "matrix.h"
#include "bs.h"
#include "pricer.h"

// sandbox cpp
int main()
{
    // matrix manipulation sandbox
    size_t nRows = 2;
    size_t nCols = 2;
    matrix mat(nRows, nCols);
    matrix invert(nRows, nCols);
    matrix testingInvert(nRows, nCols);

    // seed randomizer for double generation
    auto currentTime = std::chrono::system_clock::now().time_since_epoch().count();
    std::mt19937 rng(currentTime);
    std::uniform_real_distribution<double> dist(-1, 1);

    // randomly filling the initial matrix
    for (int i = 0; i < nRows; i++) {
        for (int j = 0; j < nCols; j++) {
            mat.m_M[i][j] = dist(rng);
            std::cout << " " << mat.m_M[i][j];
        }
        std::cout << std::endl;
    }
    
    // inverting the matrix
    invert = mat.inversion();

    // testing if inversion worked welled
    testingInvert = mat * invert;

    for (int i = 0; i < nRows; i++) {
        for (int j = 0; j < nCols; j++) {
            std::cout << " " << testingInvert.m_M[i][j];
        }
        std::cout << std::endl;
    }

    // black-scholes pricing sandbox
    double spot = 100;
    double strike = 90;
    double rate = 0.05;
    double divs = 0.00;
    double repo = 0.00;
    double vol = 0.25;
    double matu = 1;
    double price_bs;
    double price_pde;

    // creating an option object and pricing the corresponding call
    blackScholes option(spot, strike, rate, divs, repo, vol, matu);
    price_bs = option.callOptionPrice();

    double multiplier = 0.1;
    matrix terminalCondition(10, 1);
    matrix boundaryConditions(15, 2);

    std::function<double(double, double)> a = [](double t, double x) -> double {
        return x - t;
    };

    std::function<double(double, double)> b = [](double t, double x) -> double {
        return x + t;
    };

    std::function<double(double, double)> c = [](double t, double x) -> double {
        return x * t;
    };

    // creating a pricing object using pde and implicit finite difference methods
    pricerPDE pricer(strike, matu, vol, rate, divs, repo, multiplier, terminalCondition, boundaryConditions, a, b, c);
    price_pde = pricer.callOptionPrice();

    std::cout << price_pde << std::endl;
}