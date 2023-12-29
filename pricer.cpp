#include <iostream>
#include <vector>
#include <cmath>
#include "pricer.h"
#include "matrix.h"

// PDE pricer object constructor
pricerPDE::pricerPDE(double spot, double strike, double matu, double vol, 
                     double rate, double divs, double repo, double multiplier, 
                     const matrix& terminalCondition, 
                     const matrix& lowerBoundaries, const matrix& upperBoundaries,
                     coeffFunction a, coeffFunction b, coeffFunction c, coeffFunction d)
{
    p_strike = strike;
    p_matu = matu;
    p_vol = vol;
    p_rate = rate;
    p_divs = divs;
    p_repo = repo;
    p_multiplier = multiplier;
    p_m = static_cast<int>(terminalCondition.m_M.size());
    p_n = static_cast<int>(lowerBoundaries.m_M.size());
    p_lowerBoundaries = lowerBoundaries;
    p_upperBoundaries = upperBoundaries;
    p_terminalCondition = terminalCondition;
    p_a = a;
    p_b = b;
    p_c = c;
    p_d = d;
    
    setupGrid();
}

// outputs the final price
double pricerPDE::callOptionPrice()
{
    // application of the algorithm to fill the price grid
    applyCrankNicholson();

    // now that the solution is in the grid, we extract it (directly or from interpo)
    // double price = extractPrice();

    for (size_t i = 0; i < p_m - 1; i++) {
        std::cout << " " << this->p_priceGrid.m_M[i][0];
    }

    return 0.0;
}

// extracts the final price of the option
double pricerPDE::extractPrice()
{
    // wrong indices, work on the logic behind
    int tIdx = findClosestIdx(this->p_timeGrid, 0.0);
    int xIdx = findClosestIdx(this->p_spotGrid, this->p_spot);
    
    return interpo(tIdx, xIdx);
}

// finds the closest index from the one provided on the grid
int pricerPDE::findClosestIdx(const matrix& grid, double value)
{
    double disFromElement = std::abs(value - grid.m_M[0][0]);
    int resIdx;

    for (int i = 1; i < grid.m_M.size(); i++) {
        if (disFromElement > std::abs(value - grid.m_M[i][0])) {
            resIdx = i;
            disFromElement = std::abs(value - grid.m_M[i][0]);
        }
    }
    return resIdx;
}

void pricerPDE::applyCrankNicholson()
{
    // initializing the variable grids sizes (time and spot discretization)
    size_t p_n = this->p_n;
    size_t p_m = this->p_m;

    // output values of the given four functions across time
    double iterSpot, iterTime;
    double iterA, iterB, iterC, iterD;

    double theta = 0.0;

    // creating iterable matrix U corresponding to price with boundary values
    matrix U(p_m - 1,1);
    double uFirst, uLast;
    double prev_uFirst, prev_uLast;

    // initalizing terminal values
    uFirst = this->p_terminalCondition.m_M[0][0];
    uLast = this->p_terminalCondition.m_M[p_m - 1][0];

    for(size_t i = 0; i < p_m - 1; i++){
        U.m_M[i][0] = this->p_terminalCondition.m_M[i][0];
    }

    // Loop through values of time backward 
    for (size_t k = p_n - 2; k > 0; k--) {
        // apply boundaries counditions
        prev_uFirst = uFirst;
        prev_uLast = uLast;
        uFirst = this->p_lowerBoundaries.m_M[k][0];
        uLast = this->p_upperBoundaries.m_M[k][0];

        matrix P(p_m - 1, p_m - 1);
        matrix Q(p_m - 1, p_m - 1);
        matrix V(p_m - 1, 1);

        // get current value of time
        iterTime = this->p_timeGrid.m_M[k][0];
         
        for(size_t i = 0; i < p_m - 1; i++){
            // get current value of spot
            iterSpot = this->p_spotGrid.m_M[i + 1][0];
            
            // compute only once values of functions a,b,c,d
            iterA = this->p_a(iterTime, iterSpot);
            iterB = this->p_b(iterTime, iterSpot);
            iterC = this->p_c(iterTime, iterSpot);
            iterD = this->p_d(iterTime, iterSpot);

            // fill iterative matrix P and Q
            for(size_t j = 0; j < p_m - 1; j++){
                // diagonal terms
                if(j == i){
                    P.m_M[i][j] = iterA - (1 / this->p_dt + 2 * theta * iterC / (std::pow(this->p_dS, 2))); 
                    Q.m_M[i][j] = 1 / this->p_dt - 2 * (1 - theta) * iterC / (std::pow(this->p_dS, 2)); 
                }
                // sup diagonal terms
                else if(j == i + 1){
                    P.m_M[i][j] = iterB / (2 * this->p_dS) + theta * iterC / (std::pow(this->p_dS, 2)); 
                    Q.m_M[i][j] = (1 - theta) * iterC / (std::pow(p_dS, 2)); 
                }
                // inf diagonal terms
                else if(j == i - 1){
                    P.m_M[i][j] = -1 * iterB / (2 * this->p_dS) + theta * iterC / (std::pow(this->p_dS, 2)); 
                    Q.m_M[i][j] = (1 - theta) * iterC / (std::pow(this->p_dS, 2)); 
                }
            }
            V.m_M[i][0] = iterD;

            // adjust V with boundary conditions
            if(i == 0 || i == this->p_m - 2){
                V.m_M[0][0] += prev_uFirst * (iterB / (2 * this->p_dS) + theta * iterC / (std::pow(this->p_dS, 2)));
                V.m_M[0][0] += uFirst * (1 - theta) * iterC / (std::pow(this->p_dS, 2));
            }
        }

        // retreive previous state of matrix
        U = P.inversion() * (Q * U + V) * -1;
    }
    // returning output into variable
    matrix finalPrice(p_m + 1, 1);

    finalPrice.m_M[0][0] = this->p_lowerBoundaries.m_M[0][0];
    finalPrice.m_M[p_m][0] = this->p_upperBoundaries.m_M[0][0];

    for(size_t i = 0; i < p_m - 1; i++){
        finalPrice.m_M[i+1][0] = U.m_M[i][0];
    }

    this->p_priceGrid = finalPrice;
}

double pricerPDE::interpo(int tIdx, int xIdx)
{
    if (tIdx < 0 || tIdx > this->p_n || xIdx < 0 || xIdx > this->p_m) {
        throw std::out_of_range("Interpo indices out of range");
    }

    // bilinear interpolation
    double tProp = (this->p_matu - this->p_timeGrid.m_M[tIdx][0]) / this->p_dt;
    double xProp = (this->p_spot - this->p_spotGrid.m_M[xIdx][0]) / this->p_dS;

    double res_1 = (1 - tProp) * (1 - xProp) * this->p_priceGrid.m_M[tIdx][xIdx];
    double res_2 = tProp * (1 - xProp) * this->p_priceGrid.m_M[tIdx + 1][xIdx];
    double res_3 = (1 - tProp) * xProp * this->p_priceGrid.m_M[tIdx][xIdx + 1];
    double res_4 = tProp * xProp * this->p_priceGrid.m_M[tIdx + 1][xIdx + 1];

    return res_1 + res_2 + res_3 + res_4;
}

// sets up the discretized time/spot/price grids
void pricerPDE::setupGrid()
{
    // initializing the variable grids
    size_t p_n = this->p_n;
    size_t p_m = this->p_m;

    matrix timeGrid(p_n+1, 1);
    matrix spotGrid(p_m+1, 1);
    matrix priceGrid(p_n+1, p_m+1);

    this->p_timeGrid = timeGrid;
    this->p_spotGrid = spotGrid;
    this->p_priceGrid = priceGrid;

    // initializing the steps
    double p_dt = this->p_matu / this->p_n;
    double p_dS = this->p_multiplier * this->p_vol * 2 * this->p_strike / this->p_m;

    // fills the time grid
    for (int i = 0; i < this->p_n+1; i++) {
        p_timeGrid.m_M[i][0] = i * p_dt;
    }

    // fills the spot grid
    for (int i = 0; i < this->p_m+1; i++) {
        p_spotGrid.m_M[i][0] = this->p_strike - this->p_multiplier * this->p_vol * this->p_strike + i * p_dS;
    }   
}