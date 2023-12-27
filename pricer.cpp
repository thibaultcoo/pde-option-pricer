#include <iostream>
#include <vector>
#include <cmath>
#include "pricer.h"
#include "matrix.h"

// PDE pricer object constructor
pricerPDE::pricerPDE(double spot, double strike, double matu, double vol, 
                     double rate, double divs, double repo, double multiplier, 
                     const matrix& terminalCondition, const matrix& boundaryConditions,
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
    p_n = static_cast<int>(boundaryConditions.m_M.size());

    setupGrid();
    setTerminalCondition(terminalCondition);
    setBoundaryConditions(boundaryConditions);
}

// outputs the final price
double pricerPDE::callOptionPrice()
{
    // iterate backwards by applying Crank-Nicholson step
    for (int t = this->p_n - 1; t > 0; --t)
        applyCrankNicholson();

    // now that the solution is in the grid, we extract it (directly or from interpo)
    double price = extractPrice();

    return price;
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

void applyCrankNicholson()
{

}

double pricerPDE::interpo(int tIdx, int xIdx)
{
    if (tIdx < 0 || tIdx > this->p_n || xIdx < 0 || xIdx > this->p_m) {
        throw std::out_of_range("Interpo indices out of range");
    }

    // bilinear interpolation
    double tProp = (this->p_matu - this->p_timeGrid.m_M[tIdx][0]) / p_dt;
    double xProp = (this->p_spot - this->p_spotGrid.m_M[xIdx][0]) / p_dS;

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

void pricerPDE::setTerminalCondition(const matrix& terminalCondition)
{
    for (int i = 0; i < this->p_m; i++) {
        p_priceGrid.m_M[this->p_n][i] = terminalCondition.m_M[i][0];
    }
}

void pricerPDE::setBoundaryConditions(const matrix& boundaryConditions)
{
    for (int i = 0; i < this->p_n; i++) {
        p_priceGrid.m_M[i][0] = boundaryConditions.m_M[i][0];
        p_priceGrid.m_M[i][this->p_m] = boundaryConditions.m_M[i][1];
    }
}