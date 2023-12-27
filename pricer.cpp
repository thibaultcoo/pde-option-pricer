#include <iostream>
#include <vector>
#include "pricer.h"
#include "matrix.h"

// PDE pricer object constructor
pricerPDE::pricerPDE(double strike, double matu, double vol, 
                     double rate, double divs, double repo,
                     double multiplier, const matrix& terminalCondition,
                     const matrix& boundaryConditions, 
                     coeffFunction a, coeffFunction b, coeffFunction c)
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
    // applyCrankNicholson();

    return 0.0;
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
    double dt = this->p_matu / this->p_n;
    double dS = this->p_multiplier * this->p_vol * 2 * this->p_strike / this->p_m;

    // fills the time grid
    for (int i = 0; i < this->p_n+1; i++) {
        p_timeGrid.m_M[i][0] = i * dt;
    }

    // fills the spot grid
    for (int i = 0; i < this->p_m+1; i++) {
        p_spotGrid.m_M[i][0] = this->p_strike - this->p_multiplier * this->p_vol * this->p_strike + i * dS;
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