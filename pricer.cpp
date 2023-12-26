#include <iostream>
#include <vector>
#include "pricer.h"

// PDE pricer object constructor
pricerPDE::pricerPDE(double strike, double matu, double vol, 
                     double rate, double divs, double repo,
                     double multiplier, int m, int n)
{
    p_multiplier = multiplier;
    p_strike = strike;
    p_matu = matu;
    p_vol = vol;
    p_rate = rate;
    p_divs = divs;
    p_repo = repo;
    p_m = m;
    p_n = n;

    setupGrid();
    applyTerminalConditions();
}

// outputs the final price
double pricerPDE::getPrice()
{
    return 0.0;
}

// sets up the discretized time/spot/price grids
void pricerPDE::setupGrid()
{
    // initializing the variable grids
    this->p_timeGrid.resize(this->p_n+1);
    this->p_spotGrid.resize(this->p_m+1);

    // initializing the steps
    double dt = this->p_matu / this->p_n;
    double dS = this->p_multiplier * this->p_vol * 2 * this->p_strike / this->p_m;

    // fills the time grid
    for (int i = 0; i < this->p_n+1; i++) {
        p_timeGrid[i] = i * dt;
    }

    // fills the spot grid
    for (int i = 0; i < this->p_m+1; i++) {
        p_spotGrid[i] = this->p_strike - this->p_multiplier * this->p_vol * this->p_strike + i * dS;
    }

    // initializing the grid that will be filled with option prices
    this->p_priceGrid.resize(this->p_n+1, std::vector<double>(this->p_m+1, 0.0));
}

// sets up the terminal condition
void pricerPDE::applyTerminalConditions()
{
    for (int i = 0; i < this->p_m+1; i++) {
        double spot = this->p_spotGrid[i];
        p_priceGrid[this->p_n][i] = std::max(spot - this->p_strike, 0.0);
    }
}