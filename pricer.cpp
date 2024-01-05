#include <iostream>
#include <vector>
#include <cmath>
#include "pricer.h"
#include "matrix.h"

// PDE pricer object constructor
pricerPDE::pricerPDE(double spot, double strike, double matu, double vol, 
                     double rate, double divs, double repo, double multiplier, 
                     matrix& terminalCondition, 
                     matrix& lowerBoundaries, matrix& upperBoundaries,
                     coeffFunction a, coeffFunction b, coeffFunction c, coeffFunction d)
{
    p_spot = spot;
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
    p_isDefaultLower = lowerBoundaries.isEmpty();
    p_isDefaultUpper = upperBoundaries.isEmpty();
    p_isDefaultTerminal = terminalCondition.isEmpty();
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
    double price = extractPrice();

    return price;
}

// extracts the final price of the option
double pricerPDE::extractPrice()
{
    int xIdx = findClosestIdx(this->p_spotGrid, this->p_spot);
    
    return interpo(xIdx);
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

// produces the linear interpolation in x
double pricerPDE::interpo(int xIdx)
{
    if (xIdx < 0 || xIdx > this->p_m) {
        throw std::out_of_range("Interpo indices out of range");
    }

    double xLower = this->p_spotGrid.m_M[xIdx][0];
    double xUpper;

    // additional security in case index is out of bounds
    if (xIdx < this->p_spotGrid.m_M.size() - 2) {
        xUpper = this->p_spotGrid.m_M[xIdx + 1][0];
    } else {
        return this->p_priceGrid.m_M[xIdx][0];
    }

    double xProp = (this->p_spot - xLower) / (xUpper - xLower);

    // interpolating the values
    double valueLower = this->p_priceGrid.m_M[xIdx + 1][0];
    double valueUpper = this->p_priceGrid.m_M[xIdx + 2][0];
    double interpolatedValue = valueLower * (1 - xProp) + valueUpper * xProp;

    return interpolatedValue;
}

// application of the chosen finite difference scheme
void pricerPDE::applyCrankNicholson()
{
    // initializing the variable grids sizes (time and spot discretization)
    size_t p_n = this->p_n;
    size_t p_m = this->p_m;

    // output values of the given four functions across time
    double iterSpot, iterTime;
    double iterA, iterB, iterC, iterD;

    // coefficient value coherent with a Crank-Nicholson scheme (half implicit/explicit resolution)
    double theta = 0.5;

    // creating iterable matrix U corresponding to price with boundary values
    matrix U(p_m - 1,1);
    double uFirst;
    double prev_uFirst;

    // initalizing terminal values
    uFirst = this->p_terminalCondition.m_M[0][0];

    for(size_t i = 0; i < this->p_m - 1; i++){
        U.m_M[i][0] = this->p_terminalCondition.m_M[i][0];
    }

    // Loop through values of time backward 
    for (size_t k = this->p_n - 2; k > 0; k--) {
        // apply boundaries conditions
        prev_uFirst = uFirst;

        setupBoundary(k / this->p_n);
        uFirst = this->p_lowerBoundaries.m_M[k][0];

        matrix P(p_m - 1, p_m - 1);
        matrix Q(p_m - 1, p_m - 1);
        matrix V(p_m - 1, 1);

        // get current value of time
        iterTime = this->p_timeGrid.m_M[k][0];
         
        for(size_t i = 0; i < this->p_m - 1; i++){
            // get current value of spot
            iterSpot = this->p_spotGrid.m_M[i + 1][0];
            
            // compute only once values of functions a,b,c,d
            iterA = this->p_a(iterTime, iterSpot);
            iterB = this->p_b(iterTime, iterSpot);
            iterC = this->p_c(iterTime, iterSpot);
            iterD = this->p_d(iterTime, iterSpot);

            // fill iterative matrix P and Q
            for(size_t j = 0; j < this->p_m - 1; j++){
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

        // retrieve previous state of matrix
        U = P.inversion() * (Q * U + V) * -1;
    }
    // returning output into variable
    this->p_priceGrid.m_M[0][0] = this->p_lowerBoundaries.m_M[0][0];
    this->p_priceGrid.m_M[p_m-1][0] = this->p_upperBoundaries.m_M[0][0];

    for(size_t i = 0; i < this->p_m - 1; i++){
        this->p_priceGrid.m_M[i+1][0] = U.m_M[i][0];
    }
}

// sets up the discretized time/spot/price grids
void pricerPDE::setupGrid()
{
    // initializing the variable grids
    size_t p_n = this->p_n;
    size_t p_m = this->p_m;

    matrix timeGrid(p_n, 1);
    matrix spotGrid(p_m, 1);
    matrix priceGrid(p_n, p_m);

    this->p_timeGrid = timeGrid;
    this->p_spotGrid = spotGrid;
    this->p_priceGrid = priceGrid;

    // initializing the steps
    this->p_dt = this->p_matu / this->p_n;
    this->p_dS = this->p_multiplier * this->p_vol * 2 * this->p_strike / this->p_m;

    // fills the time grid
    for (int i = 0; i < this->p_n; i++) {
        this->p_timeGrid.m_M[i][0] = i * this->p_dt;
    }

    // fills the spot grid
    for (int i = 0; i < this->p_m; i++) {
        this->p_spotGrid.m_M[i][0] = this->p_strike - this->p_multiplier * this->p_vol * this->p_strike + i * this->p_dS;

        // applies the natural lower boundary condition on the spot (assume to be positive or null)
        if (this->p_spotGrid.m_M[i][0] < 0) {this->p_spotGrid.m_M[i][0] = 0;};
    }

    setupBoundary(this->p_matu);
    setupTerminal();
}

// fills the lower and upper boundary conditions
void pricerPDE::setupTerminal()
{
    for (int i = 0; i < this->p_m; i++) {
        if (this->p_isDefaultTerminal) {this->p_terminalCondition.m_M[i][0] = std::max(this->p_spotGrid.m_M[i][0] - this->p_strike, 0.0);}
    }
}

// fills the terminal condition
void pricerPDE::setupBoundary(double remaining_t)
{
    for (int i = 0; i < this->p_n; i++) {
        if (this->p_isDefaultLower) {this->p_lowerBoundaries.m_M[i][0] = 0;}
        if (this->p_isDefaultUpper) {this->p_upperBoundaries.m_M[i][0] = this->p_spotGrid.m_M[i][0];}
        if (this->p_isDefaultUpper) {this->p_upperBoundaries.m_M[i][0] *= std::exp((this->p_repo - this->p_divs - this->p_rate) * (remaining_t));}
    }
}