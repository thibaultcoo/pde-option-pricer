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
    // application of the scheme that will fill the price grid
    applyCrankNicholson();

    // now that the solution is in the grid, we extract it from interpolation
    double price = extractPrice();

    return price;
}

// extracts the final price of the option
double pricerPDE::extractPrice()
{
    // we look for the closest value from the spot and retrieve its index
    int xIdx = findClosestIdx(this->p_spotGrid, this->p_spot);
    
    // having located the final price area, we thus interpolate to approximate its value
    return interpo(xIdx);
}

// finds the closest index from the one provided on the grid
int pricerPDE::findClosestIdx(const matrix& grid, double value)
{
    for (int i = 0; i < grid.m_M.size(); i++) {
        // as soon as the spot seen on the grid is bigger than our spot, we know what our boundaries are
        if (value - grid.m_M[i][0] < 0) {
            return i - 1;
        }
    }
    // corner case: the spot is outside the grid
    return grid.m_M.size() - 2;
}

// produces the linear interpolation in x
double pricerPDE::interpo(int xIdx)
{
    // isolating the spot boundaries for the interpolation
    double xLower = this->p_spotGrid.m_M[xIdx][0];
    double xUpper = this->p_spotGrid.m_M[xIdx + 1][0];

    // isolating the corresponding prices
    double priceLower = this->p_priceGrid.m_M[xIdx][0];
    double priceUpper = this->p_priceGrid.m_M[xIdx + 1][0];

    // computing the resulting interpolated price
    double xProp = (this->p_spot - xLower) / (xUpper - xLower);
    double interpolatedPrice = priceLower * (1 - xProp) + priceUpper * xProp;

    return interpolatedPrice;
}

// application of the chosen finite difference scheme
void pricerPDE::applyCrankNicholson()
{
    // initializing the variable grids sizes (time and spot discretization respectively)
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
            if (i == 0) {
                V.m_M[i][0] += uFirst * (-iterB / (2 * this->p_dS) + theta * iterC / (std::pow(this->p_dS, 2))) + prev_uFirst * (1 - theta) * iterC / (std::pow(this->p_dS, 2));
            }
            else if (i == this->p_m - 2) {
                V.m_M[i][0] += uFirst * (iterB / (2 * this->p_dS) + theta * iterC / (std::pow(this->p_dS, 2))) + prev_uFirst * (1 - theta) * iterC / (std::pow(this->p_dS, 2));
            }
        }

        // retrieve previous state of matrix
        U = P.inversion() * (Q * U + V) * -1;

        // lets us know what percentage of the code has compiled
        std::cout << "Code is at " << (1.0 - static_cast<double>(k) / static_cast<double>(this->p_n)) * 100 << "%" << std::endl;

    }
    // cleans output when compilation is done
    for (int i = 0; i < 20; ++i) {std::cout << "\n";}

    for (size_t i = 0; i < this->p_m - 1; i++) {
        this->p_priceGrid.m_M[i][0] = U.m_M[i][0];
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
    matrix priceGrid(p_m, p_n);

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

    // setting up the boundary conditions (will be set at their default value if input is empty)
    setupBoundary(this->p_matu);
    setupTerminal();
}

// fills the lower and upper boundary conditions (default if nothing is input by Prof.)
void pricerPDE::setupTerminal()
{
    double discount = 1;

    for (int i = 0; i < this->p_m; i++) {
        if (this->p_isDefaultTerminal) {this->p_terminalCondition.m_M[i][0] = std::max(this->p_spotGrid.m_M[i][0] * discount - this->p_strike, 0.0);}
    }
}

// fills the terminal condition (default if nothing is input by Prof.)
void pricerPDE::setupBoundary(double remaining_t)
{
    for (int i = 0; i < this->p_n; i++) {
        if (this->p_isDefaultLower) {this->p_lowerBoundaries.m_M[i][0] = 0;}
        if (this->p_isDefaultUpper) {this->p_upperBoundaries.m_M[i][0] = this->p_spotGrid.m_M[i][0];}
        if (this->p_isDefaultUpper) {this->p_upperBoundaries.m_M[i][0] *= std::exp((this->p_repo + this->p_divs - this->p_rate) * (remaining_t));}
    }
}