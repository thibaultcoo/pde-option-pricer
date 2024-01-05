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
    // black-scholes pricing sandbox
    double spot = 70;
    double strike = 110;
    double rate = 0.13;
    double divs = 0.00;
    double repo = 0.00;
    double vol = 0.34;
    double matu = 2;
    double price_bs;
    double price_pde;

    // creating an option object and pricing the corresponding call
    blackScholes option(spot, strike, rate, divs, repo, vol, matu);
    price_bs = option.callOptionPrice();

    size_t priceGridSize = 50;
    size_t timeGridSize = 50;
    double multiplier = 4;

    // the conditions are to be filled with the grid (choosing custom conditions is not possible thus far)
    matrix terminalCondition(priceGridSize, 1);
    matrix lowerBoundaries(timeGridSize, 1);
    matrix upperBoundaries(timeGridSize, 1);

    // defining the coefficient functions for the algorithm resolution
    pricerPDE::coeffFunction a = [rate](double t, double x)->double {return -1 * rate;};
    pricerPDE::coeffFunction b = [rate, divs, repo](double t, double x)->double {return (rate - divs - repo) * x;};
    pricerPDE::coeffFunction c = [vol](double t, double x)->double {return 0.5 * vol * vol * x * x;};
    pricerPDE::coeffFunction d = [](double t, double x)->double {return 0.0;};

    // creating a pricing object using pde and implicit finite difference methods
    pricerPDE pricer(spot, strike, matu, vol, rate, divs, repo, multiplier, terminalCondition, lowerBoundaries, upperBoundaries, a, b, c, d);
    price_pde = pricer.callOptionPrice();

    // comparing the results following both methods to study convergence
    std::cout << " " << price_bs << std::endl;
    std::cout << " " << price_pde << std::endl;
}