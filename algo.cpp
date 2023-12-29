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

    size_t priceGridSize = 10;
    size_t timeGridSize = 15;
    double multiplier = 1;

    // the conditions are to be filled within the grid (they are given by the call payoff as well as the spot grid)
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

    std::cout << price_bs << std::endl;
    std::cout << price_pde << std::endl;
}