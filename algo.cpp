#include <iostream>
#include <cstddef>
#include <vector>
#include "matrix.h"
#include "bs.h"
#include "pricer.h"

matrix customTerminal(size_t priceGridSize);
matrix customLower(size_t timeGridSize);
matrix customUpper(size_t timeGridSize);

int main()
{
    // Prof. pricing parameters
    double spot = 50;
    double strike = 56;
    double rate = 0.08;
    double divs = 0.05;
    double repo = 0.0;
    double vol = 0.23;
    double matu = 0.70;

    // creating an option object and pricing the corresponding call
    blackScholes option(spot, strike, rate, divs, repo, vol, matu);
    double price_bs = option.callOptionPrice();

    // discretization steps for the finite difference price resolution
    size_t priceGridSize = 45;
    size_t timeGridSize = 45;

    // chosen to generate a coherent spot grid
    double multiplier = 4;

    // initializing the conditions/boundaries vectors
    matrix terminalCondition(priceGridSize, 1);
    matrix lowerBoundaries(timeGridSize, 1);
    matrix upperBoundaries(timeGridSize, 1);

    // ability for Prof. to upload custom conditions (default conditions used otherwise)
    terminalCondition = customTerminal(priceGridSize);
    lowerBoundaries = customLower(timeGridSize);
    upperBoundaries = customUpper(timeGridSize);

    // defining the coefficient functions for the algorithm resolution
    pricerPDE::coeffFunction a = [rate](double t, double x)->double {return -1 * rate;};
    pricerPDE::coeffFunction b = [rate, divs, repo](double t, double x)->double {return (rate - divs - repo) * x;};
    pricerPDE::coeffFunction c = [vol](double t, double x)->double {return 0.5 * vol * vol * x * x;};
    pricerPDE::coeffFunction d = [](double t, double x)->double {return 0.0;};

    // creating a pricing object using pde and implicit finite difference methods
    pricerPDE pricer(spot, strike, matu, vol, rate, divs, repo, multiplier, terminalCondition, lowerBoundaries, upperBoundaries, a, b, c, d);
    double price_pde = pricer.callOptionPrice();

    // comparing the results following both methods to study convergence
    std::cout << "The Black-Scholes call option price is " << price_bs << "." << std::endl;
    std::cout << "The Finite-Difference call option price is " << price_pde << "." << std::endl;
}

// Prof. : to use a custom terminal condition
matrix customTerminal(size_t priceGridSize)
{
    // initializes the resulting matrix
    matrix resTerminal(priceGridSize, 1);

    // gives the custom-valued vector (replace with Prof. vector)
    std::vector<double> customTerminal(static_cast<int>(priceGridSize), 0.0);

    // returns the converted custom vector into a usable matrix
    return resTerminal.toMatrix(customTerminal);
}

// Prof. : to use custom lower boundary conditions
matrix customLower(size_t timeGridSize)
{
    // initializes the resulting matrix
    matrix resLower(timeGridSize, 1);

    // gives the custom-valued vector (replace with Prof. vector)
    std::vector<double> customLower(static_cast<int>(timeGridSize), 0.0);

    // returns the converted custom vector into a usable matrix
    return resLower.toMatrix(customLower);
}

// Prof. : to use custom upper boundary conditions
matrix customUpper(size_t timeGridSize)
{
    // initializes the resulting matrix
    matrix resUpper(timeGridSize, 1);

    // gives the custom-valued vector (replace with Prof. vector)
    std::vector<double> customUpper(static_cast<int>(timeGridSize), 0.0);

    // returns the converted custom vector into a usable matrix
    return resUpper.toMatrix(customUpper);
}