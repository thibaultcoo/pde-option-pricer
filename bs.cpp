#include <iostream>
#include <vector>
#include <cmath>
#include "bs.h"

// black-Scholes object constructor
blackScholes::blackScholes(double spotOpt, double strikeOpt, double rateOpt, double divsOpt, double volOpt, double matuOpt)
{
    double spot = spotOpt;
    double strike = strikeOpt;
    double rate = rateOpt;
    double divs = divsOpt;
    double vol = volOpt;
    double matu = matuOpt;
}

double blackScholes::d1()
{
    double upper = std::log(this->spot / this->strike) + (this->rate - this->divs + 0.5 * this->vol * this->vol) * this->matu;
    double lower = this->vol * std::sqrt(this->matu);

    return upper / lower;
}

double blackScholes::d2()
{
    return d1() - this->vol * std::sqrt(this->matu);
}

// approximation by Abramowitz and Stegun; Handbook of Mathematical Functions
double blackScholes::normCDF(double value)
{
    const double a1 = 0.254829592;
    const double a2 = -0.284496736;
    const double a3 = 1.421413741;
    const double a4 = -1.453152027;
    const double a5 = 1.061405429;
    const double p = 0.3275911;

    // Save the sign of the value
    int sign = 1;
    if (value < 0) {sign = -1;}
    value = std::fabs(value) / std::sqrt(2.0);

    double t = 1.0 / (1.0 + p * value);
    double y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * std::exp(-value * value);

    return 0.5 * (1.0 + sign * y);
}

// computest the theoretical price of a Black-Scholes call option
double blackScholes::callOptionPrice()
{
    double left = std::exp(this->divs * this->matu) * this->spot * normCDF(d1());
    double right = std::exp(-this->rate * this->matu) * this->strike * normCDF(d2());

    return left - right;
}

// computest the theoretical price of a Black-Scholes put option
double blackScholes::callOptionPrice()
{
    double left = std::exp(-this->rate * this->matu) * this->strike * normCDF(-d2());
    double right = std::exp(this->divs * this->matu) * this->spot * normCDF(-d1());

    return left - right;
}