#include <cmath>
#include "bs.h"

// Black-Scholes object constructor
blackScholes::blackScholes(double spot, double strike, double rate, double divs, double repo, double vol, double matu)
{
    o_spot = spot;
    o_strike = strike;
    o_rate = rate;
    o_divs = divs;
    o_repo = repo;
    o_vol = vol;
    o_matu = matu;
}

double blackScholes::d1()
{
    double upper_1 = std::log(this->o_spot / this->o_strike);
    double upper_2 = (this->o_rate - this->o_divs - this->o_repo + 0.5 * this->o_vol * this->o_vol) * this->o_matu;
    double lower = this->o_vol * std::sqrt(this->o_matu);

    return (upper_1 + upper_2) / lower;
}

double blackScholes::d2()
{
    return d1() - this->o_vol * std::sqrt(this->o_matu);
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

// compute the theoretical price of a Black-Scholes call option
double blackScholes::callOptionPrice()
{
    double left = std::exp((-this->o_repo - this->o_divs) * this->o_matu) * this->o_spot * normCDF(d1());
    double right = std::exp(-this->o_rate * this->o_matu) * this->o_strike * normCDF(d2());

    return left - right;
}

// compute the theoretical price of a Black-Scholes put option
double blackScholes::putOptionPrice()
{
    double left = std::exp(-this->o_rate * this->o_matu) * this->o_strike * normCDF(-d2());
    double right = std::exp((-this->o_repo - this->o_divs) * this->o_matu) * this->o_spot * normCDF(-d1());

    return left - right;
}