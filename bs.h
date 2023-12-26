#include <cstddef>
#include <vector>

class blackScholes
{
public:
    blackScholes(double spot, double strike, double rate, double divs, double vol, double matu);

    // Methods to calculate option prices
    double callOptionPrice();
    double putOptionPrice();

private:
    // Private methods for internal calculations
    double d1();
    double d2();

    // Standard normal cumulative distribution function
    double normCDF(double value);

    // Member variables
    double spot;
    double strike;
    double rate;
    double divs;
    double vol;
    double matu;
};