#include <cstddef>
#include <vector>

class blackScholes
{
public:
    blackScholes(double spot, double strike, double rate, double divs, double vol, double matu);
    double callOptionPrice();
    double putOptionPrice();

private:
    double d1();
    double d2();
    double normCDF(double value);

    double o_spot;
    double o_strike;
    double o_rate;
    double o_divs;
    double o_vol;
    double o_matu;
};