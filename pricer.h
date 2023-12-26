#include <vector>

class pricerPDE 
{
public:
    pricerPDE(double strike, double matu, double vol, 
              double rate, double divs, double repo,
              double multiplier, int m, int n);

    void setupGrid();
    void initializeBoundaryConditions();
    void applyTerminalConditions();
    void solvePDE();
    double getPrice();

private:
    std::vector<std::vector<double>> p_priceGrid;
    std::vector<std::vector<double>> p_boundaryConditions;
    std::vector<double> p_terminalCondition;
    std::vector<double> p_timeGrid;
    std::vector<double> p_spotGrid;

    double p_multiplier;
    double p_strike;
    double p_matu;
    double p_vol;
    double p_rate;
    double p_divs;
    double p_repo;
    int p_m;
    int p_n;
};
