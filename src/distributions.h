#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <vector>

void PreCalculateLikelihoods (run_params p, vector<double>& LNPreCalc, vector<double>& OGPreCalcP);
void MakeLogNormalPreCalc (run_params p, vector<double>& LNPreCalc);
double LogNormal (double x, double mu, double sigma);
double LogNormalP (double x, const vector<double>& LNPreCalc);
void MakeOffsetGammaPreCalcP (run_params p, vector<double>& OGPreCalcP);
double OffsetGammaCDFFlex (double x, double a, double b, double o);
double OffsetGammaPreCalcP (int x, const vector<double>& OGPreCalcP);
double Poisson (int n, double lambda);
void PrecalculateSymptomDist (run_params p, vector<int>& sdist_interval, vector<double>& sdist_prob);

