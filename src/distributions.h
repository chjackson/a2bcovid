#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <vector>

double OffsetGammaCDFFlex (double x, double a, double b, double o);
double LogNormal (double x, double mu, double sigma);
double Poisson (int n, double lambda);
void PrecalculateSymptomDist (run_params p, vector<int>& sdist_interval, vector<double>& sdist_prob);

