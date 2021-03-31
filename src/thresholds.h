#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <vector>

void CalculateThresholdsFull (run_params p, const vector<double>& OGPreCalcP, const vector<double>& LNPreCalc);
void CalculateThresholdsFullExplicit (run_params p, const vector<double>& OGPreCalcP, const vector<double>& LNPreCalc);
vector<double> CalculateThresholdsSpecificDMean (run_params p, int d1, int d2d, const vector<double>& OGPreCalcP, const vector<double>& LNPreCalc);
void CalculateThresholdsNoSeq (run_params p, const vector<double>& OGPreCalcP, const vector<double>& LNPreCalc);
vector<double> CalculateThresholdsNSSpecific (run_params p, const vector<double>& OGPreCalcP, const vector<double>& LNPreCalc);
vector<double> CalculateThresholdsSpecificDExplicit (run_params p, int d1, int d2d, const vector<double>& OGPreCalcP, const vector<double>& LNPreCalc);


