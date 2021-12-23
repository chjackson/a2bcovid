#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <vector>


void GetOptions (run_params& p, int argc, const char **argv);
void GetThresholds (run_params p, vector< vector<double> >& thresholds95, vector< vector<double> >& thresholds99, double& t95NS, double& t99NS, int& error);
