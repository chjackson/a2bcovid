#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <vector>

struct stat {
	int s1;
	int s2;
	int t;
	int d1;
	int d2;
	int h;
	int b;
	int h1;
	int h2;
	double L;
};

struct cstat {
	int s1;
	int s2;
	int t;
	int d1;
	int d2;
	int h;
	int b;
	int h1;
	int h2;
	double L;
};

void CalculateThresholdsNoSeq (run_params p);
void CalculateThresholdsFull (run_params p);
void ReCalculateThresholds (run_params p);


