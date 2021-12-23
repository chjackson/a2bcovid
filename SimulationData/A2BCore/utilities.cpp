#include <vector>
#include <cmath>
#include <iostream>
#include <sstream>
#include <algorithm>

#include "basicmodel.h"
#include "utilities.h"
#include "distributions.h"
#include "likelihoods.h"
#include "io.h"


void GetPossibleTimes (int sa, int sb, vector<int>& possible_tab) {
	int ms=min(sa,sb)-15;
	int mx=max(sa,sb)+25;
	for (int d=ms;d<=mx;d++) {
		possible_tab.push_back(d);
	}
}
