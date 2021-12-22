using namespace std;
#include "basicmodel.h"
#include "likelihoods.h"

void ThresholdComparison (run_params& p, int sa, int sb, int da, int db, double lL) {
	double t95=0;
	double t99=0;
	int da_temp=da-sa+10;
	int db_temp=db-sb+10;
	if (da_temp<0||db_temp<0||da_temp>51||db_temp>51) {
		cout << "Unlikely\n";
	} else {
		t95=p.threshold95[da_temp][db_temp];
		t99=p.threshold99[da_temp][db_temp];
		if (lL>t95) {
			//cout << "Consistent " << t95 << " " << t99 << "\n";
			cout << "Consistent\n";
		} else if (lL>t99) {
			cout << "Borderline\n";
			//cout << "Borderline " << t95 << " " << t99 << "\n";
		} else {
			cout << "Unlikely\n";
			//cout << "Unlikely "  << t95 << " " << t99 << "\n";
		}
	}
}
