#include "basicmodel.h"
#include "likelihoods.h"
#include "io.h"
#include <string>

void GetOptions (run_params& p, int argc, const char **argv) {
	string p_switch;
	//Infectious potential relative to time of symptoms
	//Derived from Ashcroft with Bonhoeffer
	p.pa=97.18750;
	p.pb=0.268908;
	p.po=25.625;
	
	//Time of symptoms relative to date of infection
	//Lognormal distribution cited by He et al.
	p.smu=1.434065;
	p.ssigma=0.6612;
	
	p.delta=0;
	
	//The following are learnt from CUH data.
	//They characterise the distribution of time from positive test to date of reporting symptoms
	p.ucta=2.5932152095707406;
	p.uctb=3.7760060663975437;
	p.ucto=3.112080041460921;
	p.uct_mean=6.67992;
	
	//Can create an option to optimise over the unknown times of becoming symptomatic
	p.rate=0.0008;
	p.seq_noise=0.41369;
	p.chat=0.5;
	p.max_n=10;
	p.min_qual=0.8;
    p.ali_file="NULL";  //Genome sequence alignment
    p.pat_file="NULL";  //Basic patient data e.g. dates of symptom onset
    p.mov_file="NULL";  //HCW movements
	p.ward_file="NULL";  //Patient movements i.e. ward data
	p.pat_delim='/';
	p.mov_delim='.';
	p.diagnostic=0;
	p.calc_thresholds=0;
	p.hcw_location_default=0.5714286;
	p.pat_location_default=1;
	p.make_clusters=1;
	p.error=0;
	p.ward_format_old=0;
	p.use_all_seqs=0;
	p.symptom_uncertainty_calc=0;
	p.subset=0;
	int x=1;
	while (x < argc && (argv[x][0]=='-')) {
		p_switch=argv[x];
		if (p_switch.compare("--evo_rate")==0) {
			x++;
			p.rate=atof(argv[x]);
		} else if (p_switch.compare("--seq_noise")==0) {
			x++;
			p.seq_noise=atof(argv[x]);
		} else if (p_switch.compare("--chat")==0) {
			x++;
			p.chat=atof(argv[x]);
		} else if (p_switch.compare("--diag")==0) {
			x++;
			p.diagnostic=atoi(argv[x]);
		} else if (p_switch.compare("--ali_file")==0) {
			x++;
			p.ali_file=argv[x];
		} else if (p_switch.compare("--pat_file")==0) {
			x++;
			p.pat_file=argv[x];
		} else if (p_switch.compare("--mov_file")==0) {
			x++;
			p.mov_file=argv[x];
		} else if (p_switch.compare("--ward_file")==0) {
			x++;
			p.ward_file=argv[x];
		} else if (p_switch.compare("--sub_file")==0) {
			x++;
			p.sub_file=argv[x];
		} else if (p_switch.compare("--hcw_loc_default")==0) {
			x++;
			p.hcw_location_default=atof(argv[x]);
		} else if (p_switch.compare("--pat_loc_default")==0) {
			x++;
			p.pat_location_default=atoi(argv[x]);
		} else if (p_switch.compare("--maxn")==0) {
			x++;
			p.max_n=atoi(argv[x]);
		} else if (p_switch.compare("--delta")==0) {
			x++;
			p.delta=atof(argv[x]);
		} else if (p_switch.compare("--min_qual")==0) {
			x++;
			p.min_qual=atof(argv[x]);
		} else if (p_switch.compare("--calc_thresholds")==0) {
			x++;
			p.calc_thresholds=atoi(argv[x]);
		} else if (p_switch.compare("--use_all_seqs")==0) {
			x++;
			p.use_all_seqs=atoi(argv[x]);
		} else if (p_switch.compare("--symptom_uncertainty_calc")==0) {
			x++;
			p.symptom_uncertainty_calc=atoi(argv[x]);
		} else if (p_switch.compare("--ward_format_old")==0) {
			x++;
			p.ward_format_old=atoi(argv[x]);
		} else if (p_switch.compare("--subset")==0) {
			x++;
			p.subset=atoi(argv[x]);
//		} else if (p_switch.compare("--cluster")==0) {
//			x++;
//			p.make_clusters=atoi(argv[x]);
        } else {
			cout << "Incorrect usage " << argv[x] << "\n ";
			exit(1);
		}
		p_switch.clear();
		x++;
	}
	p.rate=(p.rate*29900)/365.25;
	/*DELTA PARAMETERS*/
	if (p.delta==1) {
		p.pa=38.4805;
		p.pb=0.468049;
		p.po=20;
		p.smu=1.39599;
		p.ssigma=0.41354;
	}
}

void GetThresholds (run_params p, vector< vector<double> >& thresholds95, vector< vector<double> >& thresholds99, double& t95NS, double& t99NS, int& error) {
	t95NS=0;
	t99NS=0;
	ifstream t95;
	ifstream t99;
	ifstream t95n;
	ifstream t99n;
	if (p.delta==0) {
		t95.open("../Data/Thresholds95.dat");//With sequences
		t99.open("../Data/Thresholds99.dat");
		t95n.open("../Data/Thresholds95NS.dat");//No sequences
		t99n.open("../Data/Thresholds99NS.dat");
	} else {
		t95.open("../Data/ThresholdsD95.dat");//With sequences
		t99.open("../Data/ThresholdsD99.dat");
		t95n.open("../Data/ThresholdsD95NS.dat");//No sequences
		t99n.open("../Data/ThresholdsD99NS.dat");
	}
	int n;
	double x;
	t95n >> x;
	t95NS=x;
	t99n >> x;
	t99NS=x;
	if (t95NS>-0.01) {
		if (p.delta==1) {
			cout << "Error reading NS threshold: File not found ../Data/ThresholdsD95NS.dat\n";
		} else {
			cout << "Error reading NS threshold: File not found ../Data/Thresholds95NS.dat\n";
		}
		cout << "Note: The missing file can be generated using the --calc_thresholds 1 flag\n";
		error=1;
	}
	if (t99NS>-0.01) {
		if (p.delta==1) {
			cout << "Error reading NS threshold: File not found ../Data/ThresholdsD99NS.dat\n";
		} else {
			cout << "Error reading NS threshold: File not found ../Data/Thresholds99NS.dat\n";
		}
		cout << "Note: The missing file can be generated using the --calc_thresholds 1 flag\n";
		error=1;
	}
	int index1=-10;
	vector<double> t;
	for (int i=0;i<1000000;i++) {
		if (!(t95 >> n)) break;
		if (!(t95 >> n)) break;
		if (!(t95 >> x)) break;
		t.push_back(x);
		index1++;
		if (index1>40) {
			index1=-10;
			thresholds95.push_back(t);
			t.clear();
		}
	}
	index1=-10;
	for (int i=0;i<1000000;i++) {
		if (!(t99 >> n)) break;
		if (!(t99 >> n)) break;
		if (!(t99 >> x)) break;
		t.push_back(x);
		index1++;
		if (index1>40) {
			index1=-10;
			thresholds99.push_back(t);
			t.clear();
		}
	}
	if (thresholds95.size()<1) {
		if (p.delta==1) {
			cout << "Error reading thresholds.  Can't find file ../Data/ThresholdsD95.dat: Check this\n";
		} else {
			cout << "Error reading thresholds.  Can't find file ../Data/Thresholds95.dat: Check this\n";
		}
		cout << "Note: The missing file can be generated using the --calc_thresholds 1 flag\n";
		error=1;
	}
	if (thresholds99.size()<1) {
		if (p.delta==1) {
			cout << "Error reading thresholds.  Can't find file ../Data/ThresholdsD99.dat: Check this\n";
		} else {
			cout << "Error reading thresholds.  Can't find file ../Data/Thresholds99.dat: Check this\n";
		}
		cout << "Note: The missing file can be generated using the --calc_thresholds 1 flag\n";
		error=1;
	}

}

