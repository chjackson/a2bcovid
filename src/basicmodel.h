#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <vector>
using namespace std;

//#include <gsl/gsl_randist.h>
//#include <gsl/gsl_rng.h>
//#include <gsl/gsl_cdf.h>

struct run_params {
	double pa;	//Parameters for the infectious potential relative to time of symptoms
	double pb;
	double po;
	double smu;	//Parameters for the time of symptoms relative to date of infection
	double ssigma;
	double ucta; //Parameters for the distribution between becoming symptomatic and being tested
	double uctb; //For use when the time of becoming symptomatic is unknown
	double ucto;
	double uct_mean;  //Mean of the above distribution
	double rate;	//Substitutions per site per year
	double seq_noise; //Mean number of errors in sequencing
	double chat; //Prior probability of co-location
	int missing_spatial; //Strategy for missing spatial data.  0=remove datapoint.  1=constant probability.  2=Data-driven based on Addenbrookes data we have so far
	string ali_file;	//Name of alignment file
	string pat_file;	//Name of file with patient data i.e. code, date of symptoms, HCW status, etc.
	string mov_file;	//Name of file with HCW movement data
	string ward_file;	//Name of file with patient movement data
	string sub_file;	//Name of file with subset data
	int turner;	//Don't estimate uncertainty in tree
	int seqvar; //Flag to account for uncertainty over sequence composition
	int max_n; //Maximum number of ambiguous sites to allow per sequence.  Delete any sequence with more than this
	double min_qual; //Minimum sequence quality by coverage of the genome
	double t95NS;  //95% threshold, no sequence data
	double t99NS;  //99% threshold, no sequence data
	vector< vector<double> > threshold95; //95% thresholds by {D_A,D_B}
	vector< vector<double> > threshold99; //99% thresholds by {D_A,D_B}
	char pat_delim; //Delimiter for date information in the patient data file
	char mov_delim; //Delimiter for date information in the HCW movement data file
	int diagnostic; //Flag to run various diagnostics - prints out various data along the way
	int alt_like; //Flag for alternative likelihood model.  This increases the likelihood of transmission for each additional day of contact
	int calc_thresholds; //Calculate thresholds of the likelihood function
	double hcw_location_default;  //Probability a HCW is there given no other information
	double pat_location_default;  //Probability a HCW is there given no other information
	int ward_format_old; //Flag to use a previous formatting style for the ward_file data.  Redundant in this version; implemented in the R code
	int use_all_seqs; //Flag to use all sequences.  Don't discard duplicate sequences.  Find the maximum likelihood value across all sequences for each pair
	int error; //General error flag
	int symptom_uncertainty_calc;
};


struct loc {
	string ward;
	int date;
	double prob;
};

struct allele {
	int loc;
	char nuc;
	int count;
	double freq;
};

struct ipair {
	int loc;
	double q;
};

struct idat {
	int loc;
	double q;
	char allele;
};

struct tprob {
	int time;
	double weight;
};

struct ijlike {
	double lL_tot;
	double ns_lL_tot;
	int min; //Minimum time
	int max; //Maximum time
	vector<int> da; //Time of sequencing relative to s1.  Vector because the symptom time is a vector
	vector<int> db; //Time of sequencing relative to s2
	vector<int> contact_times;
	vector<double> contact_likes;
	vector<double> noseq_likes;
};

struct ts { //Structure to vectorise the time of becoming symptomatic
	vector<int> time;
	vector<double> prob;
	int most_likely;  //Use as a proxy for outputs
};

struct pat {
	string code;
	string seq;
	string code_match;
	int hcw; //Flag to indicate health care worker (as opposed to patient) 1 = HCW; 0 = Patient
	int type; //Flag for type of infection.  1 = Community.  2 = Patient.  3 = HCW
	ts time_s;  //Time of becoming symptomatic or time of testing positive
	int time_seq; //Time of collecting sample used in sequencing
	int dtime_s;  //If the time of becoming symptomatic is optimised, this gives the change in that time from the default
	int time_s_cert; //Certainty in time_s.  Flag 1 = Treat time_s as accurate.  0 = Asymptomatic or unknown.  time_s is the date of positive test
	int seq_n; //Number of N nucleotides in sequence
	vector<int> location; //Dates on which the individual is able to transmit or be infected due to their location
	vector<loc> locat; //Dates on which the individual is able to transmit or be infected due to their location
	vector<idat> seq_uncertainty; //Sites at which the sequence is not known
};

struct sparseseq {
    vector<int> locus;
    vector<char> allele;
};

struct branch {
	int from;
	int nbelow; //Number of individuals in clade
	vector<allele> vars;
	vector<int> indivs;
};


struct treestore {
	vector<sparseseq> variants;
	vector<string> individuals;
	vector<int> origins;
	vector<int> dests;
	vector<double> weights;
	double logL;
};

struct treestore_plus {
	vector<treestore> ts;
	int total_inc;
};

struct ed { //Directed from i to j
	int dest; //destination
	double weight;
	int origin; //Where this came from
};

struct treedatstore {
	int max_vertex;
	vector<int> best_edges;
	vector<ed> edgelist;
	vector< vector<double> > likelihoods;
	double logL_max;
	double sampling_LL;
};

struct sdist {
    int i;
    int j;
    int dist;
};

struct tpair {
	int from;
	int to;
};

struct tpairs {
	int root;
	vector<tpair> pair;
};

struct vertex {
	int index;
	vector<ed> edge;
};
