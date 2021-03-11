#include <Rcpp.h>

using namespace Rcpp;

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>

#include "basicmodel.h"
#include "utilities.h"
#include "diagnostics.h"
#include "distributions.h"
#include "likelihoods.h"
#include "variants.h"
#include "thresholds.h"
#include "io.h"

run_params DefineParams(List params);

// [[Rcpp::export]]
DataFrame mainC(List params) {
  Rcpp::Function msg("message");
  
  run_params p;
  DataFrame out;

  p = DefineParams(params);

  vector<int> sdist_interval;
  vector<double> sdist_prob;
  PrecalculateSymptomDist(p,sdist_interval,sdist_prob);
	
  if (p.calc_thresholds==1) {
    CalculateThresholdsNoSeq(p);
    CalculateThresholdsFull(p);
//    return List::create();
  }


  //Code to read in CSV format
  vector<pat> pdat;
  if (p.ali_file.compare("")==0) {
	  ReadPatFromCSVNoSeq(p,pdat,sdist_interval,sdist_prob);  //No sequence data
    //Remove duplicate individuals by code
    RemoveDuplicatesNoSeq(pdat);
    if (pdat.size()<2) {
      Rcout << "Not enough individuals to assess transmission\n";
      return 0;
    }
  } else {
	  ReadPatFromCSV(p,pdat,sdist_interval,sdist_prob);
  }
  if (p.error==1) {
	return 0;
  }


  vector<string> names;
  vector<string> seqs;
  vector<string> removed;
  if (p.ali_file.compare("")!=0) {
    ReadFastaAli(p,names,seqs); //Read in large file containing genome sequences
	  Rcpp::Rcout << "Now here\n";
	CheckBaseCase(seqs);
    //Correct names >
	  Rcpp::Rcout << "Now here\n";
    CorrectNames(names);  //Convert the names from the fasta file into the format used for input
    //Add sequence data to records
	  Rcpp::Rcout << "Go to Incorporate\n";

    IncorporateSequenceData(p,pdat,names,seqs);
    RemoveIndividualsNoSequence(p,pdat,removed);
    //Calculate number of Ns in each sequence.  Use to discriminate between same day samples
    CountNs(pdat);
	if (p.use_all_seqs==0) {
      //Remove duplicate records for individuals with more than one sequence.  Take the earliest, best quality
	  RemoveRepeatPatients(p,pdat); //N.B. This currently takes the earliest regardless of quality.  Should impose a minimum quality...
	}
	RemoveIndividualsSequenceQuality(p,pdat,removed);
  }

  vector <vector<int> > repeatpatients;
  if (p.use_all_seqs==1) {
	FindRepeatPatients (p,pdat,repeatpatients);
  }

  // TODO sprintf this then msg
  std::ostringstream outstr; 
  outstr << "Have complete data for " << pdat.size() << " individuals";
  msg(outstr.str());

  //Code to read in location data
  int pd=0;
  if (p.ward_file.compare("")!=0) {
	  ReadWardMovFromCSV(p,pdat);
	  pd=1;
  }
  if (p.mov_file.compare("")!=0) {
	  ReadHCWMovFromCSV(p,pdat);
	  EditHCWMovData(pdat); //12 hour window of uncertainty - days with probability 0.5
	  pd=1;
  }
  if (p.error==1) {
    return 0;
  }

  if (p.diagnostic==1&&pd==1) {
	  PrintPdat(pdat);
  }

  //Deal with individuals with no sequence data - assumed present
  vector<string> fixed;
  FixIndivudualsNoLocation (p,pdat,fixed);
  if (p.diagnostic==1) {
    PrintPdat(pdat);
  }

  //Generate sequence variant data
  vector<sparseseq> variants;
  if (p.ali_file.compare("")==0) {
    FindVariantsNoSeq (variants,pdat);
  } else {
    //Convert sequences to variants
    string consensus;
    //Consensus of all sequences
    string all_consensus;
    FindConsensus(all_consensus,seqs);
    //Variants with respect to general consensus.
    FindVariants (variants,all_consensus,pdat);
    //Find sequences with an 'N' at sites with called variants
    vector<allele> allvar;
    ListAllVariantPositions(variants,allvar);

  outstr.str("");
  outstr << "Number of variants is " << variants.size(); 
  msg(outstr.str());

    if (p.diagnostic==1) {
	  PrintVariantLoci(allvar);
	  PrintVariants (variants,pdat);
    }
    vector<int> nloc_count;
    FindAmbiguousVarPositions (allvar,pdat,nloc_count);
    //Delete sequences with too many Ns
    RemoveIndividualsMultipleN(p,nloc_count,pdat,variants,removed);
  }

  //Find distances between sequences
  vector< vector<int> > seqdists;
  if (p.ali_file.compare("")==0) {
    FindPairwiseDistancesNoSeq (p,seqdists,variants,pdat);
  } else {
    FindPairwiseDistances (p,seqdists,variants,pdat);
  }

  //Find pairwise consensus distances
  vector< vector<tpair> > seqdists_c; //Distance to pairwise consensus
  if (p.ali_file.compare("")==0) {
    FromConsensusDistancesNoSeq (variants,seqdists_c);
  } else {
    FromConsensusDistances (variants,seqdists_c);
  }

  if (p.diagnostic==1) {
    PrintSequenceDistances(seqdists);
  }
  //Calculate pairwise likelihoods through time
  vector< vector<ijlike> > like_trans;
  CalculateTDLikelihoods (p,pdat,seqdists,seqdists_c,like_trans);

  vector<int> ordered;
  FindOrdering (like_trans,ordered);
  //Find best likelihoods for each indiviudal.  Redundant with --use_all_seqs 0
  FindBestLikelihoods (repeatpatients,like_trans,pdat);
  //Data output
  out = LikelihoodOutputR(p,ordered,pdat,like_trans);
  FinalOutput(p,removed,fixed);

  return out;
}


// Rcpp port of GetOptions in io.cpp
// "params" is created in R

run_params DefineParams(List params)
{
  run_params p;

  string p_switch;

  //Infectious potential relative to time of symptoms
  //Derived from Ashcroft with Bonhoeffer
  // TODO other comments
  p.pa = as<NumericVector>(params["pa"])[0];
  p.pb = as<NumericVector>(params["pb"])[0];
  p.po = as<NumericVector>(params["po"])[0];
  p.smu = as<NumericVector>(params["smu"])[0];
  p.ssigma = as<NumericVector>(params["ssigma"])[0];

  p.ucta = as<NumericVector>(params["ucta"])[0];
  p.uctb = as<NumericVector>(params["uctb"])[0];
  p.ucto = as<NumericVector>(params["ucto"])[0];
  p.uct_mean = as<NumericVector>(params["uct_mean"])[0];
  p.rate = as<NumericVector>(params["rate"])[0];
  p.seq_noise = as<NumericVector>(params["seq_noise"])[0];
  p.threshold = as<NumericVector>(params["threshold"])[0];
  p.thresholdns = as<NumericVector>(params["threshold_ns"])[0];
  p.max_n = as<NumericVector>(params["max_n"])[0];
  p.min_qual = as<NumericVector>(params["min_qual"])[0];

  p.ali_file = as<CharacterVector>(params["ali_file"])[0];
  p.pat_file = as<CharacterVector>(params["pat_file"])[0];
  p.mov_file = as<CharacterVector>(params["mov_file"])[0];
  p.ward_file = as<CharacterVector>(params["ward_file"])[0];

  p.pat_delim='/'; // Hard coded for now.
  p.mov_delim='.';

  p.diagnostic = as<NumericVector>(params["diagnostic"])[0];
  p.calc_thresholds = as<NumericVector>(params["calc_thresholds"])[0];
  p.hcw_location_default = as<NumericVector>(params["hcw_location_default"])[0];
  p.pat_location_default  = as<NumericVector>(params["pat_location_default"])[0];
  p.use_all_seqs  = as<NumericVector>(params["use_all_seqs"])[0];
  p.symptom_uncertainty_calc  = as<NumericVector>(params["symptom_uncertainty_calc"])[0];

  p.rate=(p.rate*29900)/365.25;
  SetThreshold(p);
  p.error = 0; 

  return p;
}
