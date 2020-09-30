#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <vector>

void RemoveDuplicatesNoSeq (vector<pat>& pdat);
void CorrectNames (vector<string>& names);
void MatchSequencetoPatient (vector<pat>& pdat);
void IncorporateSequenceData (run_params p, vector<pat>& pdat, vector<string>& names, vector<string>& seqs);
void CountNs (vector<pat>& pdat);

string FindMostCommonWard(vector<pat>& pdat);

void FindConsensus (string& consensus, vector<string>& seqs);

void FindConsensusPat (int& seq_len, string& consensus, vector<pat> pdat);
void FindSequenceLength(int& seq_len, string& consensus, vector<pat> pdat);


void FindVariants (vector<sparseseq>& variants, string& consensus, vector<pat>& pdat);
void FindPairwiseDistances (run_params p, vector< vector<int> >& seqdists, vector<sparseseq>& variants, vector<pat>& pdat);
void FindPairwiseDistancesNoSeq (run_params p, vector< vector<int> >& seqdists, vector<sparseseq>& variants, vector<pat>& pdat);

void FromConsensusDistances (const vector<sparseseq>& variants, vector< vector<tpair> >& seqdists_c);
void FromConsensusDistancesNoSeq (const vector<sparseseq>& variants, vector< vector<tpair> >& seqdists_c);

void ProcessSymptomUncertainty (run_params& p, vector<pat>& pdat);
void RemoveIndividualsNoSequence(const run_params p, vector<pat>& pdat, vector<string>& removed);
void RemoveIndividualsSequenceQuality(const run_params p, vector<pat>& pdat, vector<string>& removed);
void RemoveIndividualsMultipleN (const run_params p, vector<int>& nloc_count, vector<pat>& pdat, vector<sparseseq>& variants, vector<string>& removed);
void FixIndivudualsNoLocation (const run_params p, vector<pat>& pdat, vector<string>& fixed);
void RemoveRepeatPatients (run_params p, vector<pat>& pdat);
void RepairSequence(int i, int j, vector<pat>& pdat);
