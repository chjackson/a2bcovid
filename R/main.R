##' Sentence of what this function does here
##'
##' Longer description of what this function does here
##'
##' @param data_type Specifes the kind of data which is being used to perform the calculation.  Options are 0: Only times of symptom onset.  No sequence data is needed.  The code will assume that all individuals are in the ward at the time.  1: Times of symptom onset and genome sequence data.  Again the code will assume that all individuals are in the ward throughout the period in question. (Relevant command : --ali_file) 2: As for 1, but with patient location data.  (Relevant command --ward_file) 3: As for 2, but with staff location data.  (Relevant command --mov_file)
##' @param pat_file (Required) Specify the name of a file containing the basic data for each individual.  This should be a comma separated (.csv) file with data in columns: 1.  Individual ID (A code or identified corresponding to the individual)  2.  Onset date : The date at which the individual first experienced symptoms.  Date format should be dd/mm/yyyy.  3.  Onset date source : Equal to 1 if the date of onset is known.   Equal to 2 if the infection was asymptomatic.  In this case the onset date is the date on which the first positive swab was collected.  Equal to 3 if data is missing or unknown.  In this case the onset date is the date on which the first positive swab was collected.  If the onset date is anything other than 1 the true onset date is estimated by the code using data collected from Cambridge University hospitals. (Relevant command --uct_mean).  4.  Infection type : Equal to 1 if the individual is a patient and a community case (i.e. who could not have been infected by others in the dataset but who could potentially transmit the virus to others).  This was defined as being positive for the virus 48 hours before admission to hospital with no healthcare contact in the previous 14 days prior to admission).  Equal to 2 if the individual is a patient and not a community case (i.e. who could potentially transmit and receive infection).  Equal to 3 if the individual is a healthcare worker.  5.  Sequence ID : A code used to link the individual to genome sequence information.  This should match the header of the sequence corresponding to the individual in the accompanying .fasta file (see --ali_file for this).  6.  Date of sample collection : Used in evolutionary calculations.  Date format should be dd/mm/yyyy.  7.  Sample received date : Currently not used in the calculation.
##' @param ali_file  Specify the name of a file in FASTA format containing genome sequences.  This file must contain all required sequences, specified by the sequence ID in the data of pat_file.
##' @param ward_file  Specify the name of a file containing the location of patients over time.  This should be a comma separated (.csv) file.  The format of the file is designed to be compatible with local information on patient movements.  The first line is a header.  Subsequent lines are in columns as follows: 1.  Individual ID (same as for --pat_file)  2.  Cluster ID e.g. the name of the ward being studied.  3.  Infection type e.g. 'patient' or 'HCW' for health care worker.  4.  Availability of data e.g. 'patient_moves_available' 5 onwards.  Data of the location of a patient, in sets of three columns.  These specify in turn: i)   The name of the location of the individual e.g. WARD_01. ii)  The start date of the individual being in that location. iii) The end date of the individual being in that location.  In practice only the first column, and columns from 5 onwards are used.
##' @param mov_file  Data describing when specific health care workers were on the ward in question.  The first line is a header line with column names.  The first two of these are labels, while those from the third column onwards describe dates, specified in dd.mm.yyyy format.  After the first line, the data is specified in columns as follows:  1.  Individual ID (same as for --pat_file)  2.  Cluster ID e.g. the name of the ward in question.  3 onwards.  Presence/absence data.  A 'Y' indicates that the health care worker was on the ward on the date specified for that column in the first row.  An 'N' indicates that the health care worker was not present on the ward on that date.  Either 'Y' or 'N' should be specified for each date.
##' @param evo_rate Rate of evolution of the virus, specified in nucleotide substitutions per locus per year.
##' @param seq_noise An estimate of the number of mutations separating two genome sequences that arises from sequencing noise.  The default parameter was estimated from data collected by Cambridge University Hospitals within single hosts, using the criteria that at least 90% of the reported nucleotides were unambiguous.
##' @param min_qual Minimum sequence quality for a sequence to be included, measured as a fraction of genome coverage (e.g. 0.8 would indicate that at least 80% of the genome must have been specified by a sequence
##' @param max_n Maximum number of ambiguous nucleotides tolerated in a sequence counted at positions in the sequence data for which there is a polymorphism.  This parameter deals with a case of a sequence of generally high quality in which the missing coverage of the genome is all at critical sites
##' @param pat_default Default probability of a patient being present on the ward on a given day if no location information is specified for that individual.  Default 1.
##' @param hcw_default Default probability of a health care worker being present on the ward on a given day if no location information is specified for that individual.  Default is 4/7.
##' @param uct_mean Mean time between an individual becoming symptomatic for coronavirus infection and testing positive.  This value is used to estimate times of individuals becoming symptomatic in the case that no symptom dates are available
##' @param ucta Alpha parameter for a gamma distribution of the times bewterrn becomining symptomatic and testing positive.  Currently not used.
##' @param uctb Beta parameter for a gamma distribution of the times bewterrn becomining symptomatic and testing positive.  Currently not used.
##' @param ucto Offset parameter for a gamma distribution of the times bewterrn becomining symptomatic and testing positive.  Currently not used.
##' @param pa Alpha parameter for the gamma distribution for the infectious potential of an individual.
##' @param pb Beta parameter for the gamma distribution for the infectious potential of an individual.
##' @param po Offset parameter for the gamma distribution for the infectious potential of an individual.
##' @param smu Mu parameter for the lognormal distribution for the time from infection to becoming symptomatic.
##' @param ssigma Sigma parameter for the lognormal distribution for the time from infection to becoming symptomatic.
##' @param diagnostic Flag to enable extensive diagnostic output from the function.
##' @param noseq  documentme
##'
##'
##' @return Document the return value here
##'
##' @author Our names
##'
##' @references Our paper
##'
##' @importFrom Rcpp evalCpp
##' @importFrom shiny runApp
##' @importFrom rstudioapi viewer
##' @useDynLib a2bcovid, .registration = TRUE
##'
##' @export
mainR <- function(data_type = 2,
                  pa=97.18750, pb=0.268908, po=25.625,
                  smu=1.434065, ssigma=0.6612,
                  ucta=2.5932152095707406, uctb=3.7760060663975437, ucto=3.112080041460921,
                  uct_mean=6.67992,
                  rate  = 0.0008,
                  seq_noise = 0.772469,
                  threshold=0, threshold_ns=0,
                  max_n = 10,
                  min_qual = 0.8,
                  ali_file = "Seqs_editN20_manali_plus.fa",
                  pat_file = "data1.csv",
                  mov_file = "data2.csv",
                  ward_file = "ward_movement_network_edit_anonymised_20200811_NoPII.csv",
                  noseq = 0,
                  calc_thresholds=FALSE,
                  diagnostic =FALSE,
                  hcw_location_default = 0.5714286,
                  pat_location_default = 1
)
{
  params <- list(data_type=data_type,
                 pa=pa, pb=pb, po=po, smu=smu, ssigma=ssigma,
                 ucta=ucta, uctb=uctb, ucto=ucto, uct_mean=uct_mean,
                 rate=rate, seq_noise=seq_noise,
                 threshold=threshold, threshold_ns=threshold_ns,  max_n=max_n, min_qual=min_qual,
                 ali_file=ali_file, pat_file=pat_file,  mov_file=mov_file, ward_file=ward_file,
                 noseq=noseq,
                 calc_thresholds=calc_thresholds,
                 diagnostic=diagnostic,
                 hcw_location_default=hcw_location_default, pat_location_default=pat_location_default)
  .Call(`_a2bcovid_mainC`, params);
}
