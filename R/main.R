##' Sentence of what this function does here
##'
##' Longer description of what this function does here
##'
##' @param data_type Document this argument here
##'
##' @param pa Document this argument here
##'
##' @param pb Document this argument here
##'
##' @param po Etc...
##' @param smu  documentme
##' @param ssigma  documentme
##' @param ucta  documentme
##' @param uctb  documentme
##' @param ucto  documentme
##' @param uct_mean  documentme
##' @param rate  documentme
##' @param seq_noise  documentme
##' @param threshold  documentme
##' @param threshold_ns  documentme
##' @param max_n  documentme
##' @param min_qual  documentme
##' @param ali_file  documentme
##' @param pat_file  documentme
##' @param mov_file  documentme
##' @param ward_file  documentme
##' @param noseq  documentme
##' @param calc_thresholds  documentme
##' @param diagnostic  documentme
##' @param hcw_location_default  documentme
##' @param pat_location_default  documentme
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
