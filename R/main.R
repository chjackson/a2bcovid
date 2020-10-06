##' Sentence of what this function does here
##'
##' Longer description of what this function does here
##'
##' @param pat_file (Required) Specify the name of a file containing the basic data for each individual.  This should be a comma separated (.csv) file with data in columns: 1.  Individual ID (A code or identified corresponding to the individual)  2.  Onset date : The date at which the individual first experienced symptoms.  Date format should be dd/mm/yyyy.  3.  Onset date source : Equal to 1 if the date of onset is known.   Equal to 2 if the infection was asymptomatic.  In this case the onset date is the date on which the first positive swab was collected.  Equal to 3 if data is missing or unknown.  In this case the onset date is the date on which the first positive swab was collected.  If the onset date is anything other than 1 the true onset date is estimated by the code using data collected from Cambridge University hospitals. (Relevant command --uct_mean).  4.  Infection type : Equal to 1 if the individual is a patient and a community case (i.e. who could not have been infected by others in the dataset but who could potentially transmit the virus to others).  This was defined as being positive for the virus 48 hours before admission to hospital with no healthcare contact in the previous 14 days prior to admission).  Equal to 2 if the individual is a patient and not a community case (i.e. who could potentially transmit and receive infection).  Equal to 3 if the individual is a healthcare worker.  Whether or not an individual is a healthcare worker is set by this parameter.  5.  Sequence ID : A code used to link the individual to genome sequence information.  This should match the header of the sequence corresponding to the individual in the accompanying .fasta file (see --ali_file for this).  6.  Date of sample collection : Used in evolutionary calculations.  Date format should be dd/mm/yyyy.  7.  Sample received date : Currently not used in the calculation.
##' @param mov_file  Data describing when specific health care workers were on the ward in question.  The first line is a header line with column names.  The first two of these are labels, while those from the third column onwards describe dates, specified in dd.mm.yyyy format.  After the first line, the data is specified in columns as follows:  1.  Individual ID (same as for --pat_file)  2.  Cluster ID e.g. the name of the ward in question.  3 onwards.  Presence/absence data.  A 'Y' indicates that the health care worker was on the ward on the date specified for that column in the first row.  An 'N' indicates that the health care worker was not present on the ward on that date.  Either 'Y' or 'N' should be specified for each date.
##' @param ali_file  Specify the name of a file in FASTA format containing genome sequences.  This file must contain all required sequences, specified by the sequence ID in the data of pat_file.
##' @param ward_file  Specify the name of a file containing the location of patients over time.  This should be a comma separated (.csv) file.  The format of the file is designed to be compatible with local information on patient movements.  The first line is a header.  Subsequent lines are in columns as follows: 1.  Individual ID (same as for --pat_file)  2.  Cluster ID e.g. the name of the ward being studied.  3.  Infection type e.g. 'patient' or 'HCW' for health care worker.  4.  Availability of data e.g. 'patient_moves_available' 5 onwards.  Data of the location of a patient, in sets of three columns.  These specify in turn: i)   The name of the location of the individual e.g. WARD_01. ii)  The start date of the individual being in that location. iii) The end date of the individual being in that location.  In practice only the first column, and columns from 5 onwards are used.
##' @param data_type Specifes the kind of data which is being used to perform the calculation.  Options are 0: Only times of symptom onset.  No sequence data is needed.  The code will assume that all individuals are in the ward at the time.  1: Times of symptom onset and genome sequence data.  Again the code will assume that all individuals are in the ward throughout the period in question. (Relevant command : --ali_file) 2: As for 1, but with patient location data.  (Relevant command --ward_file) 3: As for 2, but with staff location data.  (Relevant command --mov_file)
##' @param evo_rate Rate of evolution of the virus, specified in nucleotide substitutions per locus per year.
##' @param seq_noise An estimate of the number of mutations separating two genome sequences that arises from sequencing noise.  The default parameter was estimated from data collected by Cambridge University Hospitals within single hosts, using the criteria that at least 90\% of the reported nucleotides were unambiguous.
##' @param min_qual Minimum sequence quality for a sequence to be included, measured as a fraction of genome coverage (e.g. 0.8 would indicate that at least 80\% of the genome must have been specified by a sequence
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
##' @param make_clusters Flag to sort output matrix into machine-derived clusters.  Alternative is that the input order of individuals is mirrored in the output
##'
##' @param threshold documentme
##' @param threshold_ns documentme
##' @param calc_thresholds documentme
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
a2bcovid <- function(
  pat_file,
  mov_file ,
  ali_file ,
  ward_file ,
  data_type = 2,
                  pa=97.18750, pb=0.268908, po=25.625,
                  smu=1.434065, ssigma=0.6612,
                  ucta=2.5932152095707406, uctb=3.7760060663975437, ucto=3.112080041460921,
                  uct_mean=6.67992,
                  evo_rate  = 0.0008,
                  seq_noise = 0.772469,
                  threshold=0, threshold_ns=0,
                  max_n = 10,
                  min_qual = 0.8,
                  noseq = 0,
                  calc_thresholds=FALSE,
                  diagnostic =FALSE,
                  hcw_default = 0.5714286,
                  pat_default = 1,
                  make_clusters = 1
)
{
  check_file(ali_file)
  check_file(pat_file)
  check_file(mov_file)
  check_file(ward_file)
  params <- list(data_type=data_type,
                 pa=pa, pb=pb, po=po, smu=smu, ssigma=ssigma,
                 ucta=ucta, uctb=uctb, ucto=ucto, uct_mean=uct_mean,
                 rate=evo_rate, seq_noise=seq_noise,
                 threshold=threshold, threshold_ns=threshold_ns,  max_n=max_n, min_qual=min_qual,
                 ali_file=ali_file, pat_file=pat_file,  mov_file=mov_file, ward_file=ward_file,
                 noseq=noseq,
                 calc_thresholds=calc_thresholds,
                 diagnostic=diagnostic,
				 hcw_location_default=hcw_default,
				 pat_location_default=pat_default,
				 make_clusters=make_clusters)
  res <- .Call(`_a2bcovid_mainC`, params)
  res$consistency <- ordered(res$consistency,
                             levels=c("Unlikely","Borderline","Consistent"))
  res
}

check_file <- function(filename){
  if (!file.exists(filename)) stop(sprintf("File `%s` not found", filename))
}


##' Plot results of an a2bcovid analysis
##'
##' Plots a grid of colours indicating likelihood of transmission paths between
##' each pair of individuals
##'
##'
##' @param x Data frame returned by \code{\link{a2bcovid}}.
##'
##' @param hi_from Character string, naming a variable in the dataframe
##'   indicating "from" individual IDs to be highlighted in the plot. If not
##'   supplied, then no IDs will be highlighted.
##'
##' @param hi_to Character string indicating "to" individual IDs to be
##'   highlighted, similarly.
##'
##' @param hi_col Colour to use to highlight individual IDs.
##'
##' @param palette Colour palette, passed to
##'   \code{\link[ggplot2]{scale_fill_brewer}}.
##'
##' @param direction Direction of colours. Defaults to 1.  Change to -1 to
##'   reverse the order of colours.
##'
##' @return A \pkg{ggplot2} plot object.
##'
##' @import ggplot2
##'
##' @export
plot_a2bcovid <- function(x, hi_from, hi_to, hi_col="red",
                          direction = 1){
  cbpallette3=c("#3C87C8","#FCF9DA","#D64E47")
  cbpallette21=c("#3C87C8","#FCF9DA")
  cbpallette22=c("#3C87C8","#D64E47")
  cbpallette23=c("#FCF9DA","#D64E47")
  cbpallette11=c("#3C87C8")
  cbpallette12=c("#FCF9DA")
  cbpallette13=c("#D64E47")
  col<-list()
  col[[1]]=cbpallette3
  col[[2]]=cbpallette21
  col[[3]]=cbpallette22
  col[[4]]=cbpallette23
  col[[5]]=cbpallette11
  col[[6]]=cbpallette12
  col[[7]]=cbpallette13
  x$from<-factor(x$from,unique(x$from))
  x$to<-factor(x$to,unique(x$to))
  if (!missing(hi_from)) {
    if (!(hi_from %in% names(x))) stop(sprintf("`%s` not found in `x`", hi_from))
    x_from <- x[!duplicated(x$from),]
    xhi_from <- x_from[match(levels(factor(x$from)), x_from$from), hi_from]
    x_cols <- ifelse(xhi_from, hi_col, "black")
  } else x_cols <- "black"
  if (!missing(hi_to)) {
    if (!(hi_to %in% names(x))) stop(sprintf("`%s` not found in `x`", hi_from))
    x_to <- x[!duplicated(x$to),]
    xhi_to <- x_to[match(levels(factor(x$to)), x_to$to), hi_to]
    y_cols <- ifelse(xhi_to, hi_col, "black")
  } else y_cols <- "black"
  cdat=unique((x$consistency))
  val=1
  if (length(cdat)==2&&cdat[1]=='Unlikely'&&cdat[2]=='Borderline'){
    val=2
  }
  if (length(cdat)==2&&cdat[1]=='Unlikely'&&cdat[2]=='Consistent'){
    val=3
  }
  if (length(cdat)==2&&cdat[1]=='Borderline'&&cdat[2]=='Consistent'){
    val=4
  }
  if (length(cdat)==1&&cdat[1]=='Unlikely'){
     val=5
  }
  if (length(cdat)==1&&cdat[1]=='Borderline'){
    val=6
  }
  if (length(cdat)==1&&cdat[1]=='Consistent'){
    val=7
  }
  ggplot(x, aes_string(y="from", x="to")) +
    geom_raster(aes_string(fill="consistency")) +
    theme(axis.text.x = ggtext::element_markdown(angle = 90, vjust=0.5, colour = x_cols),
          axis.text.y = ggtext::element_markdown(colour = y_cols),
          legend.title = element_blank()) +
    scale_fill_manual(values=col[[val]])
}
