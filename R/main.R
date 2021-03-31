##' Identifying clusters of COVID-19 infections from genome sequence and
##' individual location data
##'
##' A tool which combines genome sequence and the locations of infected
##' individuals, using a statistical and evolutionary model, to estimate the
##' likelihood that transmission occurred between particular individuals, and
##' then to identify clusters of infections.  It is currently designed to apply
##' to COVID-19 infection dynamics on hospital wards.
##'
##'
##' @param pat_file (Required) A character string with the path to a file
##'   containing the basic data for each individual. This should be a comma
##'   separated (.csv) file with data in columns:
##'
##'   1. Individual ID (A code or identifier corresponding to the individual)
##'
##'   2. Onset date.  The date at which the individual first experienced
##'   symptoms. Date format should be dd/mm/yyyy.
##'
##'   3.  Onset date source : Equal to 1 if the date of onset is known.   Equal
##'   to 2 if the infection was asymptomatic.  In this case the onset date is
##'   the date on which the first positive swab was collected. Equal to 3 if
##'   data is missing or unknown.  In this case the onset date is the date on
##'   which the first positive swab was collected.  If the onset date is
##'   anything other than 1 the true onset date is estimated by the code using
##'   data collected from Cambridge University hospitals (argument
##'   \code{uct_mean}).
##'
##'   4.  Infection type : Equal to 1 if the individual is a patient and a
##'   community case (i.e. who could not have been infected by others in the
##'   dataset but who could potentially transmit the virus to others).  This was
##'   defined as being positive for the virus 48 hours before admission to
##'   hospital with no healthcare contact in the previous 14 days prior to
##'   admission).  Equal to 2 if the individual is a patient and not a community
##'   case (i.e. who could potentially transmit and receive infection).  Equal
##'   to 3 if the individual is a healthcare worker.  Whether or not an
##'   individual is a healthcare worker is set by this parameter.
##'
##'   5. Sequence ID : A code used to link the individual to genome sequence
##'   information.  This should match the header of the sequence corresponding
##'   to the individual in the accompanying .fasta file (arcument
##'   \code{ali_file}).
##'
##'   6.  Date of sample collection : Used in evolutionary calculations. Date
##'   format should be dd/mm/yyyy.
##'
##'   7.  Sample received date : Currently not used in the calculation.
##'
##'   An example is given with the installed package.  The path to the example
##'   file can be shown by the R command \code{system.file("extdata",
##'   "Example_genetic_temporal_data.csv", package="a2bcovid")}
##'
##'
##'
##' @param hcw_loc_file  A character string with the path to a file of data
##'   describing when specific health care workers were on the ward in question.
##'   If this argument is omitted or set to an empty string, then this kind of
##'   data is not used in the calculation.
##'
##'   The first line is a header line with column names. The first two of these
##'   are labels, while those from the third column onwards describe dates,
##'   specified in dd.mm.yyyy format.  After the first line, the data is
##'   specified in columns as follows:
##'
##'   1.  Individual ID (same as for \code{pat_file}).
##'
##'   2.  Cluster ID e.g. the name of the ward in question.
##'
##'   3 onwards:  Presence/absence data.  A 'Y' indicates that the health care
##'   worker was on the ward on the date specified for that column in the first
##'   row.  An 'N' indicates that the health care worker was not present on the
##'   ward on that date.  Either 'Y' or 'N' should be specified for each date.
##'
##'   An example is given with the installed package.  The path to the example
##'   file can be shown by the R command \code{system.file("extdata",
##'   "Example_movement_file.csv", package="a2bcovid")}
##'
##' @param ali_file  A character string with the path to a file in FASTA format
##'   containing genome sequence alignments.  This file must contain all required
##'   sequences, specified by the sequence ID in the data of \code{pat_file}. If
##' this argument is omitted or set to an empty string, then genomic data are
##' not used in the calculation.
##'
##' An example is given with the installed package.  The path to the example
##' file can be shown by the R command \code{system.file("extdata",
##' "Example_sequences.fa", package="a2bcovid")}
##'
##' @param pat_loc_file  A character string with the path to a file containing
##'   the location of patients over time. If this argument is omitted or set to
##'   an empty string, then this kind of data is not used in the calculation.
##'
##'   This should be a comma separated (.csv) file.  Two alternative formats are
##'   accepted, "wide" and "long" formats.    These are based on the formats in
##'   use in the hospital setting where the package was developed.
##'
##'   The first line of the file should be a header with variable names.  If
##'   there is a variable called \code{"start_date"} then the file is assumed to
##'   be in long format.  Or if there is a variable called \code{"StartDate_0"}
##'   then the file is assumed to be in wide format.  Otherwise the variable
##'   names are ignored.   The columns should appear in the specified order.
##'
##'   The names don't matter, but the columns should appear in the specified
##'   order.
##'
##'   In wide format, each row represents a different  patient.
##'   The file should have the following columns:
##'
##'   1. Individual ID (same as for \code{pat_file}).
##'
##'   2. Cluster ID e.g. the name of the ward being studied.
##'
##'   3. Infection type e.g. 'patient' or 'HCW' for health care worker.
##'
##'   4. Availability of data e.g. 'patient_moves_available'.
##'
##'   5 onwards.  Data of the location of a patient, in sets of three columns.
##'   These specify in turn: i)   The name of the location of the individual
##'   e.g. WARD_01. ii)  The start date of the individual being in that
##'   location. iii) The end date of the individual being in that location.  In
##'   practice only the first column, and columns from 5 onwards are used.

##'   An example is given with the installed package.  The path to the example
##'   file can be shown by the R command \code{system.file("extdata",
##'   "Example_pat_loc_file_wide.csv", package="a2bcovid")}
##'
##'
##'   In long format, each row represents a single
##'   stay on a specific ward for a specific patient.  The file should have the
##'   following columns:
##'
##'   1. Individual ID
##'
##'   2. Cluster ID, typically the name of the ward
##'
##'   3. Start date/time for the ward stay, in d/m/Y format (optionally in d/m/Y
##'   H:M format, but the time is currently ignored)
##'
##'   4. Name of the ward the patient went to next (or "Discharge") if they were
##'   discharged
##'
##'   5. End date/time for the ward stay, in d/m/Y format (optionally in d/m/Y
##'   H:M format, but the time is currently ignored).
##'
##'   An example is given with the installed package.  The path to the example
##'   file can be shown by the R command \code{system.file("extdata",
##'   "Example_pat_loc_file_wide.csv", package="a2bcovid")}
##'
##' @param evo_rate Rate of evolution of the virus, specified in nucleotide
##'   substitutions per locus per year.
##'
##' @param seq_noise An estimate of the number of mutations separating two
##'   genome sequences that arises from sequencing noise.  The default parameter
##'   was estimated from data collected by Cambridge University Hospitals within
##'   single hosts, using the criteria that at least 90\% of the reported
##'   nucleotides were unambiguous.
##'
##' @param chat Prior estimate of the probability of any two individuals being
##'   in contact on any given day, conditional on transmission between the two
##'   individuals having taken place.
##'
##' @param min_qual Minimum sequence quality for a sequence to be included,
##'   measured as a fraction of genome coverage (e.g. 0.8 would indicate that at
##'   least 80\% of the genome must have been specified by a sequence
##'
##' @param max_n Maximum number of ambiguous nucleotides tolerated in a sequence
##'   counted at positions in the sequence data for which there is a
##'   polymorphism.  This parameter deals with a case of a sequence of generally
##'   high quality in which the missing coverage of the genome is all at
##'   critical sites
##'
##' @param pat_default Default probability of a patient being present on the
##'   ward on a given day if no location information is specified for that
##'   individual.  Default 1.
##'
##' @param hcw_default Default probability of a health care worker being present
##'   on the ward on a given day if no location information is specified for
##'   that individual.  Default is 4/7.
##'
##' @param uct_mean Mean time between an individual becoming symptomatic for
##'   coronavirus infection and testing positive.  This value is used to
##'   estimate times of individuals becoming symptomatic in the case that no
##'   symptom dates are available
##'
##' @param ucta Alpha parameter for a gamma distribution of the times bewterrn
##'   becomining symptomatic and testing positive.  Currently not used.
##'
##' @param uctb Beta parameter for a gamma distribution of the times bewterrn
##'   becomining symptomatic and testing positive.  Currently not used.
##'
##' @param ucto Offset parameter for a gamma distribution of the times bewterrn
##'   becomining symptomatic and testing positive.  Currently not used.
##'
##' @param pa Alpha parameter for the gamma distribution for the infectious
##'   potential of an individual.
##'
##' @param pb Beta parameter for the gamma distribution for the infectious
##'   potential of an individual.
##'
##' @param po Offset parameter for the gamma distribution for the infectious
##'   potential of an individual.
##'
##' @param smu Mu parameter for the lognormal distribution for the time from
##'   infection to becoming symptomatic.
##'
##' @param ssigma Sigma parameter for the lognormal distribution for the time
##'   from infection to becoming symptomatic.
##'
##' @param diagnostic Flag to enable extensive diagnostic output from the
##'   function.
##'
##' @param use_all_seqs Flag to use multiple sequences from an individual, rather
##'   than simply the first collected.  Reports the maximum likelihood calculated
##'   across all sequences from an individual.
##'
##' @symptom_uncertainty_calc Flag to use a complete offset gamma distribution,
##' specified by the parameters ucta, uctb, and ucto, to model the undertainty
##' in the date of onset of symptom
##'
##' @param calc_thresholds documentme
##'
##'
##' @return A data frame with the following columns
##'
##'   \code{from}
##'
##'   \code{to}
##'
##'   \code{hcw_from}
##'
##'   \code{hcw_to}
##'
##'   \code{ordered_i}
##'
##'   \code{ordered_j}
##'
##'   \code{likelihood}
##'
##'   \code{consistency}
##'
##'   \code{under_threshold}
##'
##' @author Chris Illingworth \email{chris.illingworth@mrc.bsu.cam.ac.uk}, Chris Jackson
##' \email{chris.jackson@mrc.bsu.cam.ac.uk}.
##'
##' @references "A2B-Covid: A method for evaluating potential Covid-19 transmission events".
##' Illingworth C., Hamilton W., Jackson C. et al. Under preparation.
##'
##' @examples
##'
##' ## Example data supplied with the package
##' pat_file <- system.file("extdata", "Example_genetic_temporal_data.csv", package="a2bcovid")
##' hcw_loc_file <- system.file("extdata", "Example_movement_file.csv", package="a2bcovid")
##' ali_file <- system.file("extdata", "Example_sequences.fa", package="a2bcovid")
##' pat_loc_file <- system.file("extdata", "Example_pat_loc_file.csv", package="a2bcovid")
##'
##' res <- a2bcovid(pat_file = pat_file, hcw_loc_file = hcw_loc_file, ali_file = ali_file,
##'                 pat_loc_file = pat_loc_file)
##' plot_a2bcovid(res, hi_from="from_hcw", hi_to="to_hcw")
##'
##' @seealso \code{\link{plot_a2bcovid}}
##'
##' @importFrom Rcpp evalCpp
##' @importFrom shiny runApp
##' @importFrom rstudioapi viewer
##' @importFrom utils read.csv write.csv
##' @importFrom stats approx pnorm qnorm
##' @importFrom magrittr %>%
##' @useDynLib a2bcovid, .registration = TRUE
##'
##' @export
a2bcovid <- function(
  pat_file,
  hcw_loc_file = "",
  ali_file = "",
  pat_loc_file = "",
  pa=97.18750, pb=0.268908, po=25.625,
  smu=1.434065, ssigma=0.6612,
  ucta=2.5932152095707406, uctb=3.7760060663975437, ucto=3.112080041460921,
  uct_mean=6.67992,
  evo_rate  = 0.0008,
  seq_noise = 0.41369,
  chat = 0.5,
  max_n = 10,
  min_qual = 0.8,
  calc_thresholds=FALSE,
  diagnostic =FALSE,
  hcw_default = 0.5714286,
  pat_default = 1,
  use_all_seqs = 0,
  symptom_uncertainty_calc = 0
)
{
  ali_file <- check_file(ali_file)
  pat_file <- check_file(pat_file)
  hcw_loc_file <- check_file(hcw_loc_file)
  pat_loc_file <- check_file(pat_loc_file)
  if (pat_loc_file != "") {
    pat_loc_format <- guess_loc_format(pat_loc_file)
    if (pat_loc_format=="long")
      pat_loc_file <- long_to_wide(pat_loc_file)
  }
  params <- list(pa=pa, pb=pb, po=po, smu=smu, ssigma=ssigma,
                 ucta=ucta, uctb=uctb, ucto=ucto, uct_mean=uct_mean,
                 rate=evo_rate, seq_noise=seq_noise,
				 chat=chat,
                 max_n=max_n, min_qual=min_qual,
                 ali_file=ali_file, pat_file=pat_file,
                 mov_file=hcw_loc_file, ward_file=pat_loc_file,
                 calc_thresholds=calc_thresholds,
                 diagnostic=diagnostic,
				 hcw_location_default=hcw_default,
				 pat_location_default=pat_default,
				 use_all_seqs=use_all_seqs,
				 symptom_uncertainty_calc=symptom_uncertainty_calc)
  res <- .Call(`_a2bcovid_mainC`, params)
  ## default factor levels as in data order, not sorted alphabetically
  res$from <- factor(res$from, unique(res$from))
  res$to <- factor(res$to, unique(res$to))
  ## converts from 0 based to 1 based index
  res$ordered_i <- res$ordered_i + 1
  res$ordered_j <- res$ordered_j + 1
  res$likelihood[res$likelihood==-1e+10] <- -Inf
  thresh <- if (ali_file=="") thresholds_noseq else thresholds_seq
  suppressWarnings(res$p <- approx(x = thresh$lik, y = thresh$p, xout=res$likelihood, yleft=1, yright=0)$y)
  res$p <- 1 - res$p
  res$consistency <- ordered(res$consistency,
                             levels=c("Unlikely","Borderline","Consistent"))
  res <- res[res$from != res$to,]
  attr(res, "files") <- list(pat=pat_file, ali=ali_file,  pat_loc=pat_loc_file, hcw_loc=hcw_loc_file)
  class(res) <- c("a2bcovid", class(res))
  res
}

##' @export
plot.a2bcovid <- function(x,...){
  plot_a2bcovid(x, ...)
}

guess_loc_format <- function(filename){
  dat <- read.csv(filename)
  if ("start_date" %in% names(dat)) return("long")
  else if ("StartDate_0" %in% names(dat)) return("wide")
  else stop("Could not determine patient location file format. No column named either `start_date` or `StartDate_0`")
}

## Check for file existence

check_file <- function(filename){
  if (!(filename=="")){
    if (!file.exists(filename))
      stop(sprintf("File `%s` not found", filename))
    filename <-  normalizePath(filename)
  }
  filename
}


##' Plot results of an a2bcovid analysis
##'
##' Plots a grid of colours indicating likelihood of transmission paths between
##' each pair of individuals.
##'
##' @param x Data frame returned by \code{\link{a2bcovid}}.
##'
##' @param cluster If \code{TRUE} (the default) then the individual IDs are rearranged  using
##'   a heuristic clustering method so that apparent clusters of infections
##'   appear contiguously in the plot.   If \code{FALSE} then
##'   individuals are arranged in their original order in the data.
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
##' @param hi_lab Legend to describe which individuals are highlighted.  By default this is "Healthcare workers".
##'
##' @param palette Colour palette, passed to
##'   \code{\link[ggplot2]{scale_fill_brewer}}.  If omitted, a default will be chosen
##'
##' @param continuous  If \code{TRUE} then the p-values are plotted on a
##'   continuous colour scale.  Currently only implemented with the default
##'   colour palette.  If \code{FALSE} (the default), then the p-values are
##'   classified into ranges and plotted as a categorical variable.
##'
##' @param direction Direction of colours in the brewer palettes. Defaults to 1.  Change to -1 to
##'   reverse the order of colours.
##'
##' @return A \pkg{ggplot2} plot object.
##'
##' @seealso \code{\link{a2bcovid}}
##'
##' @import ggplot2
##' @importFrom scales squish
##'
##' @export
plot_a2bcovid <- function(x, cluster = TRUE,
                          hi_from="from_hcw", hi_to="to_hcw", hi_col="red",
                          hi_lab=NULL,
                          palette=NULL,
                          continuous=FALSE,
                          direction = 1){
  midp <- pnorm((qnorm(0.05) + qnorm(0.01))/2)
  if (continuous) {
    scale_chosen <- scale_fill_gradientn(
      colours = c( "#D64E47", "#FCF9DA", "#3C87C8"),
      values = scales::rescale(c(0.1, midp, 0)),
      breaks=c(0.1, 0.05, 0.01, 0),
      limits=c(0, 0.1),
      oob=scales::squish,
      guide="legend"
    )
    leg_lab <- "p-value"
  } else {
    if (is.null(palette)) {
      cols_default <- c("Consistent"="#D64E47",
                        "Borderline"="#FCF9DA",
                        "Unlikely"="#3C87C8")
      scale_chosen <- scale_fill_manual(values = cols_default, breaks=names(cols_default))
    } else {
      scale_chosen <- scale_fill_brewer(palette = palette,  direction=direction)
    }
    leg_lab <- ""
  }
  ## TODO match with rownames(RColorBrewer::brewer.pal.info)

  if (cluster) {
    x$from <- factor(x$from, levels = levels(x$from)[x$ordered_i[!duplicated(x$from)]])
    x$to <- factor(x$to, levels = levels(x$to)[x$ordered_j[!duplicated(x$to)]])
  }

  xhi_from <- xhi_to <- NULL
  if (!is.null(hi_from)) {
    if (!(hi_from %in% names(x))) stop(sprintf("`%s` not found in `x`", hi_from))
    x_from <- x[!duplicated(x$from),]
    xhi_from <- x_from[match(levels(factor(x$from)), x_from$from), hi_from]
    from_cols <- ifelse(xhi_from, hi_col, "black")
  } else from_cols <- "black"
  if (!is.null(hi_to)) {
    if (!(hi_to %in% names(x))) stop(sprintf("`%s` not found in `x`", hi_from))
    x_to <- x[!duplicated(x$to),]
    xhi_to <- x_to[match(levels(factor(x$to)), x_to$to), hi_to]
    to_cols <- ifelse(xhi_to, hi_col, "black")
  } else to_cols <- "black"
  if (is.null(hi_lab) & ((length(xhi_from)>0)||(length(xhi_to)>0)))
    hi_lab <- sprintf("Healthcare workers labelled in %s",hi_col)

  fill_var <- if (continuous) "p" else "consistency"
  ggplot(x, aes_string(y="from", x="to")) +
    geom_tile(aes_string(fill=fill_var), colour = "black") +
    theme(axis.text.x = ggtext::element_markdown(angle = 90, vjust=0.5, colour = to_cols),
          axis.text.y = ggtext::element_markdown(colour = from_cols),
          plot.subtitle = ggtext::element_markdown(colour=hi_col)
          ) +
    labs(subtitle=hi_lab,
         fill=leg_lab) +
    scale_chosen +
    theme(panel.background = element_rect(fill = 'white', colour = 'white'))
}
