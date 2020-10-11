check_date <- function(mydate) {
  valdate <- tryCatch(!is.na(as.Date(mydate, "%d/%m/%Y")),
                      error = function(err) {FALSE})
  bad_dates <- which(!valdate)
  if (any(bad_dates)) {
    ret <- sprintf("%s not in dd/mm/YYYY date format", bad_row_str(bad_dates))
  } else ret <- TRUE
  ret
}

bad_row_str <- function(badid){
  if (length(badid) == 1) errstr <- sprintf("row %s", badid)
  else {
    if (length(badid) == 2) errstr <- sprintf("rows %s,%s",badid[1],badid[2])
    if (length(badid) > 2) errstr <- sprintf("rows %s,%s and others",badid[1],badid[2])
  }
  errstr
}

validate_patfile <- function(pat_file){
  csvtry <- try(f <- read.csv(pat_file))
  if (inherits(csvtry, "try-error")){
    message(sprintf("Reading the patient CSV file returned the following error:\n`%s`",
                    as.character(csvtry)))
  }
  val2 <- check_date(f[,2])
  if (!isTRUE(valdate))
    message(sprintf("Invalid dates found in column 2: %s", valdate))
  bad_entries <- which(!(f[,3] %in% 1:3))
  if (length(bad_entries) > 1) {
    message("Invalid entries found in column 3: %s should have the value 1, 2 or 3", bad_row_str(bad_entries))
  }
  bad_entries <- which(!(f[,4] %in% 1:3))
  if (length(bad_entries) > 1) {
    message("Invalid entries found in column 4: %s should have the value 1, 2 or 3", bad_row_str(bad_entries))
  }
  val6 <- check_date(f[,6])
  if (!isTRUE(valdate))
    message(sprintf("Invalid dates found in column 6: %s", valdate))

  ## col 1 indiv ID. OK can be anything

  ## col 5 sequence ID code ?? any restriction?   should cross validate with seq file

  ## can it have extra cols
}

#pat_file <- system.file("extdata", "Example_genetic_temporal_data.csv", package="a2bcovid")
#read.csv()
#csvtry <- try(read.csv("foo"))

#validate_patfile("foo")
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
