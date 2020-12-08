##' Convert patient location data for an a2bcovid analysis from long to wide
##' format
##'
##' @param long_file A path name to a CSV file in long format.
##'
##' @return A path name to a temporary file containing the equivalent data in
##'   wide format. This can be read with \code{\link{read.csv}}.
##'
##'  The names of the columns in the wide and long formats are both documented
##'  in the \code{\link{a2bcovid}} help page, argument \code{pat_loc_file}.
##'
##' @export
long_to_wide <- function(long_file){

    dat <- read.csv(long_file, col.names = c("patient", "from_ward", "start_date",
                                               "to_ward", "end_date"),
                    stringsAsFactors = FALSE)
    #dat$start_date <- as.Date(strptime(dat$start_date, format="%d/%m/%Y"))
    #dat$end_date <- as.Date(strptime(dat$end_date, format="%d/%m/%Y"))
    #dat$patient <- factor(dat$patient, levels=unique(dat$patient))

    # dat <- dat[order(dat$patient, dat$start_date),]
    # Find first and last observation for each ward stay
    # then pick start date from first and end date from last
    # Should be left with one row per ward stay
    ptward <- paste(dat$patient, dat$from_ward, sep="_")
    n <- length(ptward)
    ptward_prev <- c(NA, ptward[1:(n-1)])
    ptward_next <- c(ptward[2:n], NA)
    wardobsid_first <- (ptward != ptward_prev) | is.na(ptward_prev)
    wardobsid_last <- (ptward != ptward_next) | is.na(ptward_next)
    enddates <- dat$end_date[wardobsid_last]
    dat <- dat[wardobsid_first,]
    dat$end_date <- enddates
    #dat$start_date <- format(dat$start_date, format="%d/%m/%Y")
    #dat$end_date <- format(dat$end_date, format="%d/%m/%Y")
    dat$stayid <- sequence(table(dat$patient)[unique(dat$patient)])

    locs <- tidyr::pivot_wider(dat[,c("patient","from_ward","stayid")],
                               names_from="stayid", names_prefix="LocationName",
                               values_from="from_ward")
    starts <- tidyr::pivot_wider(dat[,c("patient","start_date","stayid")],
                                 names_from="stayid", names_prefix="StartDate",
                                 values_from="start_date")
    ends <- tidyr::pivot_wider(dat[,c("patient","end_date","stayid")],
                               names_from="stayid", names_prefix="EndDate",
                               values_from="end_date")

    # One row per patient, with a triple of columns for each ward stay
    widedata <- data.frame(
        patient_study_id = locs$patient,
        ward_cluster_network = "A",
        hcw_status = "patient",
        patient_movement_data_available = "patient_moves_available"
    )
    firstcols <- colnames(widedata)
    widedata <- cbind(widedata, locs[,-1], starts[,-1], ends[,-1])
    nstays <- max(dat$stayid)
    cols <- paste0(rep(c("LocationName","StartDate","EndDate"), nstays),
                  rep(1:nstays, each=3))
    ret <- widedata[,c(firstcols, cols)]
    wide_file <- tempfile()
    write.csv(ret, wide_file, row.names = FALSE)
    wide_file
}

##' Convert patient location data for an a2bcovid analysis from wide to long
##' format
##'
##' @param wide_file A path name to a CSV file in wide format.
##'
##' @return A path name to a temporary file containing the equivalent data in
##'   long format. This can be read with \code{\link{read.csv}}.
##'
##'  The names of the columns in the wide and long formats are both documented
##'  in the \code{\link{a2bcovid}} help page, argument \code{pat_loc_file}.
##'
##' @export
wide_to_long <- function(wide_file){
    dat <- read.csv(wide_file, stringsAsFactors = FALSE, header=TRUE)
    firstcols <- c("patient_study_id","ward_cluster_network","hcw_status",
                   "patient_movement_data_available")
    triples <- setdiff(names(dat), firstcols)
    stayno <- as.numeric(gsub("[[:alpha:]]+_([0-9]+)","\\1", triples))
    stayname <- gsub("([[:alpha:]]+)_[0-9]+","\\1", triples)
    badnames <- !stayname %in% c("LocationName","StartDate","EndDate")
    if (any(badnames)) {
        firstbadname <- stayname[which(badnames)[1]]
        stop("Column name %s doesn't match any of the required names", firstbadname)
    }
    datlong <- tidyr::pivot_longer(dat, cols=tidyselect::all_of(triples),
                                   names_to=c(".value", "visit"),
                                   names_sep = "_") %>%
        dplyr::select("patient_study_id", "visit", "LocationName", "StartDate", "EndDate") %>%
        dplyr::rename("from_ward" = "LocationName",
                      "start_date" = "StartDate",
                      "end_date" = "EndDate") %>%
        dplyr::filter(.data$from_ward != "") %>%
        dplyr::mutate("to_ward" = dplyr::lead(.data$from_ward))
    datlong <- datlong %>%
        dplyr::mutate("to_ward" = ifelse((!is.na(dplyr::lead(.data$patient_study_id)) &
                                     .data$patient_study_id == dplyr::lead(.data$patient_study_id)),
                                .data$to_ward, "Discharge")) %>%
        dplyr::select("patient_study_id", "from_ward", "start_date", "to_ward", "end_date")
    long_file <- tempfile()
    write.csv(datlong, long_file, row.names = FALSE)
    long_file
}


## Not currently used
check_long_dates <- function(dat){
    baddates <- which(dat$start_date > dat$end_date)
    if (length(baddates) > 0) {
        if (length(baddates) > 1) {
            extramsg <- " and others"
        }
        stop(sprintf("Start date %s after end date %s for patient %s stay on ward %s%s",
                     dat$start_date[baddates[1]], dat$end_date[baddates[1]],
                     dat$patient[baddates[1]], dat$from_ward[baddates[1]]), extramsg)
    }
}
