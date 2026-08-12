## ============================================================================
## data_ingestion.R
##
## Reads both TPS files using geomorph::readland.tps() and links them to the
## covariate spreadsheets by the (view/group/id) code embedded in each image
## filename
##
## REQUIRES: geomorph, readxl, dplyr
## install.packages(c("geomorph", "readxl", "dplyr"))
## ============================================================================

## ---- 1. Read TPS files -------------------------------------------------
## specID = "imageID" tells geomorph to name each specimen (3rd-dimension
## slice) using the IMAGE= line from the TPS file, which is what we need to
## parse out the view/group/id code. geomorph strips the directory path
## automatically
sink(tempfile())

ant_land <- readland.tps("ANTERIOR__IMAGES_SAMPLE_AND_CONTROL.TPS",
                         specID = "imageID", warnmsg = TRUE)
lat_land <- readland.tps("LATERAL_IMAGES_4_SAMPLE_AND_CONTROL.TPS",
                         specID = "imageID", warnmsg = TRUE)
sink()

## ---- 2. Parse the (view, group, id) code from specimen names -----------
parse_specimen_code <- function(names_vec) {
  m <- regmatches(names_vec,
                  regexpr("[AL][SC][0-9]{3}", names_vec, ignore.case = TRUE, perl = TRUE))
  data.frame(
    specimen = names_vec,
    view  = toupper(substr(m, 1, 1)),
    group = toupper(substr(m, 2, 2)),
    id    = as.integer(substr(m, 3, 5)),
    stringsAsFactors = FALSE
  )
}

ant_codes <- parse_specimen_code(dimnames(ant_land)[[3]])
lat_codes <- parse_specimen_code(dimnames(lat_land)[[3]])

stopifnot(sum(is.na(ant_codes$view)) == 0)
stopifnot(sum(is.na(lat_codes$view)) == 0)

check_ids <- function(codes, label) {
  for (g in c("S", "C")) {
    ids <- sort(codes$id[codes$group == g])
    missing <- setdiff(1:80, ids)
    dup <- ids[duplicated(ids)]
    if (length(missing) > 0) warning(label, " group ", g, " missing IDs: ", paste(missing, collapse = ","))
    if (length(dup) > 0) warning(label, " group ", g, " duplicate IDs: ", paste(dup, collapse = ","))
  }
}
check_ids(ant_codes, "Anterior")
check_ids(lat_codes, "Lateral")

## ---- 3. Read and link covariate sheets ----------------------------------
read_covariates <- function(path, group_code) {
  df <- read_excel(path, sheet = "Sheet1")
  names(df) <- trimws(names(df))
  df <- df[!is.na(df[["S/N(ID)"]]), ]
  df[["id"]] <- as.integer(df[["S/N(ID)"]])
  df[["group"]] <- group_code
  chr_cols <- sapply(df, is.character)
  df[chr_cols] <- lapply(df[chr_cols], trimws)
  df
}

ctrl <- read_covariates("Data_Sheet_of_Control.xlsx", "C")
samp <- read_covariates("Data_sheet_of_the_sample.xlsx", "S")
common_cols <- intersect(names(ctrl), names(samp))
covar <- bind_rows(ctrl[common_cols], samp[common_cols])
covar$CPD <- as.integer(covar$group == "S")
