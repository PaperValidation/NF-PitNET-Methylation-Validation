# Centralized dataset paths loaded from the project root .env.
#
# Source of truth for DATASET_PATH (raw/permanent, G:/CMCC Sync) and
# DATASET_PATH_TEMP (intermediate/temp, E:/Data Temp). Scripts should source
# this file and read paths from the returned list instead of hardcoding.
#
# Usage (from an R script inside scripts/):
#   source("C:/Users/tidus/PycharmProjects/Neurosurgery-NF-PitNET/scripts/paths.R")
#   BASE  <- paths$DATASET_PATH_WSL          # /mnt/g/...   (WSL form)
#   TEMP  <- paths$DATASET_PATH_TEMP_WSL     # /mnt/e/Data Temp/NF-PitNET

paths <- local({
  proj <- "C:/Users/tidus/PycharmProjects/Neurosurgery-NF-PitNET"

  env <- new.env()
  con <- tryCatch(file(file.path(proj, ".env"), "r"), error = function(e) NULL)
  if (!is.null(con)) {
    on.exit(close(con), add = TRUE)
    while (TRUE) {
      line <- readLines(con, n = 1)
      if (length(line) == 0) break
      line <- trimws(line)
      if (nchar(line) == 0 || substr(line, 1, 1) == "#") next
      m <- regexpr("^([A-Za-z0-9_]+)\\s*=\\s*'?([^'\n]*)'?", line, perl = TRUE)
      if (m > 0) {
        key <- sub("^([A-Za-z0-9_]+).*", "\\1", line)
        val <- sub("^([A-Za-z0-9_]+)\\s*=\\s*'?(.*?)'?\\s*$", "\\2", line, perl = TRUE)
        assign(key, val, envir = env)
      }
    }
  }

  win_to_wsl <- function(w) {
    if (is.null(w) || nchar(w) < 2) return(w)
    drive <- tolower(substr(w, 1, 1))
    rest <- sub("^.:/?", "", gsub("\\\\", "/", w))
    paste0("/mnt/", drive, "/", rest)
  }

  ds <- ifelse(is.null(env$DATASET_PATH),
               "G:/CMCC Sync/Medical Data/Neurosurgery/NF-PitNET", env$DATASET_PATH)
  dt <- ifelse(is.null(env$DATASET_PATH_TEMP),
               "E:/Data Temp/NF-PitNET", env$DATASET_PATH_TEMP)

  list(
    PROJECT_PATH          = proj,
    PROJECT_PATH_WSL      = win_to_wsl(proj),
    DATASET_PATH          = ds,
    DATASET_PATH_WSL      = win_to_wsl(ds),
    DATASET_PATH_TEMP     = dt,
    DATASET_PATH_TEMP_WSL = win_to_wsl(dt)
  )
})
