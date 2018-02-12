#' Run SiB2 
#'
#' @param file A character string with the Input file with meteorological forcing dataset 
#' @param out.dir  A character string with the path where the output file ('sib2diag<id>.txt') is saved. 
#' @param id.sim A character string to identify simulation (e.g. '_run1').
#' @param run.pars A numeric vector with model setup parameters
#'   * `date1`: 1960111601
#'   * `date2` : 1961123024
#'   * `lai.max` : 65
#'   * `aeromet` : 1
#'   * `isnow` : 0
#'   * `ipbl` : 1
#'   * `ilw` : 3
#'   * `itrunk` : 20
#'   * `ivtype` : 6
#'   * `istype` : 3
#'   * `idirr` : 0 
#' 
#' @param ini.pars initial conditions
#' @param soil.pars soil parameters
#' @param morpho.pars morphological parameters
#' @param physio.pars physiological parameters
#' @param veg.pars vegetation parameters (miscelaneous)
#' @param verbose logical, default is TRUE, to print messages
#'
#' @useDynLib sib2r, .registration = TRUE
#' @return
#' A string character with the output file name ('sib2diag<id>.txt')
#' @export
#'
#' @examples
call_sib2 <-
  function(file,
           out.dir,
           id.sim,
           run.pars,
           ini.pars,
           soil.pars,
           morpho.pars,
           physio.pars,
           veg.pars,
           verbose = TRUE
  ) {
    n1 <- nchar(file)
    n2 <- nchar(out.dir)
    n3 <- nchar(id.sim)
    
    # Call the SiB2 model
    sib2_output <- .Fortran(
      "sib2_offline_lib",
      as.integer(n1),
      as.character(file),
      as.integer(n2),
      as.character(out.dir),
      as.integer(n3),
      as.character(id.sim),
      as.integer(run.pars),
      as.double(ini.pars),
      as.double(soil.pars),
      as.double(morpho.pars),
      as.double(physio.pars),
      as.double(veg.pars)
    )

    out_file <- file.path(
      normalizePath(out.dir),
      paste0("sib2diag", id.sim, ".txt")
    )
    if(file.exists(out_file)) {
      if(verbose){
        size_file <- round(file.info(out_file)$size/1e6, 2)
        message("Output file save:", paste0(out_file, "\n", size_file))
        rm(size_file)
      }
      return(out_file)
    }
    warning( paste0("File: \n", out_file, "\not found.") )
    return(NA_character_)
  }
