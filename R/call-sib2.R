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
#' input <- system.file("extdata", "data2", package = "sib2r")
#' input
#' sib_run <- c(
#'    date1 = 1960111601,
#'    date2 = 1961123024,
#'    nlai = 65,
#'    aeromet = 1,
#'    isnow = 0,
#'    ipbl = 1,
#'    ilw = 3,
#'    itrunk = 20,
#'    ivtype = 6,
#'    istype = 3,
#'    idirr = 0
#')
#'sib_ini <- c(
#'    zlat = -20.5,
#'    zwind = 45.0,
#'    zmet = 45.0,
#'    dtt = 3600.0,
#'    tc = 298.0,
#'    tg = 298.0,
#'    td = 297.0,
#'    gwdep = 0.0,
#'    gmudmu = 1.0,
#'    app = 0.0001,
#'    bpp = 20.0,
#'    cpp = 0.9999,
#'    capac = c(0.0, 0.0),
#'    snoww = c(0.0, 0.0),
#'    www = c(0.75, 0.85, 0.98)
#')
#'sib_soil <- c(
#'    bee = 7.797,
#'    decay = 0.85,
#'    poros = 0.458,
#'    phsat = -0.2,
#'    satco = 3.5E-06,
#'    slpp = 0.0
#')
#' sib_morpho <- c(
#'     zs = 0.05,
#'     dsfc = 0.02,
#'     g1 = 1.449,
#'     ztz0 = 11.785,
#'     sodep = 2.0,
#'     slope = 0.08
#')
#'sib_physio <- c(
#'    shti = 0.3,
#'    slti = 0.2,
#'    trda = 1.3,
#'    trdm = 328.0,
#'    trop = 298.0,
#'    btheta = 0.95
#')
#'sib_parms <- c(
#'    z1 = 0.3,
#'    z2 = 1.2,
#'    zc = 0.6,
#'    chil = -0.3,
#'    leafw = 0.01,
#'    leafl = 0.3,
#'    vcover = 0.95,
#'    rootd = 1.0,
#'    phc = -200.0,
#'    tranlv = 0.07,
#'    tranln = 0.25,
#'    trandv = 0.22,
#'    trandn = 0.38,
#'    reflv = 0.11,
#'    refln = 0.58,
#'    refdv = 0.36,
#'    refdn = 0.58,
#'    sorefv = 0.12,
#'    sorefn = 0.20,
#'    effcon = 0.05,
#'    gradm = 4.0,
#'    binter = 0.04,
#'    respcp = 0.02,
#'    atheta = 0.95,
#'    hlti = 290.16,
#'    hhti = 313.16,
#'    vmax0 = 3.0e-5
#')
#' out <- call_sib2(
#'    file = input,
#'    out.dir = normalizePath("~/Desktop/"),
#'    id.sim = "_run1",
#'    run.pars = sib_run,
#'    ini.pars = sib_ini,
#'    soil.pars = sib_soil,
#'    morpho.pars = sib_morpho,
#'    physio.pars = sib_physio,
#'    veg.pars = sib_parms
#' )
#' out
#' 
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
    sib2_output
    # out_file <- file.path(
    #   normalizePath(out.dir),
    #   paste0("sib2diag", id.sim, ".txt")
    # )
    # if(file.exists(out_file)) {
    #   if(verbose){
    #     size_file <- round(file.info(out_file)$size/1e6, 2)
    #     message("Output file save:", paste0(out_file, "\n", size_file))
    #     rm(size_file)
    #   }
    #   return(out_file)
    # }
    # warning( paste0("File: \n", out_file, "\not found.") )
    # return(NA_character_)
  }
