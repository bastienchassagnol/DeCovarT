# define the number of cores to be used in an automated testing environment

define_seed_environment <- function(env = parent.frame()) {
  # your fiddly code to create a useful_thing goes here
  chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
  if (nzchar(chk) && chk == "TRUE") {
    num_cores <- 2L
  } else if (.Platform$OS.type == "unix") {
    num_cores <- parallel::detectCores() - 1
  } else {
    num_cores <- 1
  }
  # RNGkind("L'Ecuyer-CMRG"); set.seed(3)
  return(num_cores)
}
