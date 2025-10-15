.onLoad <- function(libname, pkgname) {
  # Be polite: only *use* an existing env; never create it automatically.
  try({
    if (!nzchar(Sys.getenv("RETICULATE_PYTHON", ""))) {
      envname <- getOption("dncitPaper.python_env", "dncit-paper")
      envs <- tryCatch(reticulate::conda_list()$name, error = function(e) character())
      if (envname %in% envs) reticulate::use_condaenv(envname, required = FALSE)
    }
  }, silent = TRUE)
}
