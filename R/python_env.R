#' Create the Python environment once (idempotent)
#' @param envname conda env name
#' @param cuda one of "cpu","cu121","cu128" (A100: use "cu128" for Torch 2.8, or "cu121" for Torch 2.3)
#' @param force recreate env if it exists
#' @param yml optional path to environment.yml (auto-detected if NULL)
#' @param torch_ver torch version spec (pip); e.g. "2.3.*" or "2.8.*"
#' @param monai_ver MONAI version; e.g. "1.5.1"
dncit_setup_env <- function(
  envname   = "dncit-paper",
  cuda      = c("cu121","cu128","cpu"),
  force     = FALSE,
  yml       = NULL,
  torch_ver = "2.3.*",   # change to "2.8.*" if you want the newer stack by default
  monai_ver = "1.5.1"
) {
  cuda <- match.arg(cuda)

  # locate environment.yml (installed pkg or dev tree)
  if (is.null(yml) || !file.exists(yml)) {
    yml <- system.file("python", "environment.yml", package = "dncitPaper")
    if (yml == "" || !file.exists(yml)) {
      cand <- c("inst/python/environment.yml","python/environment.yml","environment.yml")
      yml <- cand[file.exists(cand)][1]
    }
  }
  if (is.null(yml) || !file.exists(yml)) stop("environment.yml not found.", call. = FALSE)

  # ensure conda/miniconda
  if (is.null(tryCatch(reticulate::conda_binary(), error = function(e) NULL))) {
    message("Installing Miniconda (first-time only)…")
    reticulate::install_miniconda()
  }

  # (re)create env
  envs <- tryCatch(reticulate::conda_list()$name, error = function(e) character())
  if (force && (envname %in% envs)) {
    message("Removing existing env '", envname, "' (force=TRUE)…")
    reticulate::conda_remove(envname = envname)
    envs <- setdiff(envs, envname)
  }
  if (!(envname %in% envs)) {
    message("Creating env '", envname, "' from ", yml, " …")
    reticulate::conda_create(envname = envname, yaml = yml)
  } else {
    message("Env '", envname, "' already exists; skipping creation.")
  }

  # core scientific stack (safe to rerun)
  reticulate::conda_install(
    envname  = envname,
    packages = c(
      "numpy","pandas","scipy","scikit-learn",
      "nibabel","tqdm","pillow","matplotlib",
      "simpleitk","xlrd","openpyxl"
    ),
    channel  = "conda-forge"
  )

  # pick the PyTorch wheel index for the requested CUDA runtime
  idx_url <- switch(cuda,
    "cpu"   = "https://download.pytorch.org/whl/cpu",
    "cu121" = "https://download.pytorch.org/whl/cu121",
    "cu128" = "https://download.pytorch.org/whl/cu128"
  )

  # Torch (pip in conda env) – force/upgrade so we can switch between 2.3 and 2.8 cleanly
  reticulate::conda_install(
    envname     = envname,
    packages    = sprintf("torch==%s", torch_ver),
    pip         = TRUE,
    pip_options = c("--index-url", idx_url, "--upgrade", "--force-reinstall")
  )

  # MONAI (compatible with Torch up to 2.8 per 1.5.1 notes)
  reticulate::conda_install(
    envname  = envname,
    packages = sprintf("monai==%s", monai_ver),
    pip      = TRUE
  )

  invisible(TRUE)
}

# Activate env (respects RETICULATE_PYTHON if set)
dncit_use_env <- function(envname = "dncit-paper", required = FALSE) {
  rp <- Sys.getenv("RETICULATE_PYTHON", unset = "")
  if (nzchar(rp)) { reticulate::use_python(rp, required = required); return(invisible(TRUE)) }
  reticulate::use_condaenv(envname, required = required); invisible(TRUE)
}

# Quick sanity
dncit_python_ok <- function() {
  try({
    print(reticulate::py_config())
    if (!reticulate::py_module_available("torch"))  warning("torch not found")
    if (!reticulate::py_module_available("monai"))  warning("monai not found")
    if (!reticulate::py_module_available("nibabel")) warning("nibabel not found")
  }); invisible(TRUE)
}
