% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_gen.R
\name{data_gen}
\alias{data_gen}
\title{Title}
\usage{
data_gen(
  seed,
  idx_sample = NULL,
  n_sample = NULL,
  idx_beta2 = NULL,
  beta2s = NULL,
  n = NULL,
  post_non_lin = 4,
  eps_sigmaX = 0,
  eps_sigmaY = 1,
  eps_sigmaZ = 0,
  embedding_orig = "fastsurfer",
  embedding_obs = "fastsurfer",
  confounder = "AS",
  response = "simulated"
)
}
\arguments{
\item{seed}{random seed for current data generation}

\item{idx_sample}{index of current sample size}

\item{n_sample}{list of sample sizes}

\item{idx_beta2}{index of current beta2}

\item{beta2s}{list of beta2s}

\item{n}{sample size}

\item{post_non_lin}{post non-linearity function}

\item{eps_sigmaX}{noise level for X}

\item{eps_sigmaY}{noise level for Y}

\item{eps_sigmaZ}{noise level for Z}

\item{embedding_orig}{embedding map used to obtain original feature representations (used to simulate conditional association or CI)}

\item{embedding_obs}{embedding map used to obtain observed feature representations (used for CI testing)}

\item{confounder}{confounder used to simulate Y}

\item{response}{response variable}
}
\value{
list of X_obs, Y, Z
}
\description{
Title
}
