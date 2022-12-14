% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Clustering.R
\name{z2_score}
\alias{z2_score}
\title{Calculate Chi-squared statistics based on two Normal Distributions}
\usage{
z2_score(mu1, se1, mu2, se2)
}
\arguments{
\item{mu1}{mean of the first normal distribution}

\item{se1}{standard error of the first normal distribution}

\item{mu2}{mean of the second normal distribution}

\item{se2}{standard error of the second normal distribution}
}
\description{
`z2_score` computes and returns Chi-squared statistics based on two normal distributions
}
\details{
Given \eqn{X \sim N(\mu_x,\sigma_x^2)} and \eqn{Y \sim N(\mu_y,\sigma_y^2)},
I calculate \deqn{(\frac{(\mu_x - \mu_y)}{(\sigma_x^2 +\sigma_y^2)})^2}
}
