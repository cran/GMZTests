% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DMC_test.R
\name{dmc.test}
\alias{dmc.test}
\title{Statistical test for Multiple Detrended Cross-Correlation Coefficient}
\usage{
dmc.test(N, k, m, method, nu, rep)
}
\arguments{
\item{N}{An integer value for the time series length.}

\item{k}{An integer value indicating the boundary of the division \eqn{(N/k)}.
The smallest value of \eqn{k} is \eqn{4}.}

\item{m}{an integer value or a vector of integer values indicating the size of the window for the polinomial fit.}

\item{method}{A character string indicating which correlation coefficient is to be used. If method = "rhodcca" the dmc coefficient is generated from the DCCA coefficient. If method = "dmca", the dmc coefficient is generated from the DMCA coefficient.}

\item{nu}{An integer value. See the DCCA package.}

\item{rep}{An integer value indicating the number of repetitions.}
}
\value{
An list containing "timescale", parameters of beta distribution: "shape1", "se1","shape2","se2" and confidence interval: "CI_0.90_uppper", "CI_0.95_uppper", "CI_0.99_uppper".
}
\description{
This function performs the statistical test for DMC Cross-Correlation Coefficient based in White Gaussian Noise process.
}
\details{
This function include following measures: w, timescale, dmc, rhodcca_yx1, rhodcca_yx2, rhodcca_x1x2
}
\examples{
dmc.test(N=100, k=10, m=c(4:6), method="rhodcca", nu=0, rep=10)
dmc.test(N=100, k=10, m=c(4:6), method="dmca", nu=0, rep=10)

}
\references{
SILVA-FILHO,A.M; ZEBENDE,G.; CASTRO,A.P.; GUEDES,E. Statistical test for multiple detrended cross-correlation coefficient, Physica A, v.562, 125285, 2021.

KRISTOUFEK, L. Detrending moving-average cross-correlation coefficient: Measuring cross-correlations between non-stationary series. PHYSICA A, v.406, p.169-175, 2014.
}
