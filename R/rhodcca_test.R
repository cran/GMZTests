#' @title Statistical test for detrended cross-correlation coefficient
#'
#' @description This function performs the statistical test for RHODCCA cross-correlation coefficient based in White Gaussian Noise process.
#'
#' @details This function include following measures: timescale and cross-correlation yx.
#'
#' @param N An integer value for the time series length.
#'
#' @param k An integer value indicating the boundary of the division \eqn{(N/k)}.
#'          The smallest value of \eqn{k} is \eqn{4}.
#'
#' @param m an integer value or a vector of integer values indicating the size of the window for the polinomial fit.
#'
#' @param nu An integer value. See the DCCA package.
#'
#' @param rep An integer value indicating the number of repetitions.
#'
#' @return An list containing "timescale","mean", "sd" and confidence interval: "CI_0.90", "CI_0.95", "CI_0.99".
#'
#' @examples
#' rhodcca.test(N=100, k=10, m=c(4:6), nu=0, rep=10)
#'
#' @references
#' B. Podobnik, Z.-Q. Jiang, W.-X. Zhou, H. E. Stanley, Statistical tests for power-law cross-correlated processes, Phys. Rev.,E 84, 066118, 2011.
#'
#' @importFrom stats rnorm qbeta
#' @importFrom DCCA rhodcca
#' @importFrom stats rnorm qnorm sd
#' @importFrom fitdistrplus fitdist
#'
#' @export
 rhodcca.test <- function(N,k,m,nu,rep){

    yx <- matrix(data = NA, nrow = rep, ncol = length(m), byrow = TRUE)

    for(i in 1:rep){
     y <- stats::rnorm(N, mean=0, sd=1)
     x <- stats::rnorm(N, mean=0, sd=1)

        yx[i,]  <- DCCA::rhodcca(y, x, m=m, nu=nu)$rhodcca
     }
        mean <- apply(yx, 2, mean)
          sd <- apply(yx, 2, sd)

        CI1 <- stats::qnorm(c(0.90+((1-0.90)/2)), 0,1)*sd
        CI2 <- stats::qnorm(c(0.95+((1-0.95)/2)), 0,1)*sd
        CI3 <- stats::qnorm(c(0.99+((1-0.99)/2)), 0,1)*sd

     return(list(timescale=m,
                  mean = mean,
                  sd = sd,
                  CI_0.90 = CI1,
                  CI_0.95 = CI2,
                  CI_0.99 = CI3))

 }


