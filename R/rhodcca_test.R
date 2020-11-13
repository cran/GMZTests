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
#' @param nu An integer value. See the DCCA package.
#'
#' @param rep An integer value indicating the number of repetitions.
#'
#' @return An list containing "timescale","mean", "sd" and confidence interval: "CI_0.90", "CI_0.95", "CI_0.99".
#'
#' @examples
#' rhodcca.test(N=100, k=10, nu=0, rep=10)
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
 rhodcca.test <- function(N,k,nu,rep){

     n <- 4:round(N/k,0)

    yx <- matrix(data = NA, nrow = rep, ncol = length(n), byrow = TRUE)

    for(i in 1:rep){
     y <- stats::rnorm(N, mean=0, sd=1)
     x <- stats::rnorm(N, mean=0, sd=1)

      for(j in 1:length(n)){
        yx[i,j]  <- DCCA::rhodcca(y, x, m=n[j], nu=nu)$rhodcca
       }
     }

           m <- double()
           s <- double()
         CI1 <- double()
         CI2 <- double()
         CI3 <- double()

    for(i in 1:ncol(yx)){

          m[i] <- mean(yx[,i])
          s[i] <- stats::sd((yx[,i]))

        CI1[i] <- stats::qnorm(c(0.90+((1-0.90)/2)), 0,1)*s[i]
        CI2[i] <- stats::qnorm(c(0.95+((1-0.95)/2)), 0,1)*s[i]
        CI3[i] <- stats::qnorm(c(0.99+((1-0.99)/2)), 0,1)*s[i]
    }

     return(list(timescale=n,
                  mean = m,
                  sd = s,
                  CI_0.90 = CI1,
                  CI_0.95 = CI2,
                  CI_0.99 = CI3))

 }


