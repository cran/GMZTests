#' @title Statistical test for DMCA cross-correlation coefficient.
#'
#' @description This function performs the statistical test for Detrending moving-average cross-correlation coefficient based in White Gaussian Noise process.
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
#' @param rep An integer value indicating the number of repetitions.
#'
#' @return An list containing "timescale","mean", "sd" and confidence interval: "CI_0.90", "CI_0.95", "CI_0.99".
#'
#' @examples
#' dmca.test(N=100, k=10, m=c(4:6), rep=10)
#'
#' @references
#' B. Podobnik, Z.-Q. Jiang, W.-X. Zhou, H. E. Stanley, Statistical tests for power-law cross-correlated processes, Phys. Rev.,E 84, 066118, 2011.
#'
#' @importFrom stats rnorm qnorm filter sd
#'
#' @export
 dmca.test <- function(N,k,m,rep){

    dmca <- function(x,y,m){
      xx <- cumsum(x)
      yy <- cumsum(y)

      mm <- c(rep(1,m))/m
      mm_x <- stats::filter(xx,mm)
      mm_y <- stats::filter(yy,mm)

      F2_xy <- mean((xx-mm_x)[(1+floor(m/2)):(length(xx)-floor(m/2))]*(yy-mm_y)[(1+floor(m/2)):(length(yy)-floor(m/2))])
      F2_xx <- mean((xx-mm_x)[(1+floor(m/2)):(length(xx)-floor(m/2))]*(xx-mm_x)[(1+floor(m/2)):(length(xx)-floor(m/2))])
      F2_yy <- mean((yy-mm_y)[(1+floor(m/2)):(length(yy)-floor(m/2))]*(yy-mm_y)[(1+floor(m/2)):(length(yy)-floor(m/2))])

      rho <- F2_xy/sqrt(F2_xx*F2_yy)
      return(rho)
    }

    yx <- matrix(data = NA, nrow = rep, ncol = length(m), byrow = TRUE)

    for(i in 1:rep){
     y <- stats::rnorm(N, mean=0, sd=1)
     x <- stats::rnorm(N, mean=0, sd=1)

      for(j in 1:length(m)){
        yx[i,j]  <- dmca(y, x, m[j])
       }
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


