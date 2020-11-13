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
#' @param rep An integer value indicating the number of repetitions.
#'
#' @return An list containing "timescale","mean", "sd" and confidence interval: "CI_0.90", "CI_0.95", "CI_0.99".
#'
#' @examples
#' dmca.test(N=100, k=10, rep=10)
#'
#' @references
#' B. Podobnik, Z.-Q. Jiang, W.-X. Zhou, H. E. Stanley, Statistical tests for power-law cross-correlated processes, Phys. Rev.,E 84, 066118, 2011.
#'
#' @importFrom stats rnorm qnorm filter sd
#'
#' @export
 dmca.test <- function(N,k,rep){

     n <- 4:round(N/k,0)

    yx <- matrix(data = NA, nrow = rep, ncol = length(n), byrow = TRUE)

    dmca <- function(x,y,n){
      xx <- cumsum(x)
      yy <- cumsum(y)

      mm <- c(rep(1,n))/n
      mm_x <- stats::filter(xx,mm)
      mm_y <- stats::filter(yy,mm)

      F2_xy <- mean((xx-mm_x)[(1+floor(n/2)):(length(xx)-floor(n/2))]*(yy-mm_y)[(1+floor(n/2)):(length(yy)-floor(n/2))])
      F2_xx <- mean((xx-mm_x)[(1+floor(n/2)):(length(xx)-floor(n/2))]*(xx-mm_x)[(1+floor(n/2)):(length(xx)-floor(n/2))])
      F2_yy <- mean((yy-mm_y)[(1+floor(n/2)):(length(yy)-floor(n/2))]*(yy-mm_y)[(1+floor(n/2)):(length(yy)-floor(n/2))])

      rho <- F2_xy/sqrt(F2_xx*F2_yy)
      return(rho)
    }

    for(i in 1:rep){
     y <- stats::rnorm(N, mean=0, sd=1)
     x <- stats::rnorm(N, mean=0, sd=1)

      for(j in 1:length(n)){
        yx[i,j]  <- dmca(y, x, n=n[j])
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


