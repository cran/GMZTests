#' @title Statistical test for Multiple Detrended Cross-Correlation Coefficient
#'
#' @description This function performs the statistical test for DMC Cross-Correlation Coefficient based in White Gaussian Noise process.
#'
#' @details This function include following measures: w, timescale, dmc, rhodcca_yx1, rhodcca_yx2, rhodcca_x1x2
#'
#' @param N An integer value for the time series length.
#'
#' @param k An integer value indicating the boundary of the division \eqn{(N/k)}.
#'          The smallest value of \eqn{k} is \eqn{4}.
#'
#' @param m an integer value or a vector of integer values indicating the size of the window for the polinomial fit.
#'
#' @param method A character string indicating which correlation coefficient is to be used. If method = "rhodcca" the dmc coefficient is generated from the DCCA coefficient. If method = "dmca", the dmc coefficient is generated from the DMCA coefficient.
#'
#' @param nu An integer value. See the DCCA package.
#'
#' @param rep An integer value indicating the number of repetitions.
#'
#' @return An list containing "timescale", parameters of beta distribution: "shape1", "se1","shape2","se2" and confidence interval: "CI_0.90_uppper", "CI_0.95_uppper", "CI_0.99_uppper".
#'
#' @examples
#' dmc.test(N=100, k=10, m=c(4:6), method="rhodcca", nu=0, rep=10)
#' dmc.test(N=100, k=10, m=c(4:6), method="dmca", nu=0, rep=10)
#'
#' @references
#' SILVA-FILHO,A.M; ZEBENDE,G.; CASTRO,A.P.; GUEDES,E. Statistical test for multiple detrended cross-correlation coefficient, Physica A, v.562, 125285, 2021.
#'
#' KRISTOUFEK, L. Detrending moving-average cross-correlation coefficient: Measuring cross-correlations between non-stationary series. PHYSICA A, v.406, p.169-175, 2014.
#'
#' @importFrom stats rnorm qbeta sd
#' @importFrom DCCA rhodcca
#' @importFrom fitdistrplus fitdist
#'
#' @export
 dmc.test <- function(N,k,m,method,nu,rep){

  yx1 <- matrix(data = NA, nrow = rep, ncol = length(m), byrow = TRUE)
  yx2 <- matrix(data = NA, nrow = rep, ncol = length(m), byrow = TRUE)
 x1x2 <- matrix(data = NA, nrow = rep, ncol = length(m), byrow = TRUE)
  dmc <- matrix(data = NA, nrow = rep, ncol = length(m), byrow = TRUE)

if(method =='rhodcca'){
for(i in 1:rep){
     y <- stats::rnorm(N, mean=0, sd=1)
    x1 <- stats::rnorm(N, mean=0, sd=1)
    x2 <- stats::rnorm(N, mean=0, sd=1)

        yx1[i,] <- DCCA::rhodcca(y,  x1, m, nu=nu)$rhodcca
        yx2[i,] <- DCCA::rhodcca(y,  x2, m, nu=nu)$rhodcca
       x1x2[i,] <- DCCA::rhodcca(x1, x2, m, nu=nu)$rhodcca
        dmc[i,] <- (yx1[i,]^2 + yx2[i,]^2-(2*yx1[i,]*yx2[i,]*x1x2[i,]))/(1-x1x2[i,]^2)
    }
}

  if(method =='dmca'){
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

    for(i in 1:rep){
       y <- stats::rnorm(N, mean=0, sd=1)
      x1 <- stats::rnorm(N, mean=0, sd=1)
      x2 <- stats::rnorm(N, mean=0, sd=1)

      for(j in 1:length(m)){
        yx1[i,j]  <- dmca(y,  x1, m[j])
        yx2[i,j]  <- dmca(y,  x2, m[j])
        x1x2[i,j] <- dmca(x1, x2, m[j])
        dmc[i,j]  <- (yx1[i,j]^2 + yx2[i,j]^2-(2*yx1[i,j]*yx2[i,j]*x1x2[i,j]))/(1-x1x2[i,j]^2)
      }
    }
  }

         shape1 <- double()
      ep_shape1 <- double()
         shape2 <- double()
      ep_shape2 <- double()
            CI1 <- double()
            CI2 <- double()
            CI3 <- double()

    for(i in 1:ncol(dmc)){

              model <- fitdistrplus::fitdist(dmc[,i], distr = "beta", method = "mle")
          shape1[i] <- model$estimate[1]
       ep_shape1[i] <- model$sd[1]
          shape2[i] <- model$estimate[2]
       ep_shape2[i] <- model$sd[2]

        CI1[i] <- stats::qbeta(.90, shape1=shape1[i], shape2=shape2[i], ncp = 0, lower.tail = TRUE, log.p = FALSE)
        CI2[i] <- stats::qbeta(.95, shape1=shape1[i], shape2=shape2[i], ncp = 0, lower.tail = TRUE, log.p = FALSE)
        CI3[i] <- stats::qbeta(.99, shape1=shape1[i], shape2=shape2[i], ncp = 0, lower.tail = TRUE, log.p = FALSE)
    }

     return(list(timescale=m,
                  shape1=shape1,
                  se_shape1=ep_shape1,
                  shape2=shape2,
                  se_shape2 = ep_shape2,
                  CI_0.90_uppper = CI1,
                  CI_0.95_uppper = CI2,
                  CI_0.99_uppper = CI3))
 }
