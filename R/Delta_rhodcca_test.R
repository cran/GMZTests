#' @title Statistical test for Statistical test for Delta RHODCCA cross-correlation coefficient.
#'
#' @description This function performs the statistical test for Delta RHODCCA cross-correlation coefficient from two univariate ARFIMA process.
#'
#' @details This function include following measures: timescale, rho_before, rho_affter, deltarho
#'
#' @param x A vector containing univariate time series.
#'
#' @param y A vector containing univariate time series.
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
#' @return An list containing "timescale", "mean", "sd", "rho_before", "rho_affter", "deltarho", "CI_0.90", "CI_0.95", "CI_0.99".
#'
#' @examples
#' x <- rnorm(1000)
#' y <-  rnorm(1000)
#' deltarhodcca.test(x,y,k=100,m=c(4:6),nu=0,rep=10)
#'
#' @references
#' Guedes, et al. Statistical test for DCCA cross-correlation coefficient, Physica A, v.501, 134-140, 2018.
#'
#' Guedes, et al. Statistical test for DCCA: Methods and data, Data in Brief, v. 18, 795-798, 2018.
#'
#' @importFrom DCCA rhodcca
#' @importFrom stats rnorm qnorm sd
#' @importFrom fgpt fyshuffle
#'
#' @export
deltarhodcca.test <- function(x,y,k,m,nu,rep){

  if(!(is.null(y) || is.numeric(y) || is.logical(y))){
    stop("Time series must be numeric")
  }
  if(!(is.null(x) || is.numeric(x) || is.logical(x))){
    stop("Time series must be numeric")
  }
  if(length(x) != length(y)){
    stop("Time series have different lengths")
  }
    N <- length(y)
   N1 <- length(y)/2
   N2 <- N1+1

   rho_a <- matrix(data = NA, nrow = rep, ncol = length(m), byrow = TRUE)
   rho_b <- matrix(data = NA, nrow = rep, ncol = length(m), byrow = TRUE)
deltarho <- matrix(data = NA, nrow = rep, ncol = length(m), byrow = TRUE)

     for(i in 1:rep){

         shuf1 <- fgpt::fyshuffle(1:N1)
         shuf2 <- fgpt::fyshuffle(1:N1)
         shuf3 <- fgpt::fyshuffle(1:N1)
         shuf4 <- fgpt::fyshuffle(1:N1)


            xa <- x[1:N1]
            xb <- x[N2:N]
            ya <- y[1:N1]
            yb <- y[N2:N]

            for(j in 1:length(m)){
            rho_a[i,j] <- DCCA::rhodcca(xa[shuf1], ya[shuf2], m[j], nu = nu)$rhodcca
            rho_b[i,j] <- DCCA::rhodcca(xb[shuf3], yb[shuf4], m[j], nu = nu)$rhodcca
         deltarho[i,j] <- (rho_b[i,j] - rho_a[i,j])
      	   }
         }

            mean <- apply(deltarho, 2, mean)
              sd <- apply(deltarho, 2, sd)

             CI1 <- stats::qnorm(c(0.90+((1-0.90)/2)), 0,1)*sd
             CI2 <- stats::qnorm(c(0.95+((1-0.95)/2)), 0,1)*sd
             CI3 <- stats::qnorm(c(0.99+((1-0.99)/2)), 0,1)*sd


   return(list(timescale=m,
               mean = mean,
               sd = sd,
               rho_before = rho_a,
               rho_after = rho_b,
               deltarho=deltarho,
               CI_0.90 = CI1,
               CI_0.95 = CI2,
               CI_0.99 = CI3))
}
