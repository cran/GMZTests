#' @title Statistical test for Delta DMC Multiple Detrended Cross-Correlation Coefficient
#'
#' @description This function performs the statistical test for Delta DMC cross-correlation coefficient from three univariate ARFIMA process.
#'
#' @details This function include following measures: timescale, dmc_before, dmc_after, deltadmc
#'
#' @param y A vector containing univariate time series.
#'
#' @param x1 A vector containing univariate time series.
#'
#' @param x2 A vector containing univariate time series.
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
#' @param method A character string indicating which correlation coefficient is to be used. If method = "rhodcca" the dmc coefficient is generated from the DCCA coefficient. If method = "dmca", the dmc coefficient is generated from the DMCA coefficient.
#'
#' @return An list containing "timescale" "dmc_before", "dmc_after", "deltadmc", "CI_0.90", "CI_0.95", "CI_0.99".
#'
#' @examples
#' x1 <- rnorm(100)
#' x2 <- rnorm(100)
#' y <-  rnorm(100)
#' deltadmc.test(x1,x2,y, k=10, m=c(4:6), nu=0, rep=10, method="rhodcca")
#'
#' deltadmc.test(x1,x2,y, k=10, m=c(4:6), nu=0, rep=10, method="dmca")

#' @references
#' ZEBENDE, G.F.; SILVA-FILHO, A.M. Detrended Multiple Cross-Correlation Coefficient. PHYSICA A, v.510, p.91-97, 2018.
#'
#' SILVA-FILHO,A.M; ZEBENDE,G.; CASTRO,A.P.; GUEDES,E. Statistical test for multiple detrended cross-correlation coefficient, Physica A, v.562, 125285, 2021.
#'
#' KRISTOUFEK, L. Detrending moving-average cross-correlation coefficient: Measuring cross-correlations between non-stationary series. PHYSICA A, v.406, p.169-175, 2014.
#'
#' @importFrom DCCA rhodcca
#' @importFrom stats rnorm qnorm filter sd
#' @importFrom fgpt fyshuffle
#'
#' @export
deltadmc.test <- function(x1,x2,y,k,m,nu,rep,method){

  Nx1 <- length(x1)
  Nx2 <- length(x2)
  Ny <- length(y)
  N1 <- Nx1/2
  N2 <- N1+1

  if(!(is.null(y) || is.numeric(y) || is.logical(y))){
    stop("Time series must be numeric")
  }
  if(!(is.null(x1) || is.numeric(x1) || is.logical(x1))){
    stop("Time series must be numeric")
  }
  if(!(is.null(x2) || is.numeric(x2) || is.logical(x2))){
    stop("Time series must be numeric")
  }
  if(Nx1 != Nx2){
    stop("Time series have different lengths")
  }
  if(Nx1 != Ny){
    stop("Time series have different lengths")
  }
  if(Nx2 != Ny){
    stop("Time series have different lengths")
  }

    yx1a <- matrix(data = NA, nrow = rep, ncol = length(m), byrow = TRUE)
    yx1d <- matrix(data = NA, nrow = rep, ncol = length(m), byrow = TRUE)
    yx2a <- matrix(data = NA, nrow = rep, ncol = length(m), byrow = TRUE)
    yx2d <- matrix(data = NA, nrow = rep, ncol = length(m), byrow = TRUE)
   x1x2a <- matrix(data = NA, nrow = rep, ncol = length(m), byrow = TRUE)
   x1x2d <- matrix(data = NA, nrow = rep, ncol = length(m), byrow = TRUE)
   dmcat <- matrix(data = NA, nrow = rep, ncol = length(m), byrow = TRUE)
   dmcdp <- matrix(data = NA, nrow = rep, ncol = length(m), byrow = TRUE)
deltadmc <- matrix(data = NA, nrow = rep, ncol = length(m), byrow = TRUE)

if(method =='rhodcca'){
   for(i in 1:rep){

     shuf1 <- fgpt::fyshuffle(1:N1)
     shuf2 <- fgpt::fyshuffle(1:N1)
     shuf3 <- fgpt::fyshuffle(1:N1)

       x1a <- x1[1:N1]
       x1d <- x1[N2:Nx1]

       x2a <- x2[1:N1]
       x2d <- x2[N2:Nx1]

       ya <- y[1:N1]
       yd <- y[N2:Nx1]

    for(j in 1:length(m)){
         yx1a[i,j] <- DCCA::rhodcca(ya[shuf3], x1a[shuf1], m[j], nu = nu, overlap = TRUE)$rhodcca #antes
         yx1d[i,j] <- DCCA::rhodcca(yd[shuf3], x1d[shuf1], m[j], nu = nu, overlap = TRUE)$rhodcca #depois

		     yx2a[i,j] <- DCCA::rhodcca(ya[shuf3], x2a[shuf2], m[j], nu = nu, overlap = TRUE)$rhodcca
		     yx2d[i,j] <- DCCA::rhodcca(yd[shuf3], x2d[shuf2], m[j], nu = nu, overlap = TRUE)$rhodcca

		    x1x2a[i,j] <- DCCA::rhodcca(x1a[shuf1], x2a[shuf2], m[j], nu = nu, overlap = TRUE)$rhodcca
		    x1x2d[i,j] <- DCCA::rhodcca(x1d[shuf1], x2d[shuf2], m[j], nu = nu, overlap = TRUE)$rhodcca

		    dmcat[i,j] <- (yx1a[i,j]^2 + yx2a[i,j]^2-(2*yx1a[i,j]*yx2a[i,j]*x1x2a[i,j]))/(1-x1x2a[i,j]^2)
		    dmcdp[i,j] <- (yx1d[i,j]^2 + yx2d[i,j]^2-(2*yx1d[i,j]*yx2d[i,j]*x1x2d[i,j]))/(1-x1x2d[i,j]^2)
		 deltadmc[i,j] <- (dmcdp[i,j] - dmcat[i,j])
		    }
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

      shuf1 <- fgpt::fyshuffle(1:N1)
      shuf2 <- fgpt::fyshuffle(1:N1)
      shuf3 <- fgpt::fyshuffle(1:N1)

      x1a <- x1[1:N1]
      x1d <- x1[N2:Nx1]

      x2a <- x2[1:N1]
      x2d <- x2[N2:Nx1]

      ya <- y[1:N1]
      yd <- y[N2:Nx1]

    for(j in 1:length(m)){
       yx1a[i,j] <- dmca(ya[shuf3], x1a[shuf1], m[j])
       yx1d[i,j] <- dmca(yd[shuf3], x1d[shuf1], m[j])

       yx2a[i,j] <- dmca(ya[shuf3], x2a[shuf2], m[j])
       yx2d[i,j] <- dmca(yd[shuf3], x2d[shuf2], m[j])

      x1x2a[i,j] <- dmca(x1a[shuf1], x2a[shuf2], m[j])
      x1x2d[i,j] <- dmca(x1d[shuf1], x2d[shuf2], m[j])

      dmcat[i,j] <- (yx1a[i,j]^2 + yx2a[i,j]^2-(2*yx1a[i,j]*yx2a[i,j]*x1x2a[i,j]))/(1-x1x2a[i,j]^2)
      dmcdp[i,j] <- (yx1d[i,j]^2 + yx2d[i,j]^2-(2*yx1d[i,j]*yx2d[i,j]*x1x2d[i,j]))/(1-x1x2d[i,j]^2)
   deltadmc[i,j] <- (dmcdp[i,j] - dmcat[i,j])
    }
  }
}

            mean <- apply(deltadmc, 2, mean)
              sd <- apply(deltadmc, 2, sd)

             CI1 <- stats::qnorm(c(0.90+((1-0.90)/2)), 0,1)*sd
             CI2 <- stats::qnorm(c(0.95+((1-0.95)/2)), 0,1)*sd
             CI3 <- stats::qnorm(c(0.99+((1-0.99)/2)), 0,1)*sd

   return(list(timescale=m,
               mean = mean,
               sd = sd,
               dmc_before = dmcat,
               dmc_after = dmcdp,
               deltadmc=deltadmc,
               CI_0.90 = CI1,
               CI_0.95 = CI2,
               CI_0.99 = CI3))

}
