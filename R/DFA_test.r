#' @title Statistical test for Detrended Fluctuation Analysis.
#'
#' @description This function performs the statistical test for the long-range correlation exponents obtained by the Detrended Fluctuation Analysis method.
#'
#' @details This function include following measures alpha_dfa, se_alpha_dfa, r2_alpha_dfa, min_test, max_test, mean_test, median_test, sd_test, skewness_test, kurtosis_test, Jarquebera_test_pvalue, CL_lower_test, CL_upper_test
#'
#' @param y A vector contaning univariate time series.
#'
#' @param npoints The number of different window sizes that will be used to estimate the Fluctuation function in each zone. See nonlinearTseries package.
#'
#' @param rep An integer value indicating the number of repetitions.
#'
#' @param ts.sim An logical value. If TRUE, the confidence interval for alpha_dfa is obtained from a White Gaussian Noise. If FALSE, the confidence interval for alpha_dfa is obtained from the shuffling of the original series.
#'
#' @param prob An numeric value indicating the quantile of probability to be used in estimating confidence intervals by N(0,1).
#'
#' @return An rbind matrix containing "alpha_dfa","se_alpha_dfa", "r2_alpha_dfa","min_alpha_dfa","max_test","mean_test", "median_test", "sd_test", "skewness_test", "kurtosis_test", "jarquebera_test_pvalue", and confidence interval: "CI_lower_test", "CI_upper_test".
#'
#' @examples
#' y=rnorm(1000)
#'dfa.test(y, npoints=15, rep=10,ts.sim="TRUE", prob=.95)
#'
#' @references
#' KRISTOUFEK, L. Rescaled Range Analysis and Detrended Fluctuation Analysis: Finite Sample Properties and Confidence Intervals. AUCO Czech Economic Review, v.4,n.3, p.315-329, 2010.
#'
#' @importFrom fgpt fyshuffle
#' @importFrom nonlinearTseries dfa
#' @importFrom stats coef lm rnorm qnorm sd
#' @importFrom PerformanceAnalytics skewness kurtosis
#' @importFrom tseries jarque.bera.test
#'
#' @export
dfa.test <- function(y, npoints, rep, ts.sim, prob){

  if(!(is.null(y) || is.numeric(y) || is.logical(y))){
    stop("Time series must be numeric")
  }

  alpha_dfa <- c()
  se_alpha_dfa <- c()
  r2_alpha_dfa <- c()

  if(rep > 30){
      error <- stats::qnorm(c(prob+((1-prob)/2)), 0,1)
  }
  if(rep > 0 & rep <= 30){
      error <- stats::qnorm(c(prob+((1-prob)/2)), 0,1)/rep
  }

  if(ts.sim == 'TRUE'){
    m <- matrix(data=NA, nrow=length(y), ncol=rep, byrow=F)
    for(i in 1:rep){
      m[,i] <- stats::rnorm(length(y), mean=0, sd=1)
    }
    m <- cbind(y, m)
    for(i in 1:ncol(m)){
      dfa <- nonlinearTseries::dfa(m[,i],
                                   window.size.range=c(4,round(length(y)/4,0)),
                                   npoints=npoints,
                                   do.plot=FALSE)

      model <- stats::lm(log10(dfa$fluctuation.function)~log10(dfa$window.sizes))
      alpha_dfa[i] <- stats::coef(summary(model))[2, "Estimate"]
      se_alpha_dfa[i] <- stats::coef(summary(model))[2, "Std. Error"]
      r2_alpha_dfa[i] <- summary(model)$r.squared
    }

    return(rbind(alpha_dfa = alpha_dfa[1],
                 se_alpha_dfa = se_alpha_dfa[1],
                 r2_alpha_dfa = r2_alpha_dfa[1],
                 min_test = min(alpha_dfa[2:length(alpha_dfa)]),
                 max_test = max(alpha_dfa[2:length(alpha_dfa)]),
                 mean_test = mean(alpha_dfa[2:length(alpha_dfa)]),
                 median_test = stats::median(alpha_dfa[2:length(alpha_dfa)]),
                 sd_test = stats::sd(alpha_dfa[2:length(alpha_dfa)]),
                 skewness_test = PerformanceAnalytics::skewness(alpha_dfa[2:length(alpha_dfa)], method="moment"),
                 kurtosis_test = PerformanceAnalytics::kurtosis(alpha_dfa[2:length(alpha_dfa)], method="moment"),
                 Jarquebera_test_pvalue = tseries::jarque.bera.test(alpha_dfa[2:length(alpha_dfa)])$p.value,
                 CI_lower_test = mean(alpha_dfa[2:length(alpha_dfa)]) - error*stats::sd(alpha_dfa[2:length(alpha_dfa)]),
                 CI_upper_test = mean(alpha_dfa[2:length(alpha_dfa)]) + error*stats::sd(alpha_dfa[2:length(alpha_dfa)])))
   }

         if(ts.sim == 'FALSE'){
           m <- matrix(data=NA, nrow=length(y), ncol=rep, byrow=F)
           for(i in 1:rep){
	           m[,i] <- fgpt::fyshuffle(1:length(y))
           }
		        m <- cbind(y=1:length(y), m)
		        for(i in 1:ncol(m)){
                  dfa <- nonlinearTseries::dfa(y[m[,i]],
                                               window.size.range=c(4,round(length(y)/4,0)),
                                               npoints=npoints,
                                               do.plot=FALSE)

                model <- stats::lm(log10(dfa$fluctuation.function)~log10(dfa$window.sizes))
         alpha_dfa[i] <- stats::coef(summary(model))[2, "Estimate"]
      se_alpha_dfa[i] <- stats::coef(summary(model))[2, "Std. Error"]
      r2_alpha_dfa[i] <- summary(model)$r.squared
          }

		        return(rbind(alpha_dfa = alpha_dfa[1],
		                     se_alpha_dfa = se_alpha_dfa[1],
		                     r2_alpha_dfa = r2_alpha_dfa[1],
		                     min_test = min(alpha_dfa[2:length(alpha_dfa)]),
		                     max_test = max(alpha_dfa[2:length(alpha_dfa)]),
		                     mean_test = mean(alpha_dfa[2:length(alpha_dfa)]),
		                     median_test = stats::median(alpha_dfa[2:length(alpha_dfa)]),
		                     sd_test = stats::sd(alpha_dfa[2:length(alpha_dfa)]),
		                     skewness_test = PerformanceAnalytics::skewness(alpha_dfa[2:length(alpha_dfa)], method="moment"),
		                     kurtosis_test = PerformanceAnalytics::kurtosis(alpha_dfa[2:length(alpha_dfa)], method="moment"),
		                     Jarquebera_test_pvalue = tseries::jarque.bera.test(alpha_dfa[2:length(alpha_dfa)])$p.value,
		                     CI_lower_test = mean(alpha_dfa[2:length(alpha_dfa)]) - error*stats::sd(alpha_dfa[2:length(alpha_dfa)]),
		                     CI_upper_test = mean(alpha_dfa[2:length(alpha_dfa)]) + error*stats::sd(alpha_dfa[2:length(alpha_dfa)])))

      }
}

