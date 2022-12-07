#' linearly interpolate time series with missing value
#'
#' `linear_interpolation` is used to linearly interpolate time series with missing value.
#'
#' @param dataframe a data.frame with each column representing time-series with missing value as NA. The first and the last item must exist.
#' @param n_intervals the number of intervals to be interpolated. This is used to rescale your time series on a regular interval.
#' @details `approx` function from the stats package was used to linearly inteproalte data points.
linear_interpolation <- function(dataframe,n_intervals){
  ret <- list()
  last_idx <- length(dataframe[,1])# length of the time series
  prediction_interval = seq(1, last_idx, length.out = n_intervals)
  temp <- c()
  for (j in 1:ncol(dataframe)){
    vec <- dataframe[,j]
    temp <- c(temp,stats::approx(vec,xout = prediction_interval)$y)
  }
  df_interpolated <- matrix(temp,nrow = n_intervals,byrow = F)
  ret$pred <- df_interpolated
  return(ret)
}

#' Gaussian Process regression with Matern kernel is used to interpolate time series
#'
#' `gp_interpolation` returns a data.frame with interpolated data points based on Gaussian Process regression.
#'
#' @param dataframe a data.frame with each column representing time-series with missing value as NA.
#' @param n_intervals the number of intervals to be interpolated. This is used to rescale your time series on a regular interval.
#' @details `GauPro_kernel_model` and `Matern52` function from the GauPro package was used to inteproalte data points.
gp_interpolation <- function(dataframe,n_intervals){
  ret <- list()
  last_idx <- length(dataframe[,1])# last index
  prediction_interval = seq(1, last_idx, length.out = n_intervals)
  pred <- c()
  se <- c()
  for (j in 1:ncol(dataframe)){
    vec <- dataframe[,j]
    idx <- c(1:length(vec))[!is.na(vec)]
    vec_nna <- vec[!is.na(vec)]
    kern <- GauPro::Matern52$new(0) #matern Kernel
    fit <- GauPro::GauPro_kernel_model$new(matrix(idx, ncol=1),
                                   vec_nna, kernel=kern, parallel=FALSE)
    fit_pred <- fit$predict(prediction_interval,se.fit = T)
    pred <- c(pred,fit_pred$mean)
    se <- c(se,fit_pred$se)
  }
  df_pred <- matrix(pred,nrow = n_intervals,byrow = F)
  df_se <- matrix(se,nrow = n_intervals,byrow = F)
  ret$pred <- df_pred
  ret$se <- df_se
  return(ret)
}


