#' squared Hellinger distances between two normal distributions
#'
#' `hellinger` computes and returns squared Hellinger distance between two normal distributions
#'
#' @param mu1 mean of the first normal distribution
#' @param se1 standard error of the first normal distribution
#' @param mu2 mean of the second normal distribution
#' @param se2 standard error of the second normal distribution
#' @details The squared Hellinger distance between two normal distributions, \eqn{X \sim N(\mu_x,\sigma_x^2)} and \eqn{Y \sim N(\mu_y,\sigma_y^2)} \deqn{H^2(X,Y) = 1-\sqrt{\frac{2\sigma_x\sigma_y}{\sigma_x^2+\sigma_2^2}}e^{\frac{-(\mu_x -\mu_y)^2}{4*(\sigma_x^2 + \sigma_y^2)}}}
#' @export
hellinger <- function(mu1,se1,mu2,se2){
  se1[se1==0] = 1.0e-09
  se2[se2==0] = 1.0e-09
  computed <- 1-sqrt((2*se1*se2)/(se1^2 + se2^2))*exp(-(mu1-mu2)^2/(4*(se1^2 + se2^2)))
  return(computed)
}

#' Calculate Chi-squared statistics based on two Normal Distributions
#'
#' `z2_score` computes and returns Chi-squared statistics based on two normal distributions
#'
#' @param mu1 mean of the first normal distribution
#' @param se1 standard error of the first normal distribution
#' @param mu2 mean of the second normal distribution
#' @param se2 standard error of the second normal distribution
#' @details Given \eqn{X \sim N(\mu_x,\sigma_x^2)} and \eqn{Y \sim N(\mu_y,\sigma_y^2)},
#' I calculate \deqn{(\frac{(\mu_x - \mu_y)}{(\sigma_x^2 +\sigma_y^2)})^2}
#' @export
z2_score <- function(mu1,se1,mu2,se2){
  se1[se1==0] = 1.0e-09
  se2[se2==0] = 1.0e-09
  computed <- ((mu1-mu2)/(sqrt(se1^2 + se2^2)))^2
  return(computed)
}

#' Apply k-means clustering to time series
#'
#' `kmeans` apply k-means clustering algorithm to a collection of time series and return a vector representing clusters
#'
#' @param df_timeseries a data.frame where each column represents time series.
#' @param k number of clusters
#' @param DTW whether to apply dynamic time-warping or not
#' @details when DTW is not applied `kmeans` function from the stat package is used. when DTW is applied, `dtwclust` function from tsclust package is used.
#' @export
kmeans_ts <- function(df_timeseries,k,DTW = F){
  if(DTW == F){
    predicted_cluster <- stats::kmeans(t(df_timeseries),k)$cluster
  }else{
    predicted_cluster <- dtwclust::tsclust(t(df_timeseries),k = 2,centroid = "mean")@cluster
  }
  return(predicted_cluster)
}

#' Computes distances between time-series
#'
#' `compute_dm` calculates pairwise distance between time series and returns a distance matrix.
#'
#' @param M_pred a data.frame of time-series, where each column represents a time-series
#' @param M_se a data.frame of standard errors representing the level of confidence for each interpolated data points on M_pred. use if standard error is calculated for M_pred.
#' @param DTW whether to apply dynamic time-warping or not
#' @param type either "l1" vs "hellinger", l1 distance or hellinger distance is applied
#' @export
compute_dm <- function(M_pred,M_se=NULL,DTW = F,type = "l1"){
  if (DTW == T){
    if (type == "l1"){
      distances = DTW_matrix(M_pred,type = "l1")
    }else{
      distances = DTW_matrix(M_pred,M_se,type)
    }
  }else{
    if (type == "l1"){
      distances = eucl_matrix(M_pred,type = "l1")
    }else{
      distances = eucl_matrix(M_pred,M_se,type)
    }
  }
  return(stats::as.dist(distances))
}

#' Apply k-medoid clustering to time series
#'
#' `kmedoid` apply k-medoid clustering algorithm based on a distance matrix
#'
#' @param dm is a matrix or distance matrix object
#' @param k number of clusters
#' @details this function uses `fastkmed` from kmed package.
#' @export
k_medoid <- function(dm,k){
  return(kmed::fastkmed(dm,k)$cluster)
}

#' Apply agglomerative clustering to time series
#'
#' `agglomerative` apply complete linkage agglomerative clustering algorithm based on a distance matrix
#'
#' @param dm is a matrix or distance matrix object
#' @param k number of clusters
#' @details this function uses `hclust` and `cutree` from kmed package. The complete linkage is used.
#' @export
agglomerative <- function(dm,k){
  hc1 <- stats::hclust(dm, method = "complete" )
  return(stats::cutree(hc1,k))
}


DTW_vec <- function(predA, seA = NULL, predB, seB = NULL,type = "l1"){
  n1 <- length(predA)
  n2 <- length(predB)
  DTW <- matrix(rep(Inf,(n1+1)*(n2+1)),nrow = (n1+1))
  DTW[1, 1] = 0
  for (i in 1:n1){
    for (j in 1:n2){
      if (type == "l1"){
        cost = abs(predA[i]-predB[j])
      }else if(type == "hellinger"){
        cost = hellinger(predA[i],seA[i],predB[j],seB[j])
      }else if(type == "z2"){
        cost = z2_score(predA[i],seA[i],predB[j],seB[j])
      }
      DTW[i+1,j+1] =  cost + min(DTW[i,j+1],
                                 DTW[i+1,j],
                                 DTW[i,j])
    }
  }
  return(DTW[n1+1,n2+1])
}


DTW_matrix <- function(M_pred,M_se=NULL,type = "l1"){
  n = ncol(M_pred)
  distances <- matrix(rep(0,n*n),nrow = n)
  for (i in 1:n){
    for (j in i:n){
      if (type == "l1"){
        d = DTW_vec(predA = M_pred[,i],predB = M_pred[,j])
      }else{
        d = DTW_vec(M_pred[,i],M_se[,i],M_pred[,j],M_se[,j],type = type)
      }
      distances[i,j] = d
      if (i != j){
        distances[j,i] = d
      }
    }
  }
  return(distances)
}

#euclidean matching
eucl_matrix <- function(M_pred,M_se=NULL,type="l1"){
  if(type =="l1"){
    distances = stats::dist(t(M_pred),method = "manhattan")
  }
  else{
    n = ncol(M_pred)
    distances <- matrix(rep(0,n*n),nrow = n)
    for (i in 1:n){
      for (j in i:n){
        if(type == "hellinger"){
          d = sum(hellinger(M_pred[,i],M_se[,i],
                            M_pred[,j],M_se[,j]))
        }else if(type == "z2"){
          d = sum(z2_score(M_pred[,i],M_se[,i],
                           M_pred[,j],M_se[,j]))
        }
        distances[i,j] = d
        if (i != j){
          distances[j,i] = d
        }
      }
    }
  }
  return(distances)
}



