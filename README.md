## Installing the `AsyncCluster` package

You may download the package by using the following command.

    devtools::install_github("JasonCho0914/AsyncCluster")

If you also want to see this vignette on your machine, (this will take
some time)

    devtools::install_github("JasonCho0914/AsyncCluster",
                             build_vignettes = TRUE)

The following packages were used for the analysis.

    # install.packages(c("quantmod","foreach","ggplot2))
    library(quantmod)
    library(foreach)
    library(ggplot2)
    library(AsyncCluster)

## Introduction

This package `AsyncCluster` is designed to 1) interpolate missing data,
2) compute distances and 3) perform clustering on a sparsely sampled
time series. The main contribution of this work is the introduction of
statistical distance as a distance metric on interpolated data points.
As you will see in the coming example, using statistical distance as a
distance metric greatly increased the accuracy of clustering for
sparsely sampled time-series, when we were dealing with stationary time
series with varying level of volatility in each cluster.

## Data Generation

For demonstration, daily log returns based on the adjusted closing price
of 40 companies from January 1st, 2022 to December 4th, 2022 were
gathered using `quantmod` package. 20 of the companies were big-tech
companies, such as Apple, Microsoft, Google, Oracle, and etc. The other
20 companies were companies that have shown high volatility over the
past year, most of which are small pharmaceutical companies. The
following code chunk gathers and pre-proecesses the data. The variable
`list_df` is a list of data frames. The first data frame contains the
log prices of the small pharmaceutical companies and the second data
frame contains the log prices of the big tech companies.Each column of
the data frame represents the time series of one column. (Each data
frame is n by 20, where n represents the number of days, and 20 is the
number of companies),


    # Software cos
    software_company_t <- c("AAPL","MSFT","GOOG","ORCL","ADBE","CRM",
                            "SAP","IBM","INTU","ADP","NOW",
                            "SNPS","PANW","META","SNOW",
                            "CDNS","ENPH","WDAY","AMZN","NVDA")
    getSymbols(software_company_t)
    #>  [1] "AAPL" "MSFT" "GOOG" "ORCL" "ADBE" "CRM"  "SAP"  "IBM"  "INTU" "ADP" 
    #> [11] "NOW"  "SNPS" "PANW" "META" "SNOW" "CDNS" "ENPH" "WDAY" "AMZN" "NVDA"
    cluster1 <- list(AAPL,MSFT,GOOG,ORCL,ADBE,CRM,SAP,IBM,INTU,ADP,NOW,
                     SNPS,PANW,META,SNOW,CDNS,ENPH,WDAY,AMZN,NVDA)

    # High volatility cos
    High_volatility <-c("ELOX","BIVI","CTIB","TOPS","WBEV","RMED","POAI","COSM",
                        "GLS","IVC","VIVK","AVXL","CRBP","RIGL","MEIP","VAPO","MTC",
                        "TPHS","EQRx","BMY")
    getSymbols(High_volatility)
    #>  [1] "ELOX" "BIVI" "CTIB" "TOPS" "WBEV" "RMED" "POAI" "COSM" "GLS"  "IVC" 
    #> [11] "VIVK" "AVXL" "CRBP" "RIGL" "MEIP" "VAPO" "MTC"  "TPHS" "EQRx" "BMY"
    cluster2 <- list(ELOX,BIVI,CTIB,TOPS,WBEV,RMED,POAI,COSM,
                     GLS,IVC,VIVK,AVXL,CRBP,RIGL,MEIP,VAPO,MTC,TPHS,
                     EQRX,BMY)

    sw_returns<- foreach(df = cluster1,.combine = 'cbind')%do%{
      df <- as.data.frame(df)
      dates <- as.Date(rownames(df))
      after_2022 <- (dates>"2022-01-01") & (dates < "2022-12-04")
      prices <- df[after_2022,4]
      returns <- log(prices[-1]/prices[1:(length(prices)-1)])
      return(returns)
    }

    hv_returns<- foreach(df = cluster2,.combine = 'cbind')%do%{
      df <- as.data.frame(df)
      dates <- as.Date(rownames(df))
      after_2022 <- (dates>"2022-01-01") & (dates < "2022-12-04")
      prices <- df[after_2022,4]
      returns <- log(prices[-1]/prices[1:(length(prices)-1)])
      return(returns)
    }
    list_df <- list(hv_returns,sw_returns)

The following is the visualization of the two clusters of companies, big
techs and small pharmas. We can see clear distinction between two
clusters in terms of the movements of their log returns, where the
volatility of the small pharmas are much greater than the big techs.

    i = 1
    temp_list <- list()
    for(df in list_df){
      t <- data.frame(index = seq(1,nrow(df)),df)
      t <- reshape::melt(t,id = "index")
      levels(t$variable) <- paste0(letters[i],1:length(levels(t$variable)))
      t$cluster <- i
      temp_list[[i]] <- t
      i = i+1
    }
    df_final <- do.call('rbind',temp_list)
    df_final$cluster <- as.factor(df_final$cluster)
    p <- ggplot(df_final, aes(x=index, y=value,group =variable,col = cluster))+
      geom_line()+
      guides(color=guide_legend(title="Company Type"))+
      ggtitle("Large Tech Companies vs Small Pharmaceuticals")+
      scale_color_manual(labels = c("Small Pharma", "Big Tech"), 
                         values =c("#F8766D", "#619CFF"))+
      xlab("Time")+
      ylab("Log(return)")+
      theme(plot.title = element_text(size=18),
            axis.title = element_text(size = 15),
            axis.text.x = element_text(size = 12,angle = 45, hjust = 1),
            axis.text.y = element_text(size = 12),
            legend.title = element_text(size=14),
            legend.text = element_text(size=10))
    plot(p)

![](README_files/figure-markdown_strict/visualizaiton-1.png)

`AsyncCluster` is used to cluster time series with missing values or
asynchronous time series. Data points were removed at random for
demonstration. on average 5% of the data points, not including the first
and the last data points were sampled from each time series. (Since this
is `interpolation` the first and the last data points need to exist)

    # remove data points at random
    set.seed(4000)
    rm_datapoints <- function(df_timeseries){
      ret <- list()
      j = 1
      for(df in df_timeseries){
        len = nrow(df)
        for(i in 1:ncol(df)){
          n_samples = stats::rbinom(n=1,size = len,prob = 0.05)
          picked <- sample(x = 2:(len-1),size = n_samples,replace = F)
          df[!((1:len) %in% c(1,picked,len)),i] = NA
        }
        ret[[j]] = df
        j = j+1
      }
      return(ret)
    }
    list_df_rm <- rm_datapoints(list_df)

We have two clusters of sparesely sampled time series, which we will
apply clustering on. Our goal is to see how well our clustering
algorithm is able to distinguish different clusters.

## Interpolation

`AsyncCluster` allows you to either perform linear interpolation through
`linear_interpolation` or interpolation through Gaussian Process
regression (GP interpolation) using `gp_interpolation`. The parameter
`n_intervals` is used to rescale the data into a regular interval.
`list_df_rm` were combined as one data frame and the interpolations were
made on a rescaled 100 intervals. The function outputs a list, where
`$pred` represents the dataframe representing the interpolated data
points, and `$se` represents the standard error (only for the gp
output). For comparison, we visualized one of the path:

    list_df_rm_lm <- AsyncCluster::linear_interpolation(do.call("cbind",list_df_rm),n_intervals = 100)
    list_df_rm_gp <- AsyncCluster::gp_interpolation(do.call("cbind",list_df_rm),n_intervals = 100)
    sample_path = list_df_rm[[1]][,1]
    idx = (1:length(sample_path))[!is.na(sample_path)]
    plot(y = sample_path[idx],
         x = idx,pch =20, cex = 3, 
         ylab = "log(return)", xlab = "time", main = "GP vs Linear Interpolation")
    new_intervals = seq(from = 1,to = length(sample_path),length.out = 100)
    lines(x=new_intervals,y = list_df_rm_lm$pred[,1],col = "red")
    lines(x=new_intervals,y = list_df_rm_gp$pred[,1],col = "blue")
    legend("bottomright", lty = 1, col =c("blue","red"),legend = c("GP","Linear"))

![](README_files/figure-markdown_strict/unnamed-chunk-3-1.png)

## Computing distances

We may compute pairwise distance matrix using `compute_dm()` method. You
can also apply dynamic time warping, by setting `DTW = T`. For
comparison, distance matrix based on the squared hellinger distance, on
GP interpolated data points with `DTW = F` and one with the L1 distance
on a linear interpolated data points with `DTW = F` were computed.

    # DTW =F, and hellinger distance on GP interpolated points
    dm_gp_hellinger <- compute_dm(M_pred = list_df_rm_gp$pred,M_se = list_df_rm_gp$se,DTW = F, type = "hellinger")
    #  DTW =F, and L1 distance on a linearly interpolated points
    dm_lin_l1 <- compute_dm(M_pred = list_df_rm_lm$pred,type = "l1")

    dm_lin_l1_DTW = compute_dm(list_df_rm_lm$pred,
                                DTW=T,
                                type = "l1")

## Cluster analysis

You may apply any clustering algorithm that uses distance matrix.
`AsyncCluster` currently provides `k_medoid` for the k-medoid and
`kmeans_ts()` for k-means algorithm. 4 comparisons were made in terms of
accuracy:

    # we know that the first 20 are from the same cluster and the last 20 were from the other clusters
    # Thus we have
    True_labels = c(rep(1,20),rep(2,20)) 
    accuracy <- function(prediction,gt){
      n = length(gt)
      return(max(sum(prediction == gt)/n, sum(!(prediction == gt))/n))
    }

    gp_medoid <- k_medoid(dm_gp_hellinger,k = 2)
    lin_means <- kmeans_ts(list_df_rm_lm$pred,k = 2,DTW = F)
    lin_medoid_dtw <- k_medoid(dm_lin_l1_DTW,k= 2)
    lin_medoid <- k_medoid(dm_lin_l1,k= 2)

    labels = c("GP\nHellinger\nK-medoid\nNo DTW",
                    "Linear\nL1\nK-means\nNo DTW",
                    "Linear\nL1\nK-medoid\nNo DTW",
                    "Linear\nL1\nK-medoid\nDTW")
    p <- barplot(c(accuracy(gp_medoid,True_labels),
              accuracy(lin_means,True_labels),
              accuracy(lin_medoid_dtw,True_labels),
              accuracy(lin_medoid,True_labels))*100,
            ylim = c(0,100),
            ylab = "Accuracy (%)",
            xaxt ='n',
            main = "Comaprisons between different algorithms")
    axis(1,           
         las=1, padj = 0.65,
         at=p,          
         lab=labels)

![](README_files/figure-markdown_strict/unnamed-chunk-5-1.png)

As the model suggests, K-medoid algorithm using the squared Hellinger
distance as a distance metric based on the estimated standard errors
from Gaussian Regression had the highest accuracy.
