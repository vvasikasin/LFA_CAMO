#--------------------------------------------------------------------------------------#
# R script used in 2022 paper on xxxxx 
# Script author - VV
# 19/10/2022

#--------------------------------------------------------------------------------------#

#--------------------------------------------------------------------------------------#
####                                  LIBRARIES                                     ####
#--------------------------------------------------------------------------------------#

library(EBImage)
library(dplyr)
library(hyperSpec)
library(stats)
library(baseline)

#--------------------------------------------------------------------------------------#
####                                   FUNCTIONS                                    ####
#--------------------------------------------------------------------------------------#

readLFA <- function(Imagefile){
  image <- readImage(files = Imagefile)
  data <- image@.Data
  
  #make grayscale
  data <- apply(image, c(1,2), mean)
  data <- data.frame(data)
  
  #remove artifacts by excluding <5% and >95% of intensity of each row
  data <- t(apply(data, 2, function(x)
    replace(x, x>quantile(x, probs = 0.95), NA)))
  data <- t(apply(data, 1, function(x)
    replace(x, x<quantile(x, probs = 0.05, na.rm = TRUE), NA)))
  
  #mean of each row
  data <- t(data)
  cd <- colMeans(data, na.rm = TRUE)
  cd <- data.frame(cd)
  
  ##extracting middle 5 row
  midrow <- cd %>% 
    head(n = -((nrow(cd)/2)-10)) %>% 
    tail(n = 20) %>%
    apply(1, median) %>% 
    median()
  
  ##extracting first 5 row
  firstrow <- cd %>% 
    head(n = 20) %>%
    apply(1, median) %>% 
    median()
 
  ##extracting last 5 row
  lastrow <- cd %>% 
    tail(n=20) %>%
    apply(1, median) %>% 
    median()
  
  ##normalise baseline with slope of each segment
  slopefirst = (midrow-firstrow)*2/nrow(cd)
  slopelast = (lastrow-midrow)*2/nrow(cd)
  slopelist <- c(
    ((0:(nrow(cd)/2)*slopefirst)+((nrow(cd)/2)*(slopelast-slopefirst))), 
    ((nrow(cd)/2):nrow(cd)*slopelast))
  cd2 <- cd-slopelist
  
  ##intensity signal noise reduction 
  Intensity <- 1/rowMeans(cd2)
  mov_avg_7 <- rep(1/10, 10)
  intensity2 <- stats::filter(Intensity, mov_avg_7)
  intensity2 <- data.frame(intensity2)
  intensity2[1:4,] <- intensity2[5,] 
  intensity2[(nrow(intensity2)-4):nrow(intensity2),] <- intensity2[nrow(intensity2)-5,] 
  intensity3 <- t(intensity2)
  bs.inten <- baseline::baseline.als(intensity3, lambda = 5, p = 0.01, maxit = 20)
  intensity4 <- bs.inten$corrected

  print(paste ('Max intensity first  =', 
               max(intensity4[1,(1:(ncol(intensity4)/2))])))
  print(paste ('Max intensity second  =', 
               max(intensity4[1,((ncol(intensity4)/2):(ncol(intensity4)))])))
  print(paste ('Max intensity ratio  =', 
               max(intensity4[1,(1:(ncol(intensity4)/2))]) / max(intensity4[1,((ncol(intensity4)/2):(ncol(intensity4)))])))
  
  string_vector <- c(
                    (max(intensity4[1,(1:(ncol(intensity4)/2))])), 
                    (max(intensity4[1,((ncol(intensity4)/2):(ncol(intensity4)))])), 
                    (max(intensity4[1,(1:(ncol(intensity4)/2))]) / max(intensity4[1,((ncol(intensity4)/2):(ncol(intensity4)))])),
                     "")
  paste(string_vector, collapse = '\n') %>% cat()
  
  plot(x = seq(1, nrow(intensity2), 1), 
       y = intensity4, type = "l", 
       lwd = 2, col = "blue", 
       xlab = "Pixel", ylab = "Intensity")
  }


readLFA("R/Test/T1.jpg")
