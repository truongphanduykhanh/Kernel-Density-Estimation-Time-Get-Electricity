# ========================================================================#
#                       Non Parametric - Project
#                   Huong Li Nguyen | Khanh Truong
# ========================================================================#

# Import libraries --------------------------------------------------------
library(KernSmooth)
library(ggplot2)
library(tidyverse)
library(kedd)
library(np)


# Import and clean data  --------------------------------------------------
temp <- tempfile()
url<-
  'http://api.worldbank.org/v2/en/indicator/IC.ELC.TIME?downloadformat=csv'
download.file(url,temp, mode="wb")
unzip(temp, "API_IC.ELC.TIME_DS2_en_csv_v2_10187384.csv")
data <- read.csv("API_IC.ELC.TIME_DS2_en_csv_v2_10187384.csv", skip=4, 
                 check.names = F)
data <- data[c('Country Name','Country Code', 2009:2017)] # Select columns

# Clean to get data only for countries, not continent
country <- data[
  (is.na(data$`2017`)==F) & ((data$`2017` - round(data$`2017`,1)) == 0) 
  & (is.na(data$`2016`)==F) & ((data$`2016` - round(data$`2016`,1)) == 0)
  & (is.na(data$`2015`)==F) & ((data$`2015` - round(data$`2015`,1)) == 0),]

# Remove North America as it is a continent
country <- country[-(country$`Country Name`=="North America"),] 


# Split data 2017 and 2009 ------------------------------------------------
# Data 2017
time_elec17 <- country[c('Country Name','Country Code',2017)] # 2017 data
time_elec17 <- time_elec17[complete.cases(time_elec17[,
                                        as.character(2017)]),] # Remove NA
names(time_elec17) <- c("COUNTRY_N", "COUNTRY_C", "X2017") # rename

# Data 2009
time_elec09 <- country[c('Country Name','Country Code',2009)] # 2009 data
time_elec09 <- time_elec09[complete.cases(time_elec09[,
                                        as.character(2009)]),] # Remove NA
names(time_elec09) <- c("COUNTRY_N", "COUNTRY_C", "X2009") # rename


# Histogram data 2017 -----------------------------------------------------
summary(time_elec17$X2017)
ggplot(time_elec17, aes(x=X2017))+
  geom_histogram(col="red", 
                 aes(fill=..count..)) +
  xlab("Time required in days")


# Kernel esimators data 2017 ----------------------------------------------
# Choose the optimal bandwidth
kernels <- eval(formals(h.mlcv.default)$kernel)

# Rule of thumb
h.rt <- dpik(time_elec17$X2017, scalest = "minim", kernel = "normal") 
# Least Squares Cross Validation
h.LS <- h.ucv(time_elec17$X2017, deriv.order = 0, kernel = kernels[1])$h
# Likelihood cross validation
h.ML <- h.mlcv(time_elec17$X2017, kernel = kernels[1])$h 

# Print results of the three methods
h.rt
h.ML
h.LS

# Extract (x,y) of binned kernel density estimate of density of the data
estRT <- bkde(time_elec17$X2017, bandwidth = h.rt)
estLS <- bkde(time_elec17$X2017, bandwidth=h.LS)
estML <- bkde(time_elec17$X2017, bandwidth=h.ML)

# Unlist and create a dataframe in order to use ggplot
est.LS <- data.frame(matrix(unlist(estLS), nrow = 401))
est.ML <- data.frame(matrix(unlist(estML), nrow = 401))
est.rt <- data.frame(matrix(unlist(estRT), nrow = 401))


# Plot data 2017 ----------------------------------------------------------
ggplot(time_elec17, aes(x=X2017, y=..density..))+
  geom_histogram(col="black",
                 fill = "orange",
                 alpha = 0.2) +
  geom_line(data = est.LS, aes(x = X1, y = X2, colour = "CVLS"),
            size = 0.8) +
  geom_line(data = est.ML, aes(x = X1, y = X2, colour = "CVML"),
            size = 0.8) +
  geom_line(data = est.rt, aes(x = X1, y = X2, colour = "RoT"), 
            size = 0.8, alpha = 0.6) +
  xlab("Time required in days")+
  xlim(0, 500)


# Histogram data 2009 -----------------------------------------------------
summary(time_elec09$X2009)
ggplot(time_elec09, aes(x=X2009))+
  geom_histogram(col="red", 
                 aes(fill=..count..)) +
  xlab("Time required in days")


# Kernel esimators data 2009 ----------------------------------------------
# Choose the optimal bandwidth
kernels <- eval(formals(h.mlcv.default)$kernel)

h.rt <- dpik(time_elec09$X2009, scalest = "minim", kernel = "normal")
h.LS <- h.ucv(time_elec09$X2009, deriv.order = 0, kernel = kernels[1])$h
h.ML <- h.mlcv(time_elec09$X2009, kernel = kernels[1])$h 

# Results of the three methods
h.rt
h.ML
h.LS

# Extract (x,y) of binned kernel density estimate of density of the data
estRT <- bkde(time_elec09$X2009, bandwidth = h.rt)
estLS <- bkde(time_elec09$X2009, bandwidth=h.LS)
estML <- bkde(time_elec09$X2009, bandwidth=h.ML)

# Unlist and create a dataframe in order to use ggplot
est.LS <- data.frame(matrix(unlist(estLS), nrow = 401))
est.ML <- data.frame(matrix(unlist(estML), nrow = 401))
est.rt <- data.frame(matrix(unlist(estRT), nrow = 401))


# Plot data 2009 ----------------------------------------------------------
ggplot(time_elec09, aes(x=X2009, y=..density..))+
  geom_histogram(col="black",
                 fill = "orange",
                 alpha = 0.2) +
  geom_line(data = est.LS, aes(x = X1, y = X2, colour = "CVLS"), 
            size = 0.8) +
  geom_line(data = est.ML, aes(x = X1, y = X2, colour = "CVML"), 
            size = 0.8) +
  geom_line(data = est.rt, aes(x = X1, y = X2, colour = "RoT"),
            size = 0.8, alpha = 0.6) +
  xlab("Time required in days")+
  xlim(0, 500)

# Compare variance
var(time_elec09$X2009)
var(time_elec17$X2017)


# Adaptive bandwidth ------------------------------------------------------
adap_bw <- npudensbw(time_elec17$`X2017`,bwtype='adaptive_nn',
                     bwmethod='cv.ls',ckertype="gaussian") # bandwidth
adap <- npudens(bws=adap_bw) # estimator
plot(adap)
print(adap$bw) # k in k-th nearest neighbor

# replace with following options:
# bwmethod: cv.ml, cv.ls, normal-reference
# ckertype: gaussian, epanechnikov, uniform
# bwtype: fixed, generalized_nn, adaptive_nn

hist(time_elec17$`X2017`,probability = T,main="Time to obtain electricity",
     xlab="Time required in days",ylim=c(0,0.015),breaks = 20) # Histogram

# Set different values of k
adap_bw$bw <- 5 # Set k = 5
adap_5 <- npudens(bws=adap_bw)
plot(adap_5,ylim=c(0,0.015),xlab = NULL, ylab=NULL,col="green",lwd=2)

adap_bw$bw <- 15 # Set k = 15
adap_15 <- npudens(bws=adap_bw)
plot(adap_15,ylim=c(0,0.015),xlab = NULL, ylab=NULL,col="red",lwd=2)

adap_bw$bw <- 50 # Set k = 50
adap_50 <- npudens(bws=adap_bw)
plot(adap_50,ylim=c(0,0.015),xlab = NULL, ylab=NULL,col="blue",lwd=2)


# Compare fixed and adaptive ----------------------------------------------
fix_bw <- npudensbw(time_elec17$`X2017`,  bwtype='fixed',
                    bwmethod='cv.ls',ckertype="gaussian")
fix <- npudens(bws=fix_bw)
plot(fix,ylim=c(0,0.015),xlab = NULL, ylab=NULL,col="red",lwd=2)
print(fix$bw) # value bandwidth h

adap_bw <- npudensbw(time_elec17$`X2017`, bwtype='adaptive_nn',
                     bwmethod='cv.ls',ckertype="gaussian")
adap <- npudens(bws=adap_bw)
plot(adap,ylim=c(0,0.015),xlab = NULL, ylab=NULL,col="blue",lwd=2)
print(adap$bw) # k in k-th nearest neighbor


# Compare adaptive bandwith between data 2017 and data 2009 ---------------
adap_bw <- npudensbw(time_elec09$`X2009`,bwtype='adaptive_nn',
                     bwmethod='cv.ml',ckertype="gaussian")
adap <- npudens(bws=adap_bw)
print(adap$bw) # k in k-th nearest neighbor
# replace different options to get different estimators for data 2009