#install.packages("quantmod")
library(quantmod)
#install.packages("zoo")
library(zoo)
#install.packages("xlsx")
library(xlsx)
#install.packages("openxlsx")
library(openxlsx)


rm(list=ls())

# Get data

Symbols <- c("ZURN.SW","UHR.SW","SGSN.SW","ROG.SW","CFR.SW","NOVN.SW",
          "NESN.SW","LONN.SW","LHN.SW","GIVN.SW","GEBN.SW","CSGN.SW","^SSMI")

getSymbols(Symbols) 

# Assign to dataframe "zoo"
# Get adjusted prices (6th column)
# Merging all assets together to ensure identical timeframe
prices.data <- merge.zoo(ZURN.SW[,6],UHR.SW[,6],SGSN.SW[,6],ROG.SW[,6],CFR.SW[,6],NOVN.SW[,6],NESN.SW[,6],LONN.SW[,6],LHN.SW[,6],GIVN.SW[,6],GEBN.SW[,6],CSGN.SW[,6],SSMI[,6])
prices.data <-na.omit(prices.data)
colnames(prices.data) <- Symbols
View(prices.data)

SSMI=prices.data[,-(1:12)]
prices.data= prices.data[,-13]

Dates=index(prices.data)
#Saving dates since they are lost in the dataframe conversion.
write.xlsx(as.data.frame(prices.data),"Prices_SwissPortfolio.xlsx")
write.xlsx(as.data.frame(Dates),"Dates.xlsx")
write.xlsx(as.data.frame(SSMI),"Benchmark.xlsx")