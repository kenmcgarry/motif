# AI2022_main.R
# Surprising patterns could be novel and interesting, worthy of investigation. Surprise is the difference 
# between our expectations and actual experiences. How surprised should we be based on the size of differences? 
# ...and if we are surprised by this new pattern is it valid? 

# commenced: 20/01/2021;
# submitted: 07/01/2022;
# 
library(xtable)
library(dplyr)
library(ggplot2)
library(reshape2)
library(infotheo)
library(philentropy)
library(stringr)
library(stringi)
library(stringdist)
library(pracma)
library(TraMineR)
library(PST)
library(SeqDetect)
library(tidyverse)
library(tidyr)
library(hrbrthemes)
library(xtable)
library(cluster)
library(seqHMM)
library(cluster)

setwd("C:/common_laptop/R-files/UKCI2022")
rm(list = ls()) # remove any legacy variables and data from previous analysis

source("AI2022_strfunctions.R")
source("AI2022_sequences.R")
source("AI2022_simple_example.R")




