library("dplyr")
library("writexl")
library("zoo")

# anti_Name <- c("AMI", "AMP", "ATM", "AUG", "AXO", "AZM", "CAZ", "CCV", "CEP", "CEQ", "CHL", "CIP", "COT", "CTC",
#                "CTX", "FEP", "FIS", "FOX", "GEN", "IMI", "KAN", "NAL", "PTZ", "SMX", "STR", "TET", "TIO")

##---------------------------------------------------ay##
## data cleaning
##---------------------------------------------------ay##
setwd("C:/Users/ayoung/Desktop/AMR-Linear")
suffixes <- c("Equiv", "Rslt", "Concl", "ConclPred")

common_columns <- c("Specimen.ID", "Data.Year", "Serotype", "Region.Name","Specimen.Source")

anti_Name <- c("CHL", "TIO")

data <- read.csv("./Ecoli.csv")
# group_SalT <- data %>% filter(Serotype == "I 4,[5],12:i:-") #CHL and TIO
group_SalT <- data %>% filter(Serotype == "O157:H7") #CHL and TIO
write.csv(group_SalT, "./DataBySeroAnti/SalT.csv")

data <- read.csv("./DataBySeroAnti/SalT.csv")
data <- data[-1]

for(antibiotic in anti_Name){
  
  selected_columns <- names(data)[grepl(paste0("^", antibiotic, "\\."), names(data))]
  
  
  selected_data <- data %>% select(all_of(c(common_columns, selected_columns)))
  
  
  selected_data <- selected_data %>% rename(Year = Data.Year)
  
  
  rename_columns <- function(cols, prefix) {
    gsub(paste0("^", prefix, "\\."), "", cols)
  }
  
  
  renamed_data <- selected_data %>%
    rename_with(~ rename_columns(., antibiotic), .cols = selected_columns)
  
  write.csv(renamed_data, paste0("DataBySeroAnti/SalT_", antibiotic, ".csv"), row.names = FALSE)
}

#need new column , type of censorship
# 0 - right
# 1- left
# 2- interval
##---------------------------------------------------ay##



# for (i in 1:2) {
#   # for (j in c(11, 27)) { #c(1:27)[-c(1, 6, 8, 10, 14, 16, 20, 26)] 
#   for (j in 1:2)
# print(paste0(c("SalT", "Sal4")[i], "/ ", anti_Name[j]))
# serotype <- c("SalT", "Sal4")[i]
# antibiotic <-  anti_Name[j]

serotype <- "SalT"
antibiotic <- "CHL"

min_year <- 2003
max_year <- 2023
##----------------------------------------------
## Loading data set 
##----------------------------------------------
# datFileName <- paste0("DataBySeroAnti/", serotype, "_", antibiotic, ".csv")

# dat <- read.csv(datFileName)
# # dat <-  dat[ , -1] # first column delete 
# dat <- dat %>% filter(Year >= 2002)
# dat <- as.data.frame(dat)

datFileName <- paste0("DataBySeroAnti/", serotype, "_", antibiotic, ".csv")
dat <- read.csv(datFileName, stringsAsFactors = F, header = T)
##----------------------------------------------
## remove * from Year
##----------------------------------------------
dat <- dat %>%
  mutate(Year = gsub("\\*", "", Year))
##----------------------------------------------

dat <- dat %>% filter(Year >= min_year) %>% filter(Year <= max_year)
##ay - pick only Stool
dat <- dat %>% filter(Specimen.Source == "Stool")

# # 
# dat <- dat %>%
#   group_by(Year) %>%
#   filter(n() > 20) %>%
#   ungroup()


dat <- dat[-1]
##---------------------------------------------------ay##
## order by Year
##---------------------------------------------------ay##
dat <- dat %>%
  arrange(Year)
##---------------------------------------------------ay##

##---------------------------------------------------ay##
## added l_vec, u_vec, censored
##---------------------------------------------------ay##

dat$l_vec <- NA
dat$u_vec <- NA


for (k in 1:nrow(dat)) {
  y_ij <- log2(dat$Rslt[k])
  
  if (dat$Equiv[k] == "<=" || dat$Equiv[k] == "<") {  # Left censoring
    dat$l_vec[k] <- 0.25
   
    dat$u_vec[k] <- y_ij
    
    dat$censored[k] <- 1 
    
  } else if (dat$Equiv[k] == ">") {  # Right censoring
    dat$l_vec[k] <- y_ij
    dat$u_vec[k] <- 32
    dat$censored[k] <- 0  
  }
  else {  # Interval censoring
    dat$l_vec[k] <- y_ij-1
    dat$u_vec[k] <-y_ij
    
    dat$censored[k] <- 2 
  }
  
}

##---------------------------------------------------ay##

y0 <- dat$l_vec
y1 <- dat$u_vec
censor <- dat$censored

##---------------------------------------------------ay##
## re-save csv file with l_vec, u_vec, censor data
##---------------------------------------------------ay##

write.csv(dat, paste0("DataBySeroAnti/SalT_",antibiotic, ".csv"))
##---------------------------------------------------ay##


minYear <- min(dat$Year)
yearLabel <- as.numeric(as.factor(as.numeric(dat$Year)))


uniqueYearLength <- length(unique(yearLabel))


##-------------------------------------------ay
## scatter plot by Year
##-------------------------------------------ay
# dat$Year <- as.numeric(as.character(dat$Year))
# 
# ggplot(dat, aes(x = Year)) +
#   geom_point(aes(y = Rslt, color = "Rslt"), shape = 16) +
#   labs(x = "Year", y = "Values", title = "Scatter Plot of Year") +
#   scale_color_manual(name = "Legend", values = c("Rslt" = "blue")) +
#   theme_minimal()

##------------------------------------------ay


##----------------------------------------------
## initial values : estimated from data set 
##----------------------------------------------
dat_group <- dat %>% mutate(Rslt_log2 = log(Rslt, base = 2),
                            cGroup = ifelse(Concl=="R", 1, 0)) %>% 
  group_by(Year, cGroup) %>% 
  # summarise(meanMIC = mean(Rslt_log2), 
  #           num = n(),
  #           meany0MIC = mean(l_vec, na.rm = TRUE)) %>%  
  summarise(meanMIC = mean(Rslt_log2), 
            num = n(),
            meany0MIC = mean(l_vec, na.rm = TRUE),
            .groups = "drop") %>% 
  data.frame() 

dat_totObs <- dat %>% group_by(Year) %>% summarise(totObs = n()) %>% data.frame() 
dat_full <- left_join(dat_group, dat_totObs) %>% mutate(prop = num / totObs)

## proportion of resistant 
p_tmp <-  left_join(data.frame(Year = sort(unique(dat$Year))), 
                    dat_full %>% filter(cGroup == 1) %>% dplyr::select(Year, prop))
if (any(is.na(p_tmp$prop))) {
  p_tmp$prop[is.na(p_tmp$prop)] <- 0
}
p <- p_tmp$prop
if (any(p < 1e-15)) {
  p[p < 1e-15] <- 0.000001
}
if (any(p > 1-(1e-15))) {
  p[p > 1-(1e-15)] <- 1-0.000001
}

## beta 
beta_tmp <- data.frame(Year = sort(unique(dat$Year)))

##----------------------------ay
## beta1 : pick not "R" data
## beta2 : pick "R" data
##----------------------------ay
beta_1 <- left_join(beta_tmp, dat_full %>% filter(cGroup == 0) %>% 
                      dplyr:: select(Year, meany0MIC) %>% rename(beta1 = meany0MIC)) 
beta_1_2 <- left_join(beta_1, dat_full %>% filter(cGroup == 1) %>% 
                        dplyr:: select(Year, meany0MIC) %>% rename(beta2 = meany0MIC))

if (any(is.na(beta_1_2$beta1))) {
  beta_1_2$beta1[which(is.na(beta_1_2$beta1))] <- mean(beta_1_2$beta1, na.rm = T)
}
if (any(is.na(beta_1_2$beta2))) {
  beta_1_2$beta2[which(is.na(beta_1_2$beta2))] <- mean(beta_1_2$beta2, na.rm = T)
}
beta <- beta_1_2[, c(2, 3)]


## initial value for allocation               
c <- ifelse(dat$Concl=="R", 1, 0)
## sigma
sigma <- dat %>% mutate(cGroup = ifelse(Concl=="R", 1, 0)) %>% 
  group_by(cGroup) %>% summarise(se = sd(l_vec, na.rm = TRUE)) %>% data.frame () %>% .[, "se"] 

SMALLVAL <- 0.0001

if (any(sigma == 0)) {
  sigma <- sigma + SMALLVAL 
}
if (any(is.na(sigma))) {
  sigma[is.na(sigma)] <- 0.1
}

## mu
lmod <- lm(beta$beta1 ~ c(1:nrow(beta)))

##---------------------------------------------------ay## 
## initialzie mu empty array
##---------------------------------------------------ay## 
mu <- data.frame(Intercept = numeric(nrow(beta)), Slope = numeric(nrow(beta)))
##---------------------------------------------------ay##


mu[[3]] <- lmod$coefficients
mu[[2]] <- mean(beta$beta2)
mu[[1]] <- lmod$coefficients[1] + lmod$coefficients[2] * c(1:nrow(beta))

## tau
tau <- c(summary(lmod)$sigma, apply(beta, 2, sd)[2])

## calculate logit(p)
logit <- function(x) {
  n <- length(x) 
  if (any(x < 1e-15)) {
    x[x < 1e-15] <- 1e-15
  }
  sapply(seq(1, n), function(s) { log(x[s] / (1 - x[s])) })
}

alpha <- logit(p)
muAlpha <- mean(alpha)
sigmaAlpha <- sd(alpha)

##--------------------------------------
source("Linear_mcmc.R")
if (!dir.exists(paste0("LinearModelOutput/From2002/", serotype, "_", antibiotic, "_Res2"))) {
  dir.create(paste0("LinearModelOutput/From2002/", serotype, "_", antibiotic, "_Res2"))
}

output <- paste0("LinearModelOutput/From2002/", serotype, "_", antibiotic, "_Res2/")
iterMax <- 200


model2_mcmc(y0, y1, c, p, beta, sigma, mu, tau, yearLabel, censor, 
            muAlpha, sigmaAlpha, alpha, ## initial values
            iterMax, output, prop.sd = 5, seed = 2019)    

source("Linear_res_extraction.R")  
extract(serotype, antibiotic,min_year,max_year)


