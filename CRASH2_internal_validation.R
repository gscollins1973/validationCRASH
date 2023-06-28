### Code to create Supplementary Figure 2, Supplementary Table 2
### The CRASH-2 and CRASH-3 data used in this paper are freely 
### available at https://freebird.lshtm.ac.uk. 

set.seed(1234)

suppressPackageStartupMessages({
  suppressWarnings({
    library(tidyverse)
    library(progress)
    library(rms)
    library(logistf)
    library(readxl)
  })
})

cstat <- function(prob, y) {
  n1 <- sum(!y)
  n2 <- sum(y)
  U  <- sum(rank(prob)[!y]) - n1 * (n1 + 1) / 2
  return(1 - U / n1 / n2)
}

### DEVELOPMENT DATA
CRASH2 <- read_xlsx("CRASH-2 data_0.xlsx")
CRASH2$sex <- rep(0, nrow(CRASH2))
CRASH2$sex[CRASH2$isex == 1] <- 1

### VALIDATION DATA
CRASH3.orig <- read_xlsx("CRASH-3_dataset_anonymised_for_Freebird.xlsx")
CRASH3.orig$sex2 <- rep(0, nrow(CRASH3.orig))
CRASH3.orig$sex2[CRASH3.orig$sex == 'Male'] <- 1

### impute missing data in CRASH3 and take the first imputed dataset
CRASH3.orig <- CRASH3.orig[, c("causeDeath", "age", "sex2", "systolicBloodPressure", "gcsEyeOpening", "gcsMotorResponse", "gcsVerbalResponse")]
CRASH3.MI   <- mice::mice(CRASH3.orig)
CRASH3      <- mice::complete(CRASH3.MI, 1)

### create outcome in CRASH2
CRASH2$death <- rep(0, nrow(CRASH2))
CRASH2$death[CRASH2$icause > 0] <- 1

### create outcome in CRASH3
CRASH3$death <- rep(0, nrow(CRASH3))
CRASH3$death[CRASH3$causeDeath > 0] <- 1

### pull out predictors in CRASH2, impute and take the first imputed dataset
CRASH2a <- CRASH2[, c("death", "iage", "sex", "isbp", 'igcs')]
x <- mice::mice(CRASH2a)
CRASH2a <- complete(x, 1)

### pull out predictors in CRASH2
CRASH3a <- CRASH3[, c("death", "age", "sex2", "systolicBloodPressure")]
### manually calculate GCS in CRASH3 from component scores
c3.eye    <- as.numeric(vapply(strsplit(CRASH3$gcsEyeOpening,     " ", fixed=T), "[", "", 1))
c3.motor  <- as.numeric(vapply(strsplit(CRASH3$gcsMotorResponse,  " ", fixed=T), "[", "", 1))
c3.verbal <- as.numeric(vapply(strsplit(CRASH3$gcsVerbalResponse, " ", fixed=T), "[", "", 1))
x3        <- rowSums(cbind(c3.eye, c3.motor, c3.verbal), na.rm=T)
CRASH3a$gcs <- x3

### common labels
names(CRASH2a) <- c("death", "age", "sex", "sbp", 'gcs')
names(CRASH3a) <- c("death", "age", "sex", "sbp", "gcs")

### add noise variables
CRASH2a$noise1  <- rnorm(nrow(CRASH2), 0, 1)
CRASH2a$noise2  <- rnorm(nrow(CRASH2), 0, 1)
CRASH2a$noise3  <- rnorm(nrow(CRASH2), 0, 1)
CRASH2a$noise4  <- rnorm(nrow(CRASH2), 0, 1)
CRASH2a$noise5  <- rnorm(nrow(CRASH2), 0, 1)
CRASH2a$noise6  <- rnorm(nrow(CRASH2), 0, 1)
CRASH2a$noise7  <- rnorm(nrow(CRASH2), 0, 1)
CRASH2a$noise8  <- rnorm(nrow(CRASH2), 0, 1)
CRASH2a$noise9  <- rnorm(nrow(CRASH2), 0, 1)
CRASH2a$noise10 <- rnorm(nrow(CRASH2), 0, 1)

CRASH3a$noise1  <- rnorm(nrow(CRASH3), 0, 1)
CRASH3a$noise2  <- rnorm(nrow(CRASH3), 0, 1)
CRASH3a$noise3  <- rnorm(nrow(CRASH3), 0, 1)
CRASH3a$noise4  <- rnorm(nrow(CRASH3), 0, 1)
CRASH3a$noise5  <- rnorm(nrow(CRASH3), 0, 1)
CRASH3a$noise6  <- rnorm(nrow(CRASH3), 0, 1)
CRASH3a$noise7  <- rnorm(nrow(CRASH3), 0, 1)
CRASH3a$noise8  <- rnorm(nrow(CRASH3), 0, 1)
CRASH3a$noise9  <- rnorm(nrow(CRASH3), 0, 1)
CRASH3a$noise10 <- rnorm(nrow(CRASH3), 0, 1)

CRASH2 <- CRASH2a
CRASH3 <- CRASH3a

#### FIT MODEL TO ALL DATA
fit.all         <- lrm(death ~ age + sex + sbp + gcs + 
                         noise1 + noise2 + noise3 + noise4 + noise5 +
                         noise6 + noise7 + noise8 + noise9 + noise10, data = CRASH2)
pred.CRASH3.all <- predict(fit.all, newdata = CRASH3, type = 'fitted')
CalibrationCurves::val.prob.ci.2(p = pred.CRASH3.all, y = CRASH3$death)
cstat(prob = pred.CRASH3.all, y = CRASH3$death)


### Calculate the minimum sample size based on the c-statistic of the model fit to all data
rr.ss <- pmsampsize::pmsampsize(type = 'b', 
                                prevalence = sum(CRASH2$death == 1)/nrow(CRASH2),
                                parameters = length(coef(fit.all))-1,
                                cstatistic = as.numeric(fit.all$stats['C']))


### Sample sizes of the development cohort
#N.DEV <- sort(c(200, 300, 400, 500, 1000, 5000, 10000, ceiling(rr.ss$results_table[4,1]/0.7)))
N.DEV <- sort(c(200, 300, 400, 500, 1000, 5000, 10000))

CRASH2.0 <- CRASH2[CRASH2$death == 0, ]
CRASH2.1 <- CRASH2[CRASH2$death == 1, ]
N.1     <- round(prop.table(table(CRASH2$death))[2] * N.DEV, 0)
N.0     <- N.DEV - N.1

N.SIM <- 500
c.stats                 <- matrix(ncol = length(N.DEV), nrow = N.SIM)
c.stats.boot            <- matrix(ncol = length(N.DEV), nrow = N.SIM)
c.stats.CRASH3          <- matrix(ncol = length(N.DEV), nrow = N.SIM)
slope.stats.CRASH3      <- matrix(ncol = length(N.DEV), nrow = N.SIM)
slope.stats.boot        <- matrix(ncol = length(N.DEV), nrow = N.SIM)

c.stats.D.50            <- matrix(ncol = length(N.DEV), nrow = N.SIM)
c.stats.D.70            <- matrix(ncol = length(N.DEV), nrow = N.SIM)
c.stats.V.50            <- matrix(ncol = length(N.DEV), nrow = N.SIM)
c.stats.V.70            <- matrix(ncol = length(N.DEV), nrow = N.SIM)
c.stats.D.CRASH3.50     <- matrix(ncol = length(N.DEV), nrow = N.SIM)
c.stats.D.CRASH3.70     <- matrix(ncol = length(N.DEV), nrow = N.SIM)
slope.stats.V.50        <- matrix(ncol = length(N.DEV), nrow = N.SIM)
slope.stats.V.70        <- matrix(ncol = length(N.DEV), nrow = N.SIM)
slope.stats.D.CRASH3.50 <- matrix(ncol = length(N.DEV), nrow = N.SIM)
slope.stats.D.CRASH3.70 <- matrix(ncol = length(N.DEV), nrow = N.SIM)

pb <- progress_bar$new(
  format = " simulation done [:bar] :percent eta: :eta",
  total = length(N.DEV) * N.SIM, clear = FALSE, width = 60)

for(i in 1:length(N.DEV)){
  for(j in 1:N.SIM){
    pb$tick()
    
    ## randomly sample from the development cohort
    index.0    <- sample(1:nrow(CRASH2.0), N.0[i])
    index.1    <- sample(1:nrow(CRASH2.1), N.1[i])
    CRASH2.new <- rbind(CRASH2.0[index.0, ], CRASH2.1[index.1, ])
    
    ## fit model to the development data
    fit <- lrm(death ~ age + sex + sbp + gcs + 
                 noise1 + noise2 + noise3 + noise4 + noise5 +
                 noise6 + noise7 + noise8 + noise9 + noise10, data = CRASH2.new, x = T, y = T)
    
    if(!fit$fail){ # catch the model failures (lack of convergence)
      c.stats[j, i]            <- as.numeric(fit$stats['C'])    # apparent c-statistic
      boot.fit                 <- rms::validate(fit, B = 100)
      c.stats.boot[j, i]       <- (1 + boot.fit[1, 5])/2        # bootstrap corrected c-statistic
      slope.stats.boot[j, i]   <- boot.fit[4, 5]                # bootstrap corrected calibration slope
      
      ## External validation
      pred                     <- predict(fit, newdata = CRASH3, type = 'lp')
      c.stats.CRASH3[j, i]     <- cstat(prob = plogis(pred), y = CRASH3$death)
      slope.stats.CRASH3[j, i] <- as.numeric(coef(glm(CRASH3$death~pred, family = 'binomial'))[2])
      
    } else {
      c.stats[j, i]            <- NA
      c.stats.boot[j, i]       <- NA
      c.stats.CRASH3[j, i]     <- NA
      slope.stats.boot[j, i]   <- NA
      slope.stats.CRASH3[j, i] <- NA
    }
    
    ### split sample 50:50
    index <- sample(1:nrow(CRASH2.new), N.DEV[i] * 0.5, replace = F)
    CRASH2.D <- CRASH2.new[index,]
    CRASH2.V <- CRASH2.new[(1:nrow(CRASH2.new))[!(1:nrow(CRASH2.new) %in% index)], ]
    
    ## all data to develop the model and apparent performance
    fit <- lrm(death ~ age + sex + sbp + gcs + 
                 noise1 + noise2 + noise3 + noise4 + noise5 +
                 noise6 + noise7 + noise8 + noise9 + noise10, data = CRASH2.D)
    if(!fit$fail){
      pred.D.50          <- predict(fit, newdata = CRASH2.D, type = 'lp')
      c.stats.D.50[j, i] <- cstat(prob = plogis(pred.D.50), y = CRASH2.D$death)
      
      ### internal validation
      pred.V.50              <- predict(fit, newdata = CRASH2.V, type = 'lp')
      c.stats.V.50[j, i]     <- cstat(prob = plogis(pred.V.50), y = CRASH2.V$death)
      slope.stats.V.50[j, i] <- as.numeric(coef(glm(CRASH2.V$death~pred.V.50, family = 'binomial'))[2])
      
      ### external validation
      pred.50                       <- predict(fit, newdata = CRASH3, type = 'lp')
      c.stats.D.CRASH3.50[j, i]     <- cstat(prob = plogis(pred.50), y = CRASH3$death)
      slope.stats.D.CRASH3.50[j, i] <- as.numeric(coef(glm(CRASH3$death~pred.50, family = 'binomial'))[2])
    } else {
      c.stats.D.50[j, i]            <- NA
      c.stats.V.50[j, i]            <- NA
      c.stats.D.CRASH3.50[j, i]     <- NA
      slope.stats.V.50[j, i]        <- NA
      slope.stats.D.CRASH3.50[j, i] <- NA
    }
    
    ### split sample 70:30
    index <- sample(1:nrow(CRASH2.new), N.DEV[i] * 0.7, replace = F)
    CRASH2.D <- CRASH2.new[index,]
    CRASH2.V <- CRASH2.new[(1:nrow(CRASH2.new))[!(1:nrow(CRASH2.new) %in% index)], ]
    
    ## all data to develop the model and apparent performance
    fit <- lrm(death ~ age + sex + sbp + gcs + 
                 noise1 + noise2 + noise3 + noise4 + noise5 +
                 noise6 + noise7 + noise8 + noise9 + noise10, data = CRASH2.D)
    if(!fit$fail){
      pred.D.70          <- predict(fit, newdata = CRASH2.D, type = 'lp')
      c.stats.D.70[j, i] <- cstat(prob = plogis(pred.D.70), y = CRASH2.D$death)
      
      ### internal validation
      pred.V.70              <- predict(fit, newdata = CRASH2.V, type = 'lp')
      c.stats.V.70[j, i]     <- cstat(prob = plogis(pred.V.70), y = CRASH2.V$death)
      slope.stats.V.70[j, i] <- as.numeric(coef(glm(CRASH2.V$death~pred.V.70, family = 'binomial'))[2])
      
      ### external validation
      pred.70                       <- predict(fit, newdata = CRASH3, type = 'lp')
      c.stats.D.CRASH3.70[j, i]     <- cstat(prob = plogis(pred.70), y = CRASH3$death)
      slope.stats.D.CRASH3.70[j, i] <- as.numeric(coef(glm(CRASH3$death~pred.70, family = 'binomial'))[2])
    } else {
      c.stats.D.70[j, i]            <- NA
      c.stats.V.70[j, i]            <- NA
      c.stats.D.CRASH.703[j, i]     <- NA
      slope.stats.V.70[j, i]        <- NA
      slope.stats.D.CRASH3.70[j, i] <- NA
    }
    
  }
}


### arrange performance results to plot
stats             <- reshape2::melt(c.stats)
stats.boot        <- reshape2::melt(c.stats.boot)
stats.CRASH3      <- reshape2::melt(c.stats.CRASH3)
stats.D.50        <- reshape2::melt(c.stats.D.50)
stats.D.70        <- reshape2::melt(c.stats.D.70)
stats.V.50        <- reshape2::melt(c.stats.V.50)
stats.V.70        <- reshape2::melt(c.stats.V.70)
stats.D.CRASH3.50 <- reshape2::melt(c.stats.D.CRASH3.50)
stats.D.CRASH3.70 <- reshape2::melt(c.stats.D.CRASH3.70)
 
stats.slope.CRASH3      <- reshape2::melt(slope.stats.CRASH3)
stats.slope.boot        <- reshape2::melt(slope.stats.boot)
stats.slope.V.50        <- reshape2::melt(slope.stats.V.50)
stats.slope.V.70        <- reshape2::melt(slope.stats.V.70)
stats.slope.D.CRASH3.50 <- reshape2::melt(slope.stats.D.CRASH3.50)
stats.slope.D.CRASH3.70 <- reshape2::melt(slope.stats.D.CRASH3.70)

stats             <- add_column(stats,             approach = rep("All data (apparent)",            nrow(stats)))
stats.boot        <- add_column(stats.boot,        approach = rep("Bootstrap correction",           nrow(stats.boot)))
stats.CRASH3      <- add_column(stats.CRASH3,      approach = rep("All data (external)",            nrow(stats.CRASH3)))
stats.D.50        <- add_column(stats.D.50,        approach = rep("Split sample (apparent, 50%)",   nrow(stats.D.50)))
stats.D.70        <- add_column(stats.D.70,        approach = rep("Split sample (apparent, 70%)",   nrow(stats.D.70)))
stats.V.50        <- add_column(stats.V.50,        approach = rep("Split sample (validation, 50%)", nrow(stats.V.50)))
stats.V.70        <- add_column(stats.V.70,        approach = rep("Split sample (validation, 70%)", nrow(stats.V.70)))
stats.D.CRASH3.50 <- add_column(stats.D.CRASH3.50, approach = rep("Split sample (external, 50%)",   nrow(stats.D.CRASH3.50)))
stats.D.CRASH3.70 <- add_column(stats.D.CRASH3.70, approach = rep("Split sample (external, 70%)",   nrow(stats.D.CRASH3.70)))

stats.slope.CRASH3      <- add_column(stats.slope.CRASH3,      approach = rep("All data (external)",            nrow(stats.slope.CRASH3)))
stats.slope.boot        <- add_column(stats.slope.boot,        approach = rep("Bootstrap correction",           nrow(stats.slope.boot)))
stats.slope.V.50        <- add_column(stats.slope.V.50,        approach = rep("Split sample (validation, 50%)", nrow(stats.slope.V.50)))
stats.slope.V.70        <- add_column(stats.slope.V.70,        approach = rep("Split sample (validation, 70%)", nrow(stats.slope.V.70)))
stats.slope.D.CRASH3.50 <- add_column(stats.slope.D.CRASH3.50, approach = rep("Split sample (external, 50%)",   nrow(stats.slope.D.CRASH3.50)))
stats.slope.D.CRASH3.70 <- add_column(stats.slope.D.CRASH3.70, approach = rep("Split sample (external, 70%)",   nrow(stats.slope.D.CRASH3.70)))

stats             <- add_column(stats,             approach2 = rep("Apparent", nrow(stats)))
stats.boot        <- add_column(stats.boot,        approach2 = rep("Internal", nrow(stats.boot)))
stats.CRASH3      <- add_column(stats.CRASH3,      approach2 = rep("External", nrow(stats.CRASH3)))
stats.D.50        <- add_column(stats.D.50,        approach2 = rep("Apparent", nrow(stats.D.50)))
stats.D.70        <- add_column(stats.D.70,        approach2 = rep("Apparent", nrow(stats.D.70)))
stats.V.50        <- add_column(stats.V.50,        approach2 = rep("Internal", nrow(stats.V.50)))
stats.V.70        <- add_column(stats.V.70,        approach2 = rep("Internal", nrow(stats.V.70)))
stats.D.CRASH3.50 <- add_column(stats.D.CRASH3.50, approach2 = rep("External", nrow(stats.D.CRASH3.50)))
stats.D.CRASH3.70 <- add_column(stats.D.CRASH3.70, approach2 = rep("External", nrow(stats.D.CRASH3.70)))

stats.slope.CRASH3      <- add_column(stats.slope.CRASH3,      approach2 = rep("External", nrow(stats.slope.CRASH3)))
stats.slope.boot        <- add_column(stats.slope.boot,        approach2 = rep("Internal", nrow(stats.slope.boot)))
stats.slope.V.50        <- add_column(stats.slope.V.50,        approach2 = rep("Internal", nrow(stats.slope.V.50)))
stats.slope.V.70        <- add_column(stats.slope.V.70,        approach2 = rep("Internal", nrow(stats.slope.V.70)))
stats.slope.D.CRASH3.50 <- add_column(stats.slope.D.CRASH3.50, approach2 = rep("External", nrow(stats.slope.D.CRASH3.50)))
stats.slope.D.CRASH3.70 <- add_column(stats.slope.D.CRASH3.70, approach2 = rep("External", nrow(stats.slope.D.CRASH3.70)))

stats.OUT <- bind_rows(stats, stats.boot, stats.D.50, stats.D.70, stats.V.50, stats.V.70, stats.CRASH3, stats.D.CRASH3.50, stats.D.CRASH3.70)
slope.OUT <- bind_rows(stats.slope.CRASH3, stats.slope.boot, stats.slope.V.50, stats.slope.V.70, stats.slope.D.CRASH3.50, stats.slope.D.CRASH3.70)

stats.OUT <- add_column(stats.OUT, measure = rep("c-statistic",       nrow(stats.OUT)))
slope.OUT <- add_column(slope.OUT, measure = rep("Calibration slope", nrow(slope.OUT)))

stats.OUT <- rename(stats.OUT, Sim = Var1, N = Var2)
slope.OUT <- rename(slope.OUT, Sim = Var1, N = Var2)

stats.OUT <- mutate(stats.OUT, N = factor(N, levels = seq(1:length(N.DEV)), labels = paste("N=", N.DEV, sep = '')))
slope.OUT <- mutate(slope.OUT, N = factor(N, levels = seq(1:length(N.DEV)), labels = paste("N=", N.DEV, sep = '')))

stats.OUT <- mutate(stats.OUT, approach = factor(approach, levels = c("All data (apparent)", 
                                                                      "Bootstrap correction",
                                                                      "Split sample (apparent, 70%)",
                                                                      "Split sample (apparent, 50%)",
                                                                      "Split sample (validation, 70%)",
                                                                      "Split sample (validation, 50%)",
                                                                      "All data (external)",
                                                                      "Split sample (external, 70%)",
                                                                      "Split sample (external, 50%)")))

slope.OUT <- mutate(slope.OUT, approach = factor(approach, levels = c("Bootstrap correction",
                                                                      "Split sample (validation, 70%)",
                                                                      "Split sample (validation, 50%)",
                                                                      "All data (external)",
                                                                      "Split sample (external, 70%)",
                                                                      "Split sample (external, 50%)")))

stats.OUT <- mutate(stats.OUT, approach2 = factor(approach2, levels = c("Apparent", 
                                                                        "Internal", 
                                                                        "External")))

slope.OUT <- mutate(slope.OUT, approach2 = factor(approach2, levels = c("Internal", 
                                                                        "External")))

OUT <- bind_rows(stats.OUT, slope.OUT)
OUT <- mutate(OUT, measure = factor(measure, levels = c("c-statistic", "Calibration slope")))

### plot external validation (c-statistic, calibration slope)
### Supplementary Table 2
OUT8 <- OUT %>% filter(approach2 == "External")
OUT_mean <- OUT8 %>% filter(N == 'N=10000') %>% group_by(measure) %>% dplyr::summarize(mean_val = mean(value, na.rm = T))
OUT_mean[1,2] <- as.numeric(fit.all$stats['C'])
OUT_mean[2, 2] <- 1.0

p8 <- ggplot(OUT8, aes(x = N, y = value, group = approach, shape = approach)) +
  geom_jitter(alpha = 0.2, aes(color = approach, shape = approach), 
              position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.8)) + 
  facet_grid(measure~., scales = 'free') +
  scale_colour_grey(start = 0.1, end = 0.6)+
  stat_summary(
    fun      = median, 
    geom     = "errorbar",
    aes(ymax = after_stat(y), ymin = after_stat(y)), 
    position = position_dodge(width = 0.8), 
    width    = 0.25,
    colour   = 'red') + 
  geom_hline(data = OUT_mean, aes(yintercept = mean_val), colour = 'blue') + 
  xlab("Size of available data") + 
  ylab("") + 
  theme_bw() +
  theme(legend.position = "bottom") + 
  theme(legend.title = element_blank()) + 
  theme(legend.text = element_text(size = 10)) +
  theme(axis.text = element_text(size = 8)) + 
  guides(colour = guide_legend(override.aes = list(alpha = 1), title = ""),
         shape = guide_legend(nrow = 1,
                              title = ""))
p8

### summaries (supplementary Table 2)
OUT8 %>% filter(measure == 'c-statistic') %>% 
  group_by(approach, N) %>% 
  dplyr::summarize(mean_val = mean(value, na.rm = T), 
                   sd  = sd(value, na.rm = T), 
                   min = min(value, na.rm = T),
                   max = max(value, na.rm = T),
                   mad = mad(value, na.rm = T), 
                   LQ  = quantile(value, na.rm=T, prob = 0.025), 
                   UQ  = quantile(value, na.rm=T, prob = 0.975)) %>% 
  print(n = 30)

