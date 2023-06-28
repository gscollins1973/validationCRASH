rm(list = ls())
setwd('/Users/collinsg/OneDrive - Nexus365/CSM/research/Data/CRASH data/CRASH-2//')
#setwd('/Users/gary/OneDrive - Nexus365/CSM/research/Data/CRASH data/CRASH-2//')

suppressPackageStartupMessages({
  suppressWarnings({
    library(tidyverse)
    library(progress)
    library(rms)
    library(logistf)
  })
})

# function to calculate the c-statistic
cstat <- function(prob, y) {
  n1 <- sum(!y)
  n2 <- sum(y)
  U  <- sum(rank(prob)[!y]) - n1 * (n1 + 1) / 2
  return(1 - U / n1 / n2)
}

## read in the CRASH-2 data, subset it, and add some noise N(0,1) data
CRASH <- read.csv('CRASH-2 data_0.csv')
CRASH$death <- rep(0, nrow(CRASH))
CRASH$death[CRASH$icause > 0] <- 1
CRASH$iinjurytype <- factor(CRASH$iinjurytype, levels = c(1, 2, 3))
CRASH <- CRASH[, c("death", "iage", "isex", "iinjurytype", "isbp", "irr", "icc", "ihr", "igcs")]
CRASH$noise1  <- rnorm (nrow(CRASH), 0, 1)
CRASH$noise2  <- rnorm (nrow(CRASH), 0, 1)
CRASH$noise3  <- rnorm (nrow(CRASH), 0, 1)
CRASH$noise4  <- rnorm (nrow(CRASH), 0, 1)
CRASH$noise5  <- rnorm (nrow(CRASH), 0, 1)
CRASH$noise6  <- rnorm (nrow(CRASH), 0, 1)
CRASH$noise7  <- rnorm (nrow(CRASH), 0, 1)
CRASH$noise8  <- rnorm (nrow(CRASH), 0, 1)
CRASH$noise9  <- rnorm (nrow(CRASH), 0, 1)
CRASH$noise10 <- rnorm (nrow(CRASH), 0, 1)


# Impute missing values
CRASH.MI <- mice::mice(CRASH)

# Pull out the first imputed data set (not interested in pooling over the imputed data sets)
CRASH2 <- complete(CRASH.MI, 1)

# fit full model (with all noise predictors) to the entire CRASH2 data
fit.all <- lrm(death ~ iage + isex + isbp + igcs +
                       noise1 + noise2 + noise3 + noise4 + noise5 + 
                       noise6 + noise7 + noise8 + noise9 + noise10, data = CRASH2)

# calculate the minimum sample size needed to develop the model
set.seed(876234)
rr.ss <- pmsampsize::pmsampsize(type = 'b', 
                                prevalence = sum(CRASH2$death == 1)/nrow(CRASH2),
                                parameters = length(coef(fit.all))-1,
                                cstatistic = as.numeric(fit.all$stats['C']))

# Sample sizes for model development
N.DEV <- c(200, 300, 400, 500, 1000, 5000, 10000)

# separate out those with an without the event (used in the simulation)
CRASH2.0 <- CRASH2[CRASH2$death == 0, ]
CRASH2.1 <- CRASH2[CRASH2$death == 1, ]
N.1     <- round(prop.table(table(CRASH2$death))[2] * N.DEV, 0)
N.0     <- N.DEV - N.1

N.SIM <- 500
perf.stats       <- matrix(ncol = length(N.DEV), nrow = N.SIM)
optimism         <- matrix(ncol = length(N.DEV), nrow = N.SIM)
perf.stats.boot  <- matrix(ncol = length(N.DEV), nrow = N.SIM)
perf.stats.cross <- matrix(ncol = length(N.DEV), nrow = N.SIM)
perf.stats.D     <- matrix(ncol = length(N.DEV), nrow = N.SIM)
perf.stats.V     <- matrix(ncol = length(N.DEV), nrow = N.SIM)

# set up a progress bar to monitor the simulations
pb <- progress_bar$new(
  format = " simulation done [:bar] :percent eta: :eta",
  total = length(N.DEV) * N.SIM, clear = FALSE, width = 60)

# start of simulation
for(i in 1:length(N.DEV)){
  for(j in 1:N.SIM){
    pb$tick()
    
    index.0 <- sample(1:nrow(CRASH2.0), N.0[i])
    index.1 <- sample(1:nrow(CRASH2.1), N.1[i])
    
    CRASH2.new <- rbind(CRASH2.0[index.0, ], CRASH2.1[index.1, ])
    fit        <- lrm(death ~ iage + isex + isbp + igcs +
                        noise1 + noise2 + noise3 + noise4 + noise5 + 
                        noise6 + noise7 + noise8 + noise9 + noise10, data = CRASH2.new, x = T, y = T)
   
    if(!fit$fail){
      perf.stats[j, i] <- as.numeric(fit$stats['C'])
      boot.fit   <- rms::validate(fit, B = 100)
      perf.stats.boot[j, i] <- (1 + boot.fit[1, 5])/2
      optimism[j,i] <- boot.fit[1,4]/2
      
    } else {
      perf.stats[j, i]      <- NA
      perf.stats.boot[j, i] <- NA
      optimism[j, i]        <- NA
    }
    
    ### split sample 
    index.D <- sample(1:nrow(CRASH2), N.DEV[i]*0.7, replace = F)
    CRASH2.D  <- CRASH2[index.D, ]
    index.V <- sample((1:nrow(CRASH2))[!(1:nrow(CRASH2) %in% index.D)], N.DEV[i]-length(index.D), replace = F)
    CRASH2.V  <- CRASH2[index.V, ]
    
    ## all data to develop the model and apparent performance
    ## use firth's correction (to handle sparsity in Killip class)
    fit                <- lrm(death ~ iage + isex + isbp + igcs +
                                noise1 + noise2 + noise3 + noise4 + noise5 + 
                                noise6 + noise7 + noise8 + noise9 + noise10, data = CRASH2.D)
    if(!fit$fail){
      pred.D             <- predict(fit, newdata = CRASH2.D, type = 'fitted')
      perf.stats.D[j, i] <- cstat(prob = pred.D, y = CRASH2.D$death)
      pred.V             <- predict(fit, newdata = CRASH2.V, type = 'fitted')
      perf.stats.V[j, i] <- cstat(prob = pred.V, y = CRASH2.V$death)
    } else {
      perf.stats.D[j, i] <- NA
      perf.stats.V[j, i] <- NA
    }
  }
}

### organise the results
### arrange performance results to plot
stats        <- reshape2::melt(perf.stats)
stats.boot   <- reshape2::melt(perf.stats.boot)
stats.D      <- reshape2::melt(perf.stats.D)
stats.V      <- reshape2::melt(perf.stats.V)
  
stats$approach      <- rep("All data (apparent)",            nrow(stats))
stats.boot$approach <- rep("Bootstrap correction",           nrow(stats.boot))
stats.D$approach    <- rep("Split sample (apparent, 70%)",   nrow(stats.D))
stats.V$approach    <- rep("Split sample (validation, 30%)", nrow(stats.V))
  
stats.OUT           <- bind_rows(stats, stats.boot, stats.D, stats.V)
colnames(stats.OUT) <- c("Sim", "N.DEV", "value", "approach")
stats.OUT$N.DEV     <- factor(stats.OUT$N.DEV, levels = seq(1:length(N.DEV)), labels = paste("N=", N.DEV, sep = ''))
stats.OUT$approach  <- factor(stats.OUT$approach, levels = c("All data (apparent)", 
                                                             "Bootstrap correction",
                                                             "Split sample (apparent, 70%)",
                                                             "Split sample (validation, 30%)"))
OUT_mean            <- stats.OUT %>% filter(N.DEV == 'N=10000') %>%  dplyr::summarize(mean_val = mean(value, na.rm = T))
stats.OUT$value[stats.OUT$value < 0.5] <- NA

OUT_mean$mean_val <- as.numeric(fit.all$stats['C'])

optimism.OUT           <- reshape2::melt(optimism)
colnames(optimism.OUT) <- c("Sim", "N.DEV", "value")
optimism.OUT$N.DEV     <- factor(optimism.OUT$N.DEV, levels = seq(1:length(N.DEV)), labels = paste("N=", N.DEV, sep = ''))

stats.OUT %>% group_by(approach, N.DEV) %>% 
  summarise(L = quantile(value, prob = 0.25, na.rm = T),
            M = quantile(value, prob = 0.50, na.rm = T), 
            U = quantile(value, prob = 0.75, na.rm = T)) %>% 
  print(n = 30)

#### Figure 2 in the paper ####
stats.OUT2 <- stats.OUT %>% filter(approach == "All data (apparent)")
p2 <- ggplot(stats.OUT2, aes(x = N.DEV, y = value - as.numeric(OUT_mean), group = approach)) +
  geom_jitter(alpha = 0.2) + 
  stat_summary(
    fun      = median, 
    geom     = "errorbar",
    aes(ymax = after_stat(y), ymin = after_stat(y)), 
    position = position_dodge(width = 0.8), 
    width    = 0.25,
    colour   = 'red') + 
  geom_hline(yintercept = 0, colour = 'black') + 
  #geom_hline(data = OUT_mean, aes(yintercept = mean_val), colour = 'blue') + 
  xlab("Size of available data") + 
  #ylab("Apparent c-statistic") + 
  ylab(expression(hat(c)-c[large])) + 
  theme_bw() +
  theme(legend.position = "bottom") + 
  theme(legend.title = element_blank()) + 
  theme(legend.text = element_text(size = 10)) +
  theme(axis.text = element_text(size = 12)) + 
  guides(colour = guide_legend(override.aes = list(alpha = 1)))
p2

#### Figure 3 in the paper ####
p3 <- ggplot(stats.OUT, aes(x = N.DEV, y = value - as.numeric(OUT_mean), group = approach, shape = approach)) +
  geom_jitter(alpha = 0.3, aes(color = approach, shape = approach), 
              position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.8)) + 
  scale_colour_grey(start = 0.7, end = 0.1)+
  stat_summary(
    fun      = median, 
    geom     = "errorbar",
    aes(ymax = after_stat(y), ymin = after_stat(y)), 
    position = position_dodge(width = 0.8), 
    width    = 0.25,
    colour   = 'red') + 
  geom_hline(yintercept = 0, colour = 'black') + 
  xlab("Size of available data") + 
  ylab(expression(hat(c)-c[large])) + 
  theme_bw() +
  theme(legend.position = "bottom") + 
  theme(legend.title = element_blank()) + 
  theme(legend.text = element_text(size = 10)) +
  theme(axis.text = element_text(size = 6)) + 
  guides(colour = guide_legend(override.aes = list(alpha = 1), title = ""),
         shape = guide_legend(nrow = 1,
                              title = ""))
p3
### summaries
stats.OUT %>%
  group_by(approach, N.DEV) %>% 
  dplyr::summarize(mean_val = mean(value, na.rm = T), 
                   sd  = sd(value, na.rm = T), 
                   min = min(value, na.rm = T),
                   max = max(value, na.rm = T),
                   mad = mad(value, na.rm = T), 
                   LQ  = quantile(value, na.rm=T, prob = 0.025), 
                   UQ  = quantile(value, na.rm=T, prob = 0.975)) %>% 
  print(n = 30)


mean.bias <- function(X, g){
  mean((X - g) / g, na.rm = T) * 100
}

stats.OUT %>% group_by(approach, N.DEV) %>% 
  dplyr::summarize(bias = mean.bias(value, g = 0.815)) %>% print(n=30)

############################################
### plotting optimism (not in the paper) ###
############################################
p4 <- ggplot(optimism.OUT, aes(x = N.DEV, y = value)) +
  geom_jitter(alpha = 0.2) + 
  stat_summary(
    fun      = median, 
    geom     = "errorbar",
    aes(ymax = after_stat(y), ymin = after_stat(y)), 
    position = position_dodge(width = 0.8), 
    width    = 0.25,
    colour   = 'red') + 
  xlab("Size of available data") + 
  ylab("Optimism") + 
  theme_bw() +
  theme(legend.position = "bottom") + 
  theme(legend.title = element_blank()) + 
  theme(legend.text = element_text(size = 10)) +
  theme(axis.text = element_text(size = 6)) + 
  guides(colour = guide_legend(override.aes = list(alpha = 1)))
p4