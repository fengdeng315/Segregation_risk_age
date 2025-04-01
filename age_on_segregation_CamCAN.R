library(readxl)
library(dplyr)
library(corrplot)
library(ggplot2)
library(tidyr)
library(effectsize)
library(broom)

Demographic_rest <- read_excel("Demographic1.xlsx", sheet = "Rest")
Demographic_rest<-rename(Demographic_rest,ICV=tiv_cubicmm)

# load cognitive data
Cog <- read.csv("Main_SyS_Measures_Sch_500_17Yeo.csv",header = TRUE, sep = ",")
Cog$Sex<-factor(Cog$Sex,levels=c("1","2"),labels=c("Male","Female"))

# grey matter volume
load("I:/Segregation/means_per_region_volume.Rdata")
TotalGMV <- rowSums(mean_per_region_volume[,2:35],na.rm=TRUE)
TotalGMV <- data.frame(mean_per_region_volume[,1],TotalGMV,mean_per_region_volume[,36])
colnames(TotalGMV) <- c("CCID","TotalGMV","Age")

# load imaging data
Rest_data <- read_excel("Rest_Power214_highpass_Positive_Sparsity_weighted.xlsx", sheet = "Pc")

# combined those data
Rest_All<-left_join(Demographic_rest,Cog)
Rest_All<-left_join(Rest_All,TotalGMV)
Rest_All<-left_join(Rest_All,Rest_data)

# code some data
Rest_All[Rest_All$Sex == 'Male','SexC'] = -1
Rest_All[Rest_All$Sex == 'Female','SexC'] = 1

Rest_All$ageQuad <- poly(Rest_All$Age,2) #1st linear, 2nd quad
Rest_All$agepoly_1 = datawizard::standardise(Rest_All$ageQuad[,1])
Rest_All$agepoly_2 = datawizard::standardise(Rest_All$ageQuad[,2])

# plot
Rest_All$resid_Pc <- resid(lm(Pc ~ SexC + MeanFD, data = Rest_All))
ggplot(Rest_All, aes(x = Age, y = resid_Pc)) + 
  geom_point(size = 3, alpha = 0.6, col = "black", fill = "grey80") + 
  stat_smooth(method = "lm", se = TRUE, fill = "grey60", formula = y ~ poly(x,2, raw = TRUE), linewidth = 3) + 
  scale_x_continuous(breaks = round(seq(20, max(80), by = 20),1),limits = c(15,90)) + 
  theme_bw() + 
  theme(panel.border = element_blank(), legend.position = "none", text = element_text(size=18),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black",linewidth = 1.5),
        axis.ticks = element_line(colour = "black", linewidth = 1.5)) + 
  labs(x='Age', y='Functional Segregation') + 
  theme(plot.title = element_text(hjust = 0.5))

# age effect on Pc
lm_model_1 <- lm(scale(Pc) ~ scale(ageQuad) + SexC + scale(MeanFD), data = Rest_All)
etasq(lm_model_1,anova=TRUE,partial = TRUE)
summary(lm_model_1) 
round(confint(lm_model_1), 2)

lm_model_2 <- lm(scale(Pc) ~ scale(ageQuad) + SexC + scale(TotalGMV) + scale(ICV) + scale(MeanFD), data = Rest_All)
etasq(lm_model_2,anova=TRUE,partial = TRUE)
summary(lm_model_2)
round(confint(lm_model_2), 2)

# Pc~Cognition relationship
lm_model_FI <- lm(scale(CattellPCA) ~ scale(Pc) + SexC + scale(MeanFD), data = Rest_All)
summary(lm_model_FI)
round(confint(lm_model_FI), 2)

lm_model_FI_2 <- lm(scale(CattellPCA) ~ scale(Pc) + scale(ageQuad) + SexC + scale(MeanFD), data = Rest_All)
summary(lm_model_FI_2)
round(confint(lm_model_FI_2), 2)

lm_model_Me <- lm(scale(Memory_700) ~ scale(Pc) + SexC + scale(MeanFD), data = Rest_All)
summary(lm_model_Me)
round(confint(lm_model_Me), 2)

lm_model_Me_2 <- lm(scale(Memory_700) ~ scale(Pc) + scale(ageQuad) + SexC + scale(MeanFD), data = Rest_All)
summary(lm_model_Me_2)
round(confint(lm_model_Me_2), 2)

# age impact on 10 different networks; agepoly_1 + agepoly_2
lm_Pc1<-lm(scale(Pc_1) ~ agepoly_1 + agepoly_2 + SexC + scale(MeanFD), data=Rest_All)
lm_Pc2<-lm(scale(Pc_2) ~ agepoly_1 + agepoly_2 + SexC + scale(MeanFD), data=Rest_All)
lm_Pc3<-lm(scale(Pc_3) ~ agepoly_1 + agepoly_2 + SexC + scale(MeanFD), data=Rest_All)
lm_Pc4<-lm(scale(Pc_4) ~ agepoly_1 + agepoly_2 + SexC + scale(MeanFD), data=Rest_All)
lm_Pc5<-lm(scale(Pc_5) ~ agepoly_1 + agepoly_2 + SexC + scale(MeanFD), data=Rest_All)
lm_Pc6<-lm(scale(Pc_6) ~ agepoly_1 + agepoly_2 + SexC + scale(MeanFD), data=Rest_All)
lm_Pc7<-lm(scale(Pc_7) ~ agepoly_1 + agepoly_2 + SexC + scale(MeanFD), data=Rest_All)
lm_Pc8<-lm(scale(Pc_8) ~ agepoly_1 + agepoly_2 + SexC + scale(MeanFD), data=Rest_All)
lm_Pc9<-lm(scale(Pc_9) ~ agepoly_1 + agepoly_2 + SexC + scale(MeanFD), data=Rest_All)
lm_Pc10<-lm(scale(Pc_10) ~ agepoly_1 + agepoly_2 + SexC + scale(MeanFD), data=Rest_All)

models <- list(lm_Pc1, lm_Pc2, lm_Pc3, lm_Pc4, lm_Pc5, 
               lm_Pc6, lm_Pc7, lm_Pc8, lm_Pc9, lm_Pc10)

results <- lapply(models, function(model) {
  tidy(model, conf.int = TRUE) %>%
    filter(term %in% c("agepoly_1", "agepoly_2")) %>%
    select(term, estimate, conf.low, conf.high)
})

results_df <- do.call(rbind, results) %>%
  mutate(Network = rep(c("SM1", "SM2", "VIS", "AUD", "DMN", "FPN", "VAN", "CON", "DAN", "SN"), each = 2)) %>%
  mutate(term = ifelse(term == "agepoly_1", "Linear Age Effect", "Quadratic Age Effect"))

results_df$Network<-factor(results_df$Network,levels = c(
  "SM1","CON","DAN","AUD","VAN","SM2","VIS","SN","DMN","FPN"
))

ggplot(results_df, aes(x = estimate, y = Network, xmin = conf.low, xmax = conf.high, col = term)) +
  geom_point(size = 3) +
  geom_errorbarh(height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_wrap(~term, scales = "free_x") +
  labs(x = "Effect of Age on Pc", y = "Network") +
  scale_color_manual(values = c("blue", "red")) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")

# Pc~Cognition relationship on individual networks
lm_Pc1<-lm(scale(Memory_700) ~ scale(Pc_1) + SexC + scale(MeanFD), data=Rest_All) #CattellPCA, Memory_700, + agepoly_1 + agepoly_2 + scale(Education)
lm_Pc2<-lm(scale(Memory_700) ~ scale(Pc_2) + SexC + scale(MeanFD), data=Rest_All)
lm_Pc3<-lm(scale(Memory_700) ~ scale(Pc_3) + SexC + scale(MeanFD), data=Rest_All)
lm_Pc4<-lm(scale(Memory_700) ~ scale(Pc_4) + SexC + scale(MeanFD), data=Rest_All)
lm_Pc5<-lm(scale(Memory_700) ~ scale(Pc_5) + SexC + scale(MeanFD), data=Rest_All)
lm_Pc6<-lm(scale(Memory_700) ~ scale(Pc_6) + SexC + scale(MeanFD), data=Rest_All)
lm_Pc7<-lm(scale(Memory_700) ~ scale(Pc_7) + SexC + scale(MeanFD), data=Rest_All)
lm_Pc8<-lm(scale(Memory_700) ~ scale(Pc_8) + SexC + scale(MeanFD), data=Rest_All)
lm_Pc9<-lm(scale(Memory_700) ~ scale(Pc_9) + SexC + scale(MeanFD), data=Rest_All)
lm_Pc10<-lm(scale(Memory_700) ~ scale(Pc_10) + SexC + scale(MeanFD), data=Rest_All)

models <- list(lm_Pc1, lm_Pc2, lm_Pc3, lm_Pc4, lm_Pc5, 
               lm_Pc6, lm_Pc7, lm_Pc8, lm_Pc9, lm_Pc10)

model_results <- lapply(models, function(model) {
  tidy(model, conf.int = TRUE) %>% 
    filter(grepl("scale\\(Pc_", term))
}) %>% 
  bind_rows(.id = "model") %>% 
  mutate(
    Network = factor(model, levels = 1:10, labels = c("SM1", "SM2", "VIS", "AUD", "DMN", "FPN", "VAN", "CON", "DAN", "SN")),
    significant = ifelse(p.value < 0.005, "*", "")
  )

model_results$Network<-factor(model_results$Network,levels = c(
  "SM1","CON","DAN","AUD","VAN","SM2","VIS","SN","DMN","FPN"
))

rm(lm_Pc1,lm_Pc2,lm_Pc3,lm_Pc4,lm_Pc5,lm_Pc6,lm_Pc7,lm_Pc8,lm_Pc9,lm_Pc10,models,model_results)