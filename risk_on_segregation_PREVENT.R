library(readxl)
library(dplyr)
library(psych)
library(ggplot2)
library(corrplot)
library(Rmisc)
library(lme4)
library(lmerTest)
library(effectsize)
library(car)
library(rstatix)
library(tidyr)
library(gghalves)

## Load cognition data
Cog_t1 <- read_excel("Demo_risk_cog_t1_t2.xlsx", sheet = "baseline")
Cog_t2 <- read_excel("Demo_risk_cog_t1_t2.xlsx", sheet = "follow_up")

# recode data
Cog_t1$FH1 <- factor(Cog_t1$FH1, levels=c(0,1), labels = c("FH-","FH+"))
Cog_t2$FH2 <- factor(Cog_t2$FH2, levels=c(0,1), labels = c("FH-","FH+"))

Cog_t1$apoe4 <- factor(Cog_t1$apoe4, levels=c(0,1), labels = c("non-carriers","carriers"))
Cog_t2$apoe4 <- factor(Cog_t2$apoe4, levels=c(0,1), labels = c("non-carriers","carriers"))

Cog_t1$Sex <- factor(Cog_t1$Sex, levels=c(1,2), labels = c("Male","Female"))
Cog_t2$Sex <- factor(Cog_t2$Sex, levels=c(1,2), labels = c("Male","Female"))

# Rotated PCA to extract main cognitive domains
# parrallel analysis and scree plot to decide the number of components
ParaAna_t1<-fa.parallel(Cog_t1[,9:21], fa="pc", n.iter = 500)
ParaAna_t2<-fa.parallel(Cog_t2[,9:21], fa="pc", n.iter = 500)

# rotated PCA
PCs_t1<-principal(Cog_t1[,9:21], nfactors = 3, rotate="Varimax", eps=1e-7)
PCs_t2<-principal(Cog_t2[,9:21], nfactors = 3, rotate="Varimax", eps=1e-7)

# Visualise
fa.diagram(PCs_t1)
fa.diagram(PCs_t2)
corrplot(PCs_t1$loadings,cl.ratio = 0.6,number.cex=0.8,number.digits=2, addCoef.col="black")
corrplot(PCs_t2$loadings,cl.ratio = 0.6,number.cex=0.8,number.digits=2, addCoef.col="black")

# Factor similarity
fa.congruence(PCs_t1,PCs_t2)

Cog_t1 <- cbind(Cog_t1,PCs_t1$scores) 
Cog_t2 <- cbind(Cog_t2,PCs_t2$scores) %>% rename(RC3=RC2,RC2=RC3)

rm(ParaAna_t1,ParaAna_t2,PCs_t1,PCs_t2)

Cog_t1$RC3_inv<-Cog_t1$RC3*(-1)
Cog_t2$RC3_inv<-Cog_t2$RC3*(-1)

# load functional data
Seg_t1 <- read_excel("Segregation_Power214_noGSR.xlsx", sheet = "aucParameters_T1")
Seg_t2 <- read_excel("Segregation_Power214_noGSR.xlsx", sheet = "aucParameters_T2")

Seg_t1 <- inner_join(Seg_t1,Cog_t1,by="SubjID")
Seg_t2 <- inner_join(Seg_t2,Cog_t2,by="SubjID")

rm(Cog_t1,Cog_t2)

# exclude participants with nodes less than 
Seg_t1<-Seg_t1[!(Seg_t1$No_Nodes_retained<214*0.8),]
Seg_t2<-Seg_t2[!(Seg_t2$No_Nodes_retained<214*0.8),]

# exclude participants with headmovement larger than 
Seg_t1<-Seg_t1[!(Seg_t1$meanFD>0.4),]
Seg_t2<-Seg_t2[!(Seg_t2$meanFD>0.4),]

# load brain struct data
BrainStruc_t1 <- read_excel("BrainStruc_t1_t2.xlsx", sheet = "T1")
BrainStruc_t2 <- read_excel("BrainStruc_t1_t2.xlsx", sheet = "T2")

Seg_t1 <- inner_join(Seg_t1,BrainStruc_t1,by="SubjID")
Seg_t2 <- inner_join(Seg_t2,BrainStruc_t2,by="SubjID")

rm(BrainStruc_t1,BrainStruc_t2)
Seg_t1<-Seg_t1%>%rename(TICV=EstimatedTotalIntraCranialVol)
Seg_t2<-Seg_t2%>%rename(TICV=EstimatedTotalIntraCranialVol)

# load number of retained brain nodes per network
No_Nodes_Net_t1 <-read_excel("No_Nodes_Net.xlsx", sheet = "baseline")
No_Nodes_Net_t2 <-read_excel("No_Nodes_Net.xlsx", sheet = "followup")

Seg_t1 <- inner_join(Seg_t1,No_Nodes_Net_t1,by="SubjID")
Seg_t2 <- inner_join(Seg_t2,No_Nodes_Net_t2,by="SubjID")
rm(No_Nodes_Net_t1,No_Nodes_Net_t2)

# code some data
Seg_t1[Seg_t1$Sex == 'Male','SexC'] = -1
Seg_t1[Seg_t1$Sex == 'Female','SexC'] = 1

Seg_t1[Seg_t1$FH == 'FH-','FHC'] = -1
Seg_t1[Seg_t1$FH == 'FH+','FHC'] = 1

Seg_t1_clean <- Seg_t1 %>% filter(!is.na(apoe4))
Seg_t1_clean[Seg_t1_clean$apoe4 == 'non-carriers','apoe4C'] = -1
Seg_t1_clean[Seg_t1_clean$apoe4 == 'carriers','apoe4C'] = 1

# risk effect on global segregation 
lm_model_apoe4_1 <- lm(scale(Pc_Global) ~ apoe4C + scale(Age) + SexC + scale(meanFD) + scale(No_Nodes_retained), data = Seg_t1_clean)
summary(lm_model_apoe4_1)
round(confint(lm_model_apoe4_1), 2)

lm_model_apoe4_2 <- lm(scale(Pc_Global) ~ apoe4C + scale(Age) + SexC + scale(TotalGrayVol) + scale(TICV) + scale(meanFD) + scale(No_Nodes_retained), data = Seg_t1_clean) 
round(confint(lm_model_apoe4_2), 2)

lm_model_FH <- lm(scale(Pc_Global) ~ FHC + scale(Age) + SexC + scale(meanFD) + scale(No_Nodes_retained), data = Seg_t1)
summary(lm_model_FH)
round(confint(lm_model_FH), 2)

lm_model_CAIDE <- lm(scale(Pc_Global) ~ scale(CAIDE) + scale(meanFD) + scale(No_Nodes_retained), data = Seg_t1_clean)
summary(lm_model_CAIDE)
round(confint(lm_model_CAIDE), 2)

# risk effect on network segregation
lm_Pc5<-lm(scale(Pc_5) ~ apoe4C + scale(Age) + SexC + scale(meanFD) + scale(Network5), data=Seg_t1_clean) # DMN
lm_Pc7<-lm(scale(Pc_7) ~ apoe4C + scale(Age) + SexC + scale(meanFD) + scale(Network7), data=Seg_t1_clean) # FPN
lm_Pc8<-lm(scale(Pc_8) ~ apoe4C + scale(Age) + SexC + scale(meanFD) + scale(Network8), data=Seg_t1_clean) # SN

models <- list(lm_Pc5, lm_Pc7, lm_Pc8)

results <- lapply(models, function(model) {
  tidy(model, conf.int = TRUE) %>%
    filter(term == "apoe4C") %>%
    select(term, estimate, conf.low, conf.high)
})

results_df <- do.call(rbind, results) %>%
  mutate(Network = rep(c("DMN", "FPN", "SN"), each = 1)) %>% 
  mutate(term = "apoe4C")

results_df$Network<-factor(results_df$Network,levels = c("SN","DMN","FPN"))

# PLOT baseline results
Seg_t1_clean$resid_Pc <- resid(lm(Pc_Global ~ scale(Age) + SexC + scale(meanFD) + scale(No_Nodes_retained), data = Seg_t1_clean))
Seg_t1_clean$resid_Pc_5 <- resid(lm(Pc_5 ~ scale(Age) + SexC + scale(meanFD) + scale(Network5), data = Seg_t1_clean))

Seg_t1_clean$apoe4n<-as.numeric(Seg_t1_clean$apoe4) #noncarriers=1; carriers=2
set.seed(123)
Seg_t1_clean$apoe4nj<-jitter(Seg_t1_clean$apoe4n,amount=0.09)
GPc_apoe4_SE <- summarySE(Seg_t1_clean, measurevar="resid_Pc", groupvars="apoe4n")
Pc5_apoe4_SE <- summarySE(Seg_t1_clean, measurevar="resid_Pc_5", groupvars="apoe4n")

ggplot(data=filter(Seg_t1_clean, !is.na(apoe4n)), aes(y = resid_Pc))+ 
  geom_point(data=Seg_t1_clean %>% filter(apoe4n=="1"),aes(x = apoe4nj),color = "cornflowerblue", size=2,alpha=.4)+
  geom_point(data=Seg_t1_clean %>% filter(apoe4n=="2"),aes(x = apoe4nj),color = "salmon",size=2, alpha=.6)+
  geom_errorbar(data=filter(GPc_apoe4_SE,!is.na(apoe4n)), aes(x=apoe4n, ymin=resid_Pc-ci, ymax=resid_Pc+ci), width=.1) +
  geom_point(data=filter(GPc_apoe4_SE,!is.na(apoe4n)) %>% filter(apoe4n=="1"),aes(x = apoe4n,y=resid_Pc),shape=21, color="black", fill="cornflowerblue", size=2, stroke=1)+
  geom_point(data=filter(GPc_apoe4_SE,!is.na(apoe4n)) %>% filter(apoe4n=="2"),aes(x = apoe4n,y=resid_Pc),shape=21, color="black", fill="salmon", size=2, stroke=1)+
  geom_half_violin(data=filter(Seg_t1_clean, !is.na(apoe4n)) %>% filter(apoe4n=="1"), aes(x = apoe4n, y = resid_Pc), #CHANGE
                   position=position_nudge(x=-.2),side="l",width=.5, fill="cornflowerblue",alpha=.4)+
  geom_half_violin(data=filter(Seg_t1_clean, !is.na(apoe4n)) %>% filter(apoe4n=="2"), aes(x = apoe4n, y = resid_Pc), #CHANGE
                   position=position_nudge(x=.2),side="r",width=.5, fill="salmon",alpha=.6)+
  scale_x_continuous(breaks = c(1, 2),labels=c("non-carriers","carriers"))+ #limits = c(0,3)
  ylab("Global Pc")+ 
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01))+#,limits = c(-2,2),limits = c(-0.35,-0.25)
  theme(axis.text=element_text(size=12,colour = "black"),
        axis.title=element_text(size=12,colour = "black"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black",linewidth = 0.5),
        legend.position = "none",
        axis.title.x=element_blank())

rm(GPc_apoe4_SE,Pc5_apoe4_SE)

# segregation - cognition relationships at baseline 
lm_model_RC1_G <- lm(RC1 ~ scale(Pc_Global) + scale(Age) + SexC + scale(YeaoEdu) + scale(meanFD) + scale(No_Nodes_retained), data = Seg_t1)
lm_model_RC2_G <- lm(RC2 ~ scale(Pc_Global) + scale(Age) + SexC + scale(YeaoEdu) + scale(meanFD) + scale(No_Nodes_retained), data = Seg_t1)
lm_model_RC3_G <- lm(RC3_inv ~ scale(Pc_Global) + scale(Age) + SexC + scale(YeaoEdu) + scale(meanFD) + scale(No_Nodes_retained), data = Seg_t1)
summary(lm_model_RC1_G)
summary(lm_model_RC2_G)
summary(lm_model_RC3_G)

lm_model_RC1_DMN <- lm(RC1 ~ scale(Pc_5) + scale(Age) + SexC + scale(YeaoEdu) + scale(meanFD) + scale(Network5), data = Seg_t1)
lm_model_RC2_DMN <- lm(RC2 ~ scale(Pc_5) + scale(Age) + SexC + scale(YeaoEdu) + scale(meanFD) + scale(Network5), data = Seg_t1)
lm_model_RC3_DMN <- lm(RC3_inv ~ scale(Pc_5) + scale(Age) + SexC + scale(YeaoEdu) + scale(meanFD) + scale(Network5), data = Seg_t1)
summary(lm_model_RC1_DMN)
summary(lm_model_RC2_DMN)
summary(lm_model_RC3_DMN)

# plot Pc~Cog relationship
Seg_t1 = Seg_t1[complete.cases(Seg_t1$RC1), ]
Seg_t1$RC1_resids <- resid(lm(RC1 ~  scale(Age) + SexC + scale(YeaoEdu) + scale(meanFD) + scale(No_Nodes_retained), data = Seg_t1))
Seg_t1$Pc_Global_resids <- resid(lm(Pc_Global ~ scale(Age) + SexC + scale(YeaoEdu) + scale(meanFD) + scale(No_Nodes_retained), data = Seg_t1))

Seg_t1$RC1_resids2 <- resid(lm(RC1 ~ scale(Age) + SexC + scale(YeaoEdu) + scale(meanFD) + scale(Network5), data = Seg_t1))
Seg_t1$Pc5_resids <- resid(lm(Pc_5 ~ scale(Age) + SexC + scale(YeaoEdu) + scale(meanFD) + scale(Network5), data = Seg_t1))

ggplot(data=Seg_t1,#filter(Seg_t1, !is.na(apoe4.x)
       aes(x=Pc5_resids, y = RC1_resids2))+ #, color=apoe4.x
  geom_smooth(method=lm)+
  geom_jitter(size = 2, alpha=.6)+
  xlab("DMN Pc")+
  ylab("Episodic and relational memory")+ #∆ Response time in working memory
  scale_x_continuous(labels = scales::number_format(accuracy = 0.01))+#,limits=c(0.10,0.35)
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01))+ #,limits=c(-12,16)
  theme(axis.text=element_text(size=12,colour = "black"),
        axis.title=element_text(size=12,colour = "black"),
        panel.border = element_blank(),  
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "grey",linewidth = 1),
        legend.position = "right")


# risk impact on longitudinal change of segregation 
Seg_all<-inner_join(Seg_t1_clean,Seg_t2,by="SubjID")

Pc_Global<-gather(Seg_all[,c("SubjID","apoe4C","Age.x","SexC","Pc_Global.x","Pc_Global.y")], Time, Pc_Global, c("Pc_Global.x","Pc_Global.y"))
Pc_Global$Time<-factor(Pc_Global$Time)
levels(Pc_Global$Time)[levels(Pc_Global$Time)=="Pc_Global.x"] <- "Baseline"
levels(Pc_Global$Time)[levels(Pc_Global$Time)=="Pc_Global.y"] <- "Follow_up"

Pc_5<-gather(Seg_all[,c("SubjID","Pc_5.x","Pc_5.y")], Time, Pc_5, c("Pc_5.x","Pc_5.y"))
Pc_5$Time<-factor(Pc_5$Time)
levels(Pc_5$Time)[levels(Pc_5$Time)=="Pc_5.x"] <- "Baseline"
levels(Pc_5$Time)[levels(Pc_5$Time)=="Pc_5.y"] <- "Follow_up"

meanFD<-gather(Seg_all[,c("SubjID","meanFD.x","meanFD.y")], Time, meanFD, c("meanFD.x","meanFD.y"))
meanFD$Time<-factor(meanFD$Time)
levels(meanFD$Time)[levels(meanFD$Time)=="meanFD.x"] <- "Baseline"
levels(meanFD$Time)[levels(meanFD$Time)=="meanFD.y"] <- "Follow_up"

No_Nodes_retained<-gather(Seg_all[,c("SubjID","No_Nodes_retained.x","No_Nodes_retained.y")], Time, No_Nodes_retained, c("No_Nodes_retained.x","No_Nodes_retained.y"))
No_Nodes_retained$Time<-factor(No_Nodes_retained$Time)
levels(No_Nodes_retained$Time)[levels(No_Nodes_retained$Time)=="No_Nodes_retained.x"] <- "Baseline"
levels(No_Nodes_retained$Time)[levels(No_Nodes_retained$Time)=="No_Nodes_retained.y"] <- "Follow_up"

Network5<-gather(Seg_all[,c("SubjID","Network5.x","Network5.y")], Time, Network5, c("Network5.x","Network5.y"))
Network5$Time<-factor(Network5$Time)
levels(Network5$Time)[levels(Network5$Time)=="Network5.x"] <- "Baseline"
levels(Network5$Time)[levels(Network5$Time)=="Network5.y"] <- "Follow_up"

RC1<-gather(Seg_all[,c("SubjID","RC1.x","RC1.y")], Time, RC1, c("RC1.x","RC1.y"))
RC1$Time<-factor(RC1$Time)
levels(RC1$Time)[levels(RC1$Time)=="RC1.x"] <- "Baseline"
levels(RC1$Time)[levels(RC1$Time)=="RC1.y"] <- "Follow_up"

TotalGrayVol<-gather(Seg_all[,c("SubjID","TotalGrayVol.x","TotalGrayVol.y")], Time, TotalGrayVol, c("TotalGrayVol.x","TotalGrayVol.y"))
TotalGrayVol$Time<-factor(TotalGrayVol$Time)
levels(TotalGrayVol$Time)[levels(TotalGrayVol$Time)=="TotalGrayVol.x"] <- "Baseline"
levels(TotalGrayVol$Time)[levels(TotalGrayVol$Time)=="TotalGrayVol.y"] <- "Follow_up"

TICV<-gather(Seg_all[,c("SubjID","TICV.x","TICV.y")], Time, TICV, c("TICV.x","TICV.y"))
TICV$Time<-factor(TICV$Time)
levels(TICV$Time)[levels(TICV$Time)=="TICV.x"] <- "Baseline"
levels(TICV$Time)[levels(TICV$Time)=="TICV.y"] <- "Follow_up"

Pc<-left_join(Pc_Global,Pc_5)
Pc<-left_join(Pc,meanFD)
Pc<-left_join(Pc,No_Nodes_retained)
Pc<-left_join(Pc,Network5)
Pc<-left_join(Pc,RC1)
Pc<-left_join(Pc,TotalGrayVol)
Pc<-left_join(Pc,TICV)

Pc[Pc$Time == 'Baseline','TimeC'] = -1
Pc[Pc$Time == 'Follow_up','TimeC'] = 1

model <- lmer(scale(Pc_Global) ~ TimeC*apoe4C + scale(Age.x) + SexC + scale(meanFD) + scale(No_Nodes_retained) + (1 | SubjID), data = Pc)
model_5 <- lmer(scale(Pc_5) ~ TimeC*apoe4C+ scale(Age.x) + SexC + scale(meanFD) + scale(Network5) + (1 | SubjID), data = Pc)

summary(model)
round(confint(model), 2)
summary(model_5)
round(confint(model_5), 2)

model_cog <- lmer(RC1 ~ TimeC*scale(Pc_Global) + scale(Age.x) + SexC + scale(meanFD) + scale(No_Nodes_retained) + (1 | SubjID), data = Pc)
model_cog_5 <- lmer(RC1 ~ TimeC*scale(Pc_5)+ scale(Age.x) + SexC + scale(meanFD) + scale(Network5) + (1 | SubjID), data = Pc) 

summary(model_cog)
summary(model_cog_5)

# simple effects testing
noncarriers<-filter(Seg_all,apoe4.x=="non-carriers")
carriers<-filter(Seg_all,apoe4.x=="carriers")

t.test(carriers$Pc_5.x,carriers$Pc_5.y,paired = T)#p=0.007
t.test(noncarriers$Pc_5.x,noncarriers$Pc_5.y,paired = T)#p=0.49

# plot longitudinal results
Pc$resid_Pc <- resid(lmer(Pc_Global ~ scale(Age.x) + SexC + scale(meanFD) + scale(No_Nodes_retained) + (1 | SubjID), data = Pc))
Pc$resid_Pc_5 <- resid(lmer(Pc_5 ~ scale(Age.x) + SexC + scale(meanFD) + scale(Network5) + (1 | SubjID), data = Pc))

Pc$Timen<-as.numeric(Pc$Time)
set.seed(123)
Pc$Timenj<-jitter(Pc$Timen,amount=0.09)

noncarriers<-filter(Pc,apoe4C==-1)
carriers<-filter(Pc,apoe4C==1)

library(Rmisc)
noncarriers_SE <- summarySE(noncarriers, measurevar="resid_Pc", groupvars="Timen")
carriers_SE <- summarySE(carriers, measurevar="resid_Pc", groupvars="Timen")

ggplot(data=carriers, aes(y = resid_Pc))+ #CHANGE
  geom_point(data=carriers%>%filter(Timen=="1"), aes(x = Timenj), color = "salmon", size=2, alpha=.6)+
  geom_point(data=carriers%>%filter(Timen=="2"), aes(x = Timenj), color = "salmon", size=2, alpha=.6)+
  geom_line(data=carriers, aes(x=Timenj, group=SubjID), color='gray', alpha=.5)+
  geom_errorbar(data=carriers_SE, aes(x=Timen, ymin=resid_Pc-ci, ymax=resid_Pc+ci), width=.1) +
  geom_point(data=carriers_SE, aes(x=Timen, y=resid_Pc), shape=21, color="black", fill="salmon", size=2, stroke=1) +
  geom_line(data=carriers_SE, aes(x=Timen), color='salmon', linewidth=1, alpha=.5) +
  geom_half_violin(data=carriers %>% filter(Timen=="1"), aes(x = Timen, y = resid_Pc), #CHANGE
                   position=position_nudge(x=1.3),side="r",width=.5, fill="salmon",alpha=.6)+
  geom_half_violin(data=carriers %>% filter(Timen=="2"), aes(x = Timen, y = resid_Pc), #CHANGE
                   position=position_nudge(x=.3),side="r",width=.5, fill="salmon",alpha=.6)+
  scale_x_continuous(breaks = c(1,2),labels=c("Baseline","Follow-up"),limits = c(0,3))+ #,limits = c(0.5,3)
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01),limits=c(-0.035,0.035))+ #,,limits=c(-0.35,-0.25)
  labs(title = "APOE ɛ4 carriers", y="Global Pc")+
  theme(axis.text=element_text(size=12,colour = "black"),
        axis.title=element_text(size=12,colour = "black"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black",linewidth = 0.5),
        legend.position = "none",
        axis.title.x=element_blank())

ggplot(data=noncarriers, aes(y = resid_Pc))+ #CHANGE
  geom_point(data=noncarriers%>%filter(Timen=="1"), aes(x = Timenj), color = "cornflowerblue", size=2, alpha=.4)+
  geom_point(data=noncarriers%>%filter(Timen=="2"), aes(x = Timenj), color = "cornflowerblue", size=2, alpha=.4)+
  geom_line(data=noncarriers, aes(x=Timenj, group=SubjID), color='gray', alpha=.5)+
  geom_errorbar(data=noncarriers_SE, aes(x=Timen, ymin=resid_Pc-ci, ymax=resid_Pc+ci), width=.1) +
  geom_point(data=noncarriers_SE, aes(x=Timen, y=resid_Pc), shape=21, color="black", fill="cornflowerblue", size=2, stroke=1) +
  geom_line(data=noncarriers_SE, aes(x=Timen), color='cornflowerblue', linewidth=1, alpha=.5) +
  geom_half_violin(data=noncarriers %>% filter(Timen=="1"), aes(x = Timen, y = resid_Pc), #CHANGE
                   position=position_nudge(x=-.3),side="l",width=.5, fill="cornflowerblue",alpha=.4)+
  geom_half_violin(data=noncarriers %>% filter(Timen=="2"), aes(x = Timen, y = resid_Pc), #CHANGE
                   position=position_nudge(x=-1.3),side="l",width=.5, fill="cornflowerblue",alpha=.4)+
  scale_x_continuous(breaks = c(1,2),labels=c("Baseline","Follow-up"),limits = c(0,3))+ #,limits = c(0.5,3)
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01),limits=c(-0.035,0.035))+ #,,limits=c(-0.35,-0.25) c(-0.35, -0.10)
  labs(title = "APOE ɛ4 non-carriers", y="Global Pc")+
  theme(axis.text=element_text(size=12,colour = "black"),
        axis.title=element_text(size=12,colour = "black"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black",linewidth = 0.5),
        legend.position = "none",
        axis.title.x=element_blank())
