#LAKANA INFLAMMATION ANALYSES



### Analytical script for LAKANA inflammation data

library(ggplot2)
library(lme4)
library(dplyr)
library(knitr)
library(lmerTest)
library(tidyverse)
library(zscorer)
library(margins)
library(forestplot)

## Save the data
data=CRP_mechanistic_data_2026_03_06


###############################################################################
###############################################################################

### Mixed effects models

###############################################################################
###############################################################################


#CRP ratio of means, general
data$HeelCRP_4_1
data$HeelCRP_4_2
data$VillageID
data$interventionbinary

CRP_mixed.lmer2 <- lmer(log(HeelCRP_4_2) ~ as.factor(interventionbinary) + log(HeelCRP_4_1) + (1|VillageID), data = data)
summary(CRP_mixed.lmer2)
confint(CRP_mixed.lmer2)

"Then use exp() function to get the ratio and 95%CI"



#CRP among those who had elevated CRP at baseline
elevCRPsub=data$`Elevated CRP_1`
efmodelevCRP <- lmer(log(HeelCRP_4_2) ~ as.factor(interventionbinary)*as.factor(elevCRPsub) + log(HeelCRP_4_1) + (1|VillageID), data = data)
car::deltaMethod(efmodelevCRP,paste("b1+b4"),parameterNames=c("b0","b1","b2","b3","b4"))

"use exp() to get the testimate and CIs"

"To find the two-sided p-value, use the below calculation"
z = Estimate / SE
p_value = 2 * (1-pnorm(abs(z)))



#elevated CRP, above 5
CRP_status_1=data$`Elevated CRP_1`
CRP_status_2=data$`Elevated CRP_2`
VillageID=data$VillageID
interventionbinary=data$interventionbinary

elevCRP_mixed.lmer2 <- glmer(as.factor(CRP_status_2) ~ as.factor(interventionbinary) + as.factor(CRP_status_1) + (1|VillageID), family="binomial", glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)), data = data)
summary(elevCRP_mixed.lmer2)
marginaleffects::avg_comparisons(elevCRP_mixed.lmer2, comparison = "lnratio", transform = exp)



#elevated CRP, above 10
CRP_above10_1=data$`CRP above 10_1`
CRP_above10_2=data$`CRP above 10_2`
VillageID=data$VillageID
interventionbinary=data$interventionbinary

CRPabove10_mixed.lmer2 <- glmer(as.factor(CRP_above10_2) ~ as.factor(interventionbinary) + as.factor(CRP_above10_1) + (1|VillageID), family="binomial", data = data)
summary(CRPabove10_mixed.lmer2)
marginaleffects::avg_comparisons(CRPabove10_mixed.lmer2, comparison = "lnratio", transform = exp)



#elevated CRP, above 1
CRP_above_one_1=data$`CRP above one_1`
CRP_above_one_2=data$`CRP above one_2`
VillageID=data$VillageID
interventionbinary=data$interventionbinary

CRPabove10_mixed.lmer2 <- glmer(as.factor(CRP_above_one_2) ~ as.factor(interventionbinary) + as.factor(CRP_above_one_1) + (1|VillageID), family="binomial", data = data)
summary(CRPabove10_mixed.lmer2)
marginaleffects::avg_comparisons(CRPabove10_mixed.lmer2, comparison = "lnratio", transform = exp)



################################
# Comorbidity analyses


#CRP comorbidity subgroups data separately only 2022 data
data=X2022_CRP_comorbidity_analyses_2026_04_23



#Any comorbidity
comorbsub=data$Comorbidity
efmodcomorb <- lmer(log(HeelCRP_4_2) ~ as.factor(interventionbinary)*as.factor(comorbsub) + log(HeelCRP_4_1) + (1|VillageID), data = data)
car::deltaMethod(efmodcomorb,paste("b1+b4"),parameterNames=c("b0","b1","b2","b3","b4"))

"use exp() to get the testimate and CIs"

"To find the two-sided p-value, use the below calculation"
z = Estimate / SE
p_value = 2 * (1-pnorm(abs(z)))



#Underweight and high EE score
UWandEEsub=data$`EE and UW`
efmodUWandEE <- lmer(log(HeelCRP_4_2) ~ as.factor(interventionbinary)*as.factor(UWandEEsub) + log(HeelCRP_4_1) + (1|VillageID), data = data)
car::deltaMethod(efmodUWandEE,paste("b1+b4"),parameterNames=c("b0","b1","b2","b3","b4"))

"CIs and p-value as above"



#Underweight and severe anemia
UWandSAsub=data$`SA and UW`
efmodUWandSA <- lmer(log(HeelCRP_4_2) ~ as.factor(interventionbinary)*as.factor(UWandSAsub) + log(HeelCRP_4_1) + (1|VillageID), data = data)
car::deltaMethod(efmodUWandSA,paste("b1+b4"),parameterNames=c("b0","b1","b2","b3","b4"))

"CIs and p-value as above"



#Underweight and malaria
UWandMsub=data$`M and UW`
efmodUWandM <- lmer(log(HeelCRP_4_2) ~ as.factor(interventionbinary)*as.factor(UWandMsub) + log(HeelCRP_4_1) + (1|VillageID), data = data)
car::deltaMethod(efmodUWandM,paste("b1+b4"),parameterNames=c("b0","b1","b2","b3","b4"))

"CIs and p-value as above"



#High EE score and severe anemia
EEandSAsub=data$`EE and SA`
efmodEEandSA <- lmer(log(HeelCRP_4_2) ~ as.factor(interventionbinary)*as.factor(EEandSAsub) + log(HeelCRP_4_1) + (1|VillageID), data = data)
car::deltaMethod(efmodEEandSA,paste("b1+b4"),parameterNames=c("b0","b1","b2","b3","b4"))



#High EE score and malaria
EEandMsub=data$`M and EE`
efmodEEandM <- lmer(log(HeelCRP_4_2) ~ as.factor(interventionbinary)*as.factor(EEandMsub) + log(HeelCRP_4_1) + (1|VillageID), data = data)
car::deltaMethod(efmodEEandM,paste("b1+b4"),parameterNames=c("b0","b1","b2","b3","b4"))



#Malaria and severe anemia
SAandMsub=data$`SA and M`
efmodSAandM <- lmer(log(HeelCRP_4_2) ~ as.factor(interventionbinary)*as.factor(SAandMsub) + log(HeelCRP_4_1) + (1|VillageID), data = data)
car::deltaMethod(efmodSAandM,paste("b1+b4"),parameterNames=c("b0","b1","b2","b3","b4"))



#High EE, malaria, and underweight
EEandMandUWsub=data$`EE and M and UW`
efmodEEandMandUW <- lmer(log(HeelCRP_4_2) ~ as.factor(interventionbinary)*as.factor(EEandMandUWsub) + log(HeelCRP_4_1) + (1|VillageID), data = data)
car::deltaMethod(efmodEEandMandUW,paste("b1+b4"),parameterNames=c("b0","b1","b2","b3","b4"))



#High EE, severe anemia, and underweight
EEandSAandUWsub=data$`EE and SA and UW`
efmodEEandSAandUW <- lmer(log(HeelCRP_4_2) ~ as.factor(interventionbinary)*as.factor(EEandSAandUWsub) + log(HeelCRP_4_1) + (1|VillageID), data = data)
car::deltaMethod(efmodEEandSAandUW,paste("b1+b4"),parameterNames=c("b0","b1","b2","b3","b4"))



#Malaria, severe anemia, and underweight
MandSAandUWsub=data$`SA and M and UW`
efmodMandSAandUW <- lmer(log(HeelCRP_4_2) ~ as.factor(interventionbinary)*as.factor(MandSAandUWsub) + log(HeelCRP_4_1) + (1|VillageID), data = data)
car::deltaMethod(efmodMandSAandUW,paste("b1+b4"),parameterNames=c("b0","b1","b2","b3","b4"))



#Elevated EE, malaria, severe anemia, and underweight
EEandSAandUWandMsub=data$`EE and SA and M and UW`
efmodEEandSAandUWandM <- lmer(log(HeelCRP_4_2) ~ as.factor(interventionbinary)*as.factor(EEandSAandUWandMsub) + log(HeelCRP_4_1) + (1|VillageID), data = data)
car::deltaMethod(efmodEEandSAandUWandM,paste("b1+b4"),parameterNames=c("b0","b1","b2","b3","b4"))




#################################

#Effect modifiers for CRP

#Age
AgeAbove8 = data$`ageabove8efmod`
efmodAge <- lmer(log(HeelCRP_4_2) ~ as.factor(interventionbinary)*as.factor(AgeAbove8) + log(HeelCRP_4_1) + (1|VillageID), data = data)
summary(efmodAge)

car::deltaMethod(efmodAge,paste("b1+b4"),parameterNames=c("b0","b1","b2","b3","b4"))
car::deltaMethod(efmodAge,paste("b1"),parameterNames=c("b0","b1","b2","b3","b4"))



#Sex
Male=data$Male

efmodMale <- lmer(log(HeelCRP_4_2) ~ as.factor(interventionbinary)*as.factor(Male) + log(HeelCRP_4_1) + (1|VillageID), data = data)
summary(efmodMale)

car::deltaMethod(efmodMale,paste("b1+b4"),parameterNames=c("b0","b1","b2","b3","b4"))
car::deltaMethod(efmodMale,paste("b1"),parameterNames=c("b0","b1","b2","b3","b4"))



#WAZ
WAZordinal=data$`WAZ_ordinal_low_to_high`

efmodWAZordinal <- lmer(log(HeelCRP_4_2) ~ as.factor(interventionbinary)*as.factor(WAZordinal) + log(HeelCRP_4_1) + (1|VillageID), data = data)
summary(efmodWAZordinal)

car::deltaMethod(efmodWAZordinal,paste("b1"),parameterNames=c("b0","b1","b2","b3","b4","b5","b6"))
car::deltaMethod(efmodWAZordinal,paste("b1+b6"),parameterNames=c("b0","b1","b2","b3","b4","b5","b6"))
car::deltaMethod(efmodWAZordinal,paste("b1+b5"),parameterNames=c("b0","b1","b2","b3","b4","b5","b6"))



#Season
year2022=data$`Yes 2022`

efmodtime <- lmer(log(HeelCRP_4_2) ~ as.factor(interventionbinary)*as.factor(year2022) + log(HeelCRP_4_1) + (1|VillageID), data = data)
summary(efmodtime)

car::deltaMethod(efmodtime,paste("b1+b4"),parameterNames=c("b0","b1","b2","b3","b4"))
car::deltaMethod(efmodtime,paste("b1"),parameterNames=c("b0","b1","b2","b3","b4"))



#WASH
WASHabove = data$`WASH above and equal`

efmodWASH <- lmer(log(HeelCRP_4_2) ~ as.factor(interventionbinary)*as.factor(WASHabove) + log(HeelCRP_4_1) + (1|VillageID), data = data)
summary(efmodWASH)

car::deltaMethod(efmodWASH,paste("b1+b4"),parameterNames=c("b0","b1","b2","b3","b4"))
car::deltaMethod(efmodWASH,paste("b1"),parameterNames=c("b0","b1","b2","b3","b4"))



#Asset
Assetabove = data$`Asset above and equal`

efmodAsset <- lmer(log(HeelCRP_4_2) ~ as.factor(interventionbinary)*as.factor(Assetabove) + log(HeelCRP_4_1) + (1|VillageID), data = data)
summary(efmodAsset)

car::deltaMethod(efmodAsset,paste("b1+b4"),parameterNames=c("b0","b1","b2","b3","b4"))
car::deltaMethod(efmodAsset,paste("b1"),parameterNames=c("b0","b1","b2","b3","b4"))



#Any comorbidity versus no comorbidity
comorbidity=data$Comorbidity

efmodcomorb <- lmer(log(HeelCRP_4_2) ~ as.factor(interventionbinary)*as.factor(comorbidity) + log(HeelCRP_4_1) + (1|VillageID), data = data)
summary(efmodcomorb)
coef(summary(efmodcomorb))

car::deltaMethod(efmodcomorb,paste("b1+b4"),parameterNames=c("b0","b1","b2","b3","b4"))
car::deltaMethod(efmodcomorb,paste("b1"),parameterNames=c("b0","b1","b2","b3","b4"))






###############################################################################

###############################################################################

#p-for-interactions calculations for CRP effect modification


#Age
AgeAbove8 = data$`ageabove8efmod`
mod_null = lmer(log(HeelCRP_4_2) ~ as.factor(interventionbinary) + as.factor(AgeAbove8) + log(HeelCRP_4_1) + (1|VillageID), data = data)
mod_full = lmer(log(HeelCRP_4_2) ~ as.factor(interventionbinary) * as.factor(AgeAbove8) + log(HeelCRP_4_1) + (1|VillageID), data = data)
anova(mod_null,mod_full, test = "LRT")


#Sex
Male=data$Male
mod_null = lmer(log(HeelCRP_4_2) ~ log(HeelCRP_4_1) + as.factor(interventionbinary) + as.factor(Male) + (1|VillageID), data = data)
mod_full = lmer(log(HeelCRP_4_2) ~ log(HeelCRP_4_1) + as.factor(interventionbinary) * as.factor(Male) + (1|VillageID), data = data)
anova(mod_null,mod_full, test = "LRT")


#WAZ 
WAZordinal=data$`WAZ_ordinal_low_to_high`
mod_null = lmer(log(HeelCRP_4_2) ~ log(HeelCRP_4_1) + as.factor(interventionbinary) + as.factor(WAZordinal) + (1|VillageID), data = data)
mod_full = lmer(log(HeelCRP_4_2) ~ log(HeelCRP_4_1) + as.factor(interventionbinary) * as.factor(WAZordinal) + (1|VillageID), data = data)
anova(mod_null,mod_full, test = "LRT")


#Season
year2022=data$`Yes 2022`
mod_null = lmer(log(HeelCRP_4_2) ~ log(HeelCRP_4_1) + as.factor(interventionbinary) + as.factor(year2022) + (1|VillageID), data = data)
mod_full = lmer(log(HeelCRP_4_2) ~ log(HeelCRP_4_1) + as.factor(interventionbinary) * as.factor(year2022) + (1|VillageID), data = data)
anova(mod_null,mod_full, test = "LRT")


#WASH
WASHabove = data$`WASH above and equal`
mod_null = lmer(log(HeelCRP_4_2) ~ log(HeelCRP_4_1) + as.factor(interventionbinary) + as.factor(WASHabove) + (1|VillageID), data = data)
mod_full = lmer(log(HeelCRP_4_2) ~ log(HeelCRP_4_1) + as.factor(interventionbinary) * as.factor(WASHabove) + (1|VillageID), data = data)
anova(mod_null,mod_full, test = "LRT")


#Asset
Assetabove = data$`Asset above and equal`
mod_null = lmer(log(HeelCRP_4_2) ~ log(HeelCRP_4_1) + as.factor(interventionbinary) + as.factor(Assetabove) + (1|VillageID), data = data)
mod_full = lmer(log(HeelCRP_4_2) ~ log(HeelCRP_4_1) + as.factor(interventionbinary) * as.factor(Assetabove) + (1|VillageID), data = data)
anova(mod_null,mod_full, test = "LRT")


#Comorbidity
comorbidity=data$Comorbidity
mod_null = lmer(log(HeelCRP_4_2) ~ log(HeelCRP_4_1) + as.factor(interventionbinary) + as.factor(comorbidity) + (1|VillageID), data = data)
mod_full = lmer(log(HeelCRP_4_2) ~ log(HeelCRP_4_1) + as.factor(interventionbinary) * as.factor(comorbidity) + (1|VillageID), data = data)
anova(mod_null,mod_full, test = "LRT")





###############################################################################
###############################################################################

##Density plots

###############################################################################
###############################################################################

# Libraries
library(ggplot2)
library(hrbrthemes)
library(dplyr)
library(tidyr)
library(viridis)

LogCRPplacebo <- melt(df)

# With transparency (right)
ggplot(logCRPplacebo, aes(x=value, fill=variable)) +
  geom_density(alpha=.25) + theme_classic() + ylim(0, 0.7) +
  xlim(-1.5, 5)

ggplot(logCRPazithromycin, aes(x=value, fill=variable)) +
  geom_density(alpha=.25) + theme_classic() + ylim(0, 0.7) +
  xlim(-1.5, 5)






###############################################################################
###############################################################################

## Forest plots

###############################################################################
###############################################################################


#CRP forest plot, numbers manually coded



base_data <- tibble::tibble(mean  = c(NA, 
                                      NA, 0.8056451, 0.9215926, 
                                      NA, 
                                      NA, 0.9387191, 0.7892683, 
                                      NA,
                                      NA, 0.8275217, 0.8849917, 0.8359136,
                                      NA,
                                      NA, 0.8085765, 0.9319175,
                                      NA, 
                                      NA, 0.8550973, 0.8684459,
                                      NA,
                                      NA, 0.8876764, 0.639384,
                                      NA,
                                      NA, 0.5750605, 0.9379178),
                            lower = c(NA, 
                                      NA, 0.6734537, 0.7679528, 
                                      NA, 
                                      NA, 0.7835644, 0.6605331, 
                                      NA, 
                                      NA, 0.6345622, 0.7510238, 0.6455517,
                                      NA,
                                      NA, 0.6570665, 0.7586422,
                                      NA,
                                      NA, 0.6993108, 0.7112011,
                                      NA,
                                      NA, 0.7430559, 0.4530467,
                                      NA,
                                      NA, 0.3332608, 0.7747468),
                            upper = c(NA, 
                                      NA, 0.9637725, 1.105945, 
                                      NA,
                                      NA, 1.124569, 0.9430839, 
                                      NA,
                                      NA, 1.079178, 1.042894, 1.082421,
                                      NA,
                                      NA, 0.9950125, 1.144766,
                                      NA,
                                      NA, 1.04561, 1.060457,
                                      NA,
                                      NA, 1.060457, 0.9023977,
                                      NA,
                                      NA, 0.9923296, 1.135417),
                            Variable = c("",
                                         "  Age in months (p=0.21) ","4 to 7","8 to 11",
                                         "",
                                         "  Sex (p=0.10)","Female","Male",
                                         "",
                                         "  WAZ (p=0.84)","<2-SD","≥-2SD and <0SD","≥0SD",
                                         "",
                                         "  WASH index (p=0.30)", "Above or equal to median", "Below median",
                                         "",
                                         "  Asset index (p=0.91)", "Above or equal to median", "Below median",
                                         "",
                                         "  Place and timing (p=0.09)", "Kita area, Feb-May 2022", "Kati area, Jul-Aug 2023",
                                         "",
                                         "  Comorbidities, only 2022 Kita sample (p=0.09)", "Two or more known morbidities", "No known morbidites"),
                            Control = c("",
                                        "","   1.70   ","   1.64   ",
                                        "",
                                        "","   1.61   ","   1.72   ",
                                        "",
                                        "","   1.84   ","   1.66   ","   1.52   ",
                                        "",
                                        "","   1.69   ","   1.64   ",
                                        "",
                                        "","   1.62   ","   1.71   ",
                                        "",
                                        "","   1.69   ","   1.65   ",
                                        "",
                                        "","   2.99   ","   1.59   "),
                            Azithromycin = c("",
                                             "","   1.44   ","   1.50   ",
                                             "",
                                             "","   1.54   ","   1.40   ",
                                             "",
                                             "","   1.56   ","   1.49   ","   1.34   ",
                                             "",
                                             "","   1.43   ","   1.53   ",
                                             "",
                                             "","   1.44   ","   1.50   ",
                                             "",
                                             "","   1.52   ","   1.08   ",
                                             "",
                                             "","   1.63   ","   1.51   "),
                            Effectsize = c("",
                                           "","0.81 (0.67 to 0.96)","0.92 (0.77 to 1.11)",
                                           "",
                                           "","0.94 (0.78 to 1.12)", "0.79 (0.66 to 0.94)",
                                           "",
                                           "","0.83 (0.63 to 1.08)","0.88 (0.75 to 1.04)","0.84 (0.65 to 1.08)",
                                           "",
                                           "","0.81 (0.66 to 1.00)","0.93 (0.76 to 1.14)",
                                           "",
                                           "","0.86 (0.70 to 1.05)","0.87 (0.71 to 1.06)",
                                           "",
                                           "","0.89 (0.74 to 1.06)","0.64 (0.45 to 0.90)",
                                           "",
                                           "","0.58 (0.33 to 0.99)","0.94 (0.77 to 1.14)"))
                                            
                            base_data |>
                              forestplot(labeltext = c(Variable, Control, Azithromycin, Effectsize),
                                         clip = c(0.48, 1.2),
                                         vertices = TRUE,
                                         xlog = TRUE,
                                         boxsize = 0.4,
                                         graphwidth = unit(40, 'mm'),
                                         is.summary = c(TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, 
                                                        TRUE, TRUE, FALSE, FALSE, FALSE,
                                                        TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE),
                                         xlab = "                                                                  
                               Favours azithromycin       Favors placebo",) |>
                              fp_decorate_graph(graph.pos = 4) |>
                              fp_add_lines(h_3 = gpar(lty = 1),
                                           h_31 = gpar(lwd = 1, columns = 1:5)) |>
                              fp_set_style(box = "black",
                                           line = "black") |>
                              fp_add_header(Variable = c("Variable (p-for-interactions)"),
                                            Control = c("Control"),
                                            Azithromycin = c("Azithromycin"),
                                            Effectsize = c("Ratio (95% CI)"))


                            