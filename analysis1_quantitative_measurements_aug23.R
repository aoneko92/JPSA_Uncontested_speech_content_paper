##LOAD ANALYSIS PACKAGES
library(lmtest)
library(sandwich)
library(clubSandwich)
library(multiwaycov)
library(MASS)
library(texreg)
library(alpaca)

library(pglm)
library(vtable)
library(ggsci)
library(rlang)
library(tidyverse)

require(pscl)
require(boot)
library(stargazer)


library(readr)
#### ANALYSIS 1: Quantitative measurement #### 

reg_data <- read_csv("regular_analysis_period_aug23.csv")
did_data <- read_csv("did_analysis_aug23.csv")


###quick check
length(unique(reg_data$name))
reg_data %>% summarize()

###fix some variables
reg_data$party <- ifelse(!reg_data$party %in% c("自","共","無","公","民","社","維"), "他", reg_data$party)
reg_data <- within(reg_data, party <- relevel(as.factor(party), ref = "自")) 

did_data$after <- rep(0:1, each=5, length.out = nrow(did_data))
did_data <-did_data %>% group_by(id) %>% mutate(min_reelect = min(reelect))


##load governors party
{
  gov_2 <- read.csv("gov_period.csv")
  gov_2 <- gov_2[!is.na(gov_2$period),]
  
  gov_2$gov_year
  #gov_2 <- rename(gov_2, year = period)
  reg_data <- merge(reg_data, gov_2, by = c("pref_id", "period"))
  
  reg_data$party_2 <- ifelse(!reg_data$party %in% c("自","共","無","公","民","社","維"), "他", reg_data$party)
  ###create a governor ruling party dummy
  reg_data$gov_suisen <- ifelse(is.na(reg_data$gov_suisen), "None",reg_data$gov_suisen)
  reg_data$rul_1 <- reg_data$party %in% reg_data$gov_suisen
  reg_data$ainori <- ifelse(is.na(reg_data$ainori), 0,reg_data$ainori)
  reg_data$rul_2 <- reg_data$ainori == 1 & (reg_data$party %in% c("自","無","公","民","維"))
  reg_data$anti_ldp <- ifelse(is.na(reg_data$anti_ldp), 0,reg_data$anti_ldp)
  reg_data$rul_3 <- reg_data$anti_ldp == 1 & (reg_data$party %in% c("無","民","社"))
  reg_data$rul_4 <- reg_data$anti_ldp == 0 & (reg_data$party %in% c("自","公"))
  
  
  reg_data$gov_ruling <- ifelse(reg_data$rul_2 == T| reg_data$rul_3 == T| reg_data$rul_4 == T, T, F)
  
  #gov_2 <- rename(gov_2, year = period)
  did_data <- merge(did_data, gov_2, by = c("pref_id", "period"))
  
  did_data$party_2 <- ifelse(!did_data$party %in% c("自","共","無","公","民","社","維"), "他", did_data$party)
  ###create a governor ruling party dummy
  did_data$gov_suisen <- ifelse(is.na(did_data$gov_suisen), "None",did_data$gov_suisen)
  did_data$rul_1 <- did_data$party %in% did_data$gov_suisen
  did_data$ainori <- ifelse(is.na(did_data$ainori), 0,did_data$ainori)
  did_data$rul_2 <- did_data$ainori == 1 & (did_data$party %in% c("自","無","公","民","維"))
  did_data$anti_ldp <- ifelse(is.na(did_data$anti_ldp), 0,did_data$anti_ldp)
  did_data$rul_3 <- did_data$anti_ldp == 1 & (did_data$party %in% c("無","民","社"))
  did_data$rul_4 <- did_data$anti_ldp == 0 & (did_data$party %in% c("自","公"))
  
  
  did_data$gov_ruling <- ifelse(did_data$rul_2 == T| did_data$rul_3 == T| did_data$rul_4 == T, T, F)
}

### set dependent variables
dvs <-  c("length", "N", "ent", "length_committee", "N_committee")





### regular models


mods <- c()
for (i in 1:length(dvs)){
  dv <- dvs[i]
  #if(dv %in% c("length"))
  if (i %in% c(1,4) ){
  mods[[i]]  <- lm(get(dv) ~ mutohyo + reelect + totseat + secretary + city_assembly+ gov_ruling + party + age + factor(pref_id) + factor(period)
                        , data= reg_data)
  } else {
    mods[[i]]  <- glm.nb(get(dv) ~ mutohyo + reelect + totseat + secretary + city_assembly+ gov_ruling + party + age + factor(pref_id) + factor(period)
                     , data= reg_data)
  }
  }

dvs_name <-  c("Speech length", "No. speeches", "Topic diversity", "Speech length (committees)", "No. speeches (committees)")


screenreg(list(mods[[1]],mods[[2]],mods[[3]],mods[[4]],mods[[5]]),  type = "clustered",
          cluster = ~pref_id + period , stars = c(0.01, 0.05, 0.1, 0.15),
          custom.coef.map =list('mutohyo'="Uncontested", "reelect" = "Tenure",
                                "totseat" = "Magnitude", "(Intercept)" = "Intercept" ),
          custom.gof.rows = list("Party dummies" = c("YES","YES","YES","YES","YES"),
                                 "Control" = c("YES","YES","YES","YES","YES"),
                                 "Time fixed effect" = c("YES","YES","YES","YES","YES"),
                                 "Prefecture effect" = c("YES","YES","YES","YES","YES"),
                                 " " = c(" "," "," "," "," ")),
          custom.note = paste("%stars.",
                              "Regression models 1 and 4 were estimated using OLS regression. Models 2,3 and 5 were estimated using a negative ",
                              "binomial distribution. All models use clustered standard errors by time period and prefecture. Clustered standard",
                              "errors are in parenthesis. The unit of analysis is a legislators total activity during a 4-year electoral period ",
                              "(for example 2003 - 2007). Time period was 2003 - 2019 with a total 4 distinct electoral periods included. N = 2682",
                              "in all models"),
          caption = "Table 2: Regular regression models for effect on quantitative measurements",
          head.tag = TRUE,
          custom.model.names=c(dvs_name[1], dvs_name[2], dvs_name[3], dvs_name[4], dvs_name[5]),
          include.loglik = T)
htmlreg(list(mods[[1]],mods[[2]],mods[[3]],mods[[4]],mods[[5]]),  type = "clustered",
          cluster = ~pref_id + period , stars = c(0.01, 0.05, 0.1, 0.15),
          custom.coef.map =list('mutohyo'="Uncontested", "reelect" = "Tenure",
                                "totseat" = "Magnitude", "(Intercept)" = "Intercept" ),
          custom.gof.rows = list("Party dummies" = c("YES","YES","YES","YES","YES"),
                                 "Control" = c("YES","YES","YES","YES","YES"),
                                 "Time fixed effect" = c("YES","YES","YES","YES","YES"),
                                 "Prefecture effect" = c("YES","YES","YES","YES","YES"),
                                 " " = c(" "," "," "," "," ")),
          custom.note = paste("%stars.",
                            "Regression models 1 and 4 were estimated using OLS regression. Models 2,3 and 5 were estimated using a negative ",
                            "binomial distribution. All models use clustered standard errors by time period and prefecture. Clustered standard",
                            "errors are in parenthesis. The unit of analysis is a legislators total activity during a 4-year electoral period ",
                            "(for example 2003 - 2007). Time period was 2003 - 2019 with a total 4 distinct electoral periods included. N = 2682",
                            "in all models"),
          caption = "Table 2: Regular regression models for effect on quantitative measurements",
          head.tag = TRUE,
          custom.model.names=c(dvs_name[1], dvs_name[2], dvs_name[3], dvs_name[4], dvs_name[5]),
          include.loglik = T, file = "C:/Users/robba/Dropbox/uncontested_text_analysis_paper/regression_models/table_2_sep1.html")
htmlreg(list(mods[[1]],mods[[2]],mods[[3]],mods[[4]],mods[[5]]),  type = "clustered",
        cluster = ~pref_id + period , stars = c(0.01, 0.05, 0.1, 0.15),
        custom.coef.map =list('mutohyo'="Uncontested", "reelect" = "Tenure",
                              "totseat" = "Magnitude", "(Intercept)" = "Intercept" ),
        custom.gof.rows = list("Party dummies" = c("YES","YES","YES","YES","YES"),
                               "Control" = c("YES","YES","YES","YES","YES"),
                               "Time fixed effect" = c("YES","YES","YES","YES","YES"),
                               "Prefecture effect" = c("YES","YES","YES","YES","YES"),
                               " " = c(" "," "," "," "," ")),
        custom.note = paste("%stars.",
                            "Regression models 1 and 4 were estimated using OLS regression. Models 2,3 and 5 were estimated using a negative ",
                            "binomial distribution. All models use clustered standard errors by time period and prefecture. Clustered standard",
                            "errors are in parenthesis. The unit of analysis is a legislators total activity during a 4-year electoral period ",
                            "(for example 2003 - 2007). Time period was 2003 - 2019 with a total 4 distinct electoral periods included. N = 2682",
                            "in all models"),
        caption = "Table 2: Regular regression models for effect on quantitative measurements",
        head.tag = TRUE,
        custom.model.names=c(dvs_name[1], dvs_name[2], dvs_name[3], dvs_name[4], dvs_name[5]),
        include.loglik = T, file = "C:/Users/robba/Dropbox/uncontested_text_analysis_paper/regression_models/table_2_sep1.doc")

### full regression table
htmlreg(list(mods[[1]],mods[[2]],mods[[3]],mods[[4]],mods[[5]]),  type = "clustered",
          cluster = ~pref_id + period , stars = c(0.01, 0.05, 0.1, 0.15),
          custom.note = paste("%stars.",
                              "Regression models 1 and 4 were estimated using OLS regression. Models 2,3 and 5 were estimated using a negative ",
                              "binomial distribution. All models use clustered standard errors by time period and prefecture. Clustered standard",
                              "errors are in parenthesis. The unit of analysis is a legislators total activity during a 4-year electoral period ",
                              "(for example 2003 - 2007). Time period was 2003 - 2019 with a total 4 distinct electoral periods included. N = 2682",
                              "in all models"),
          caption = "Table 2: Regular regression models for effect on quantitative measurements",
          head.tag = TRUE,
          custom.model.names=c(dvs_name[1], dvs_name[2], dvs_name[3], dvs_name[4], dvs_name[5]),
          include.loglik = T, file = "C:/Users/robba/Dropbox/uncontested_text_analysis_paper/regression_models/table_2_full_sep1.html")



###classic DiD

require(ggplot2)
require(pscl)
require(boot)
summary(lm(length_committee ~ treat*after, data= did_data))

summary(feglm.nb(N_committee ~ mutohyo + gov_ruling| id + year, data = did_data))
summary(feglm.nb(N_committee ~ mutohyo |min_reelect +  id + year, data = did_data))


library(glmmTMB)
require(broom)
require(dotwhisker)
library(broom.mixed)

eq <- length ~ factor(mutohyo) + (1|year + id)

tmb2 <- glmmTMB(eq, data=did_data,ziformula=~1, family=gaussian)

summary(tmb2)
??glmmTMB::familyglmmTMB

### generalized DiD
dvs <-  c("length", "N", "ent_high", "length_committee", "N_committee")


### perhaps try zero-inflation also 
library(fixest)
mods <- c()
for (i in 1:length(dvs)){
  print(i)
  a_data <- did_data %>% group_by(id) %>% filter(min(period) == 2)
  a_data <- a_data %>% group_by(id) %>% filter(max(year) <= min(year) + 8)
  dv <- dvs[i]
  #if(dv %in% c("length"))
  if (i %in% c(1,4) ){
    res <- lapply(dv, function(var) {
      
      #eq <- xpd(..lhs ~ mutohyo+ (1|year + id), ..lhs = var)
      #res <- glmmTMB(eq, data=a_data, ziformula=~1, family=gaussian)
      print(summary(res))
      res <- feols(xpd(..lhs ~ mutohyo+ gov_ruling| id + year, ..lhs = var), data = a_data)
      return(res)
      })  


  } else {
    res <- lapply(dv, function(var) {
      #eq <- xpd(..lhs ~ mutohyo+ (1|year + id), ..lhs = var)
      #res <- glmmTMB(eq, data=a_data, ziformula=~1, family=nbinom1)
      print(summary(res))
      res <- feglm.nb(xpd(..lhs ~ mutohyo+ gov_ruling| id + year, ..lhs = var), data = a_data)
      return(res)
          }) 
  }
  mods[[i]] <- res[[1]]
}
summary(mods[[1]])
mean()
screenreg(list(mods[[1]],mods[[2]],mods[[3]],mods[[4]],mods[[5]]),  type = "clustered", cluster = ~id , stars = c(0.01, 0.05, 0.1, 0.15),
          custom.coef.map =list('mutohyo'="Uncontested", 'gov_rulingTRUE' = "Aligned with governor"),
          custom.gof.rows = list("Time fixed effect" = c("YES","YES","YES","YES","YES"),
                                 " " = c(" "," "," "," "," ")),
          custom.note = paste("%stars.",
                              "Regression models 1 and 4 were estimated using OLS regression. Models 2,3 and 5 were estimated using a negative ",
                              "binomial distribution. All models use clustered standard errors. Clustered standard",
                              "errors are in parenthesis. Unit of analysis is the activity of one legislator in one year. To avoid estimating  effect on the same ",
                              "legislator appearing in different time periods, only one 8-year period is used in analysis,",
                              "here 2003 - 2011 where 2003 - 2007 (before election) is before the treatment and 2007 (after election) - 2011 is the ",
                              "treated period (regression results for other time periods in supplementary files). ",
                              "Groups: id = No unique legislators, Groups: year = unique years. "),
          caption = "Table 3: Difference in difference models for effect on quantitative measurements",
          head.tag = TRUE,
          custom.model.names=c(dvs_name[1], dvs_name[2], dvs_name[3], dvs_name[4], dvs_name[5]),
          include.loglik = T)
htmlreg(list(mods[[1]],mods[[2]],mods[[3]],mods[[4]],mods[[5]]),  type = "clustered", cluster = ~id , stars = c(0.01, 0.05, 0.1, 0.15),
        custom.coef.map =list('mutohyo'="Uncontested", 'gov_rulingTRUE' = "Aligned with governor"),
        custom.gof.rows = list("Time fixed effect" = c("YES","YES","YES","YES","YES"),
                               "Prefecture fixed effect" = c("YES","YES","YES","YES","YES"),
                               " " = c(" "," "," "," "," ")),
        custom.note = paste("%stars.",
                            "Regression models 1 and 4 were estimated using OLS regression. Models 2,3 and 5 were estimated using a negative ",
                            "binomial distribution. All models use clustered standard errors. Clustered standard",
                            "errors are in parenthesis. Unit of analysis is the activity of one legislator in one year. To avoid estimating  effect on the same ",
                            "legislator appearing in different time periods, only one 8-year period is used in analysis,",
                            "here 2003 - 2011 where 2003 - 2007 (before election) is before the treatment and 2007 (after election) - 2011 is the ",
                            "treated period (regression results for other time periods in supplementary files). ",
                            "Groups: id = No unique legislators, Groups: year = unique years. "),
        caption = "Table 3: Difference in difference models for effect on quantitative measurements",
        head.tag = TRUE,
        custom.model.names=c(dvs_name[1], dvs_name[2], dvs_name[3], dvs_name[4], dvs_name[5]),
        include.loglik = T, file = "C:/Users/robba/Dropbox/uncontested_text_analysis_paper/regression_models/table_3_aug24.html")
htmlreg(list(mods[[1]],mods[[2]],mods[[3]],mods[[4]],mods[[5]]),  type = "clustered", cluster = ~id , stars = c(0.01, 0.05, 0.1, 0.15),
        custom.coef.map =list('mutohyo'="Uncontested", 'gov_rulingTRUE' = "Aligned with governor"),
        custom.gof.rows = list("Time fixed effect" = c("YES","YES","YES","YES","YES"),
                               "Prefecture fixed effect" = c("YES","YES","YES","YES","YES"),
                               " " = c(" "," "," "," "," ")),
        custom.note = paste("%stars.",
                            "Regression models 1 and 4 were estimated using OLS regression. Models 2,3 and 5 were estimated using a negative ",
                            "binomial distribution. All models use clustered standard errors. Clustered standard",
                            "errors are in parenthesis. Unit of analysis is the activity of one legislator in one year. To avoid estimating  effect on the same ",
                            "legislator appearing in different time periods, only one 8-year period is used in analysis,",
                            "here 2003 - 2011 where 2003 - 2007 (before election) is before the treatment and 2007 (after election) - 2011 is the ",
                            "treated period (regression results for other time periods in supplementary files). ",
                            "Groups: id = No unique legislators, Groups: year = unique years. "),
        caption = "Table 3: Difference in difference models for effect on quantitative measurements",
        head.tag = TRUE,
        custom.model.names=c(dvs_name[1], dvs_name[2], dvs_name[3], dvs_name[4], dvs_name[5]),
        include.loglik = T, file = "C:/Users/robba/Dropbox/uncontested_text_analysis_paper/regression_models/table_3_aug24.doc")







