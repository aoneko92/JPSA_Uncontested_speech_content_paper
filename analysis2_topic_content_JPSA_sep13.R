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
library(fixest)

library(readr)
#### ANALYSIS 1: Quantitative measurement #### 

reg_data <- read_csv("regular_analysis_period_aug23.csv")
did_data <- read_csv("did_analysis_aug23.csv")
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
dvs <-  c("other_1", "other_2","other_3","other_4","genera_admin",
          "disabilities","nuclear_energy","natural_disaster",
          "forestry_fishing_environment_renewables" ,"hospitals_patients",
          "police_crime","tax_insurance","women_work_children",
          "water_sewer","education","agriculture","dams_rivers",
          "sports_events","railroads","cancer_health_patients",
          "roads_infrastructure","national_politics_war_peace",
          "airports_ports_communications","municipalities_mergers_admin")


### Regular regression models

mods <- c()
coef_mutohyo <- data.frame()
d_2 <- data.frame()
d_3 <- data.frame()
d_4 <- data.frame()
d_5 <- data.frame()


for (i in 1:length(dvs)){
  dv <- dvs[i]
  #if(dv %in% c("length"))
  res  <- lm(get(dv) ~ mutohyo + reelect + totseat + secretary + city_assembly + party + gov_ruling + age + factor(pref_id) + factor(period)
                   , data= reg_data)
  ### get mutohyo values
  coef <- res$coefficients[2]
  res_cl <- coeftest(res, vcov = vcovCL, cluster = ~pref_id + period)
  p_value <- res_cl[2,4]
  #as.data.frame(coef, p_value)
  d <- data.frame(coef = coef, p = p_value, topic = dvs[i], i = i)
  d$sig <- int(d$p <.1)
  d$sig <-ifelse(is.na(d$sig), 0 ,d$sig)
  
  coef_mutohyo <- rbind(d,coef_mutohyo)
  
  ### get reelect values
  summary(res)
  coef <- res$coefficients[3]
  res_cl <- coeftest(res, vcov = vcovCL, cluster = ~pref_id + period)
  p_value <- res_cl[3,4]
  #as.data.frame(coef, p_value)
  d <- data.frame(coef = coef, p = p_value, topic = dvs[i])
  d$sig <- int(d$p <.1)
  d$sig <-ifelse(is.na(d$sig), 0 ,d$sig)
  
  d_2 <- rbind(d,d_2)
  
  ### get totseat values
  summary(res)
  coef <- res$coefficients[4]
  res_cl <- coeftest(res, vcov = vcovCL, cluster = ~pref_id + period)
  p_value <- res_cl[4,4]
  #as.data.frame(coef, p_value)
  d <- data.frame(coef = coef, p = p_value, topic = dvs[i])
  d$sig <- int(d$p <.1)
  d$sig <-ifelse(is.na(d$sig), 0 ,d$sig)
  
  d_3 <- rbind(d,d_3)
  
  ### get r2 value
  coef <- summary(res)$adj.r.squared
  p_value <- NA
  d <- data.frame(coef = coef, p = p_value, topic = dvs[i])
  d$sig <- NA
  d_4 <- rbind(d,d_4)
  
  ### get topic mean
  coef <- mean(reg_data[,dvs[i]], na.rm = TRUE)
  p_value <- NA
  d <- data.frame(coef = coef, p = p_value, topic = dvs[i])
  d$sig <- NA
  d_5 <- rbind(d,d_5)
  
  mods[[i]] <- res
  }
coef_mutohyo$x <- 1
d_2$x <- 2
d_3$x <- 3
## x = 4 is the column for control
d_4$x <- 5
d_5$x <- 6

coef_mutohyo <- coef_mutohyo %>% arrange(coef) 
order <- c(coef_mutohyo$topic)

coef_mutohyo <- coef_mutohyo %>%mutate( max = max(abs(coef_mutohyo$coef)),
                        val = (coef/max)/2+.5  )
d_2 <- d_2 %>%mutate( max = max(abs(d_2$coef)),
                                        val = (coef/max)/2+.5  )
d_3 <- d_3 %>%mutate( max = max(abs(d_3$coef)),
                                        val = (coef/max)/2+.5  )
d_4 <- d_4 %>%mutate( max = max(abs(d_4$coef)),
                                        val = (coef/max)/2+.5  )
d_5 <- d_5 %>%mutate( max = max(abs(d_5$coef)),
                      val = (coef/max)/2+.5  )

p_data <- rbind(coef_mutohyo, d_2,d_3,d_4,d_5)

p_data$V2 <- ifelse(p_data$sig == 1, p_data$val, NA)
p_data$V3 <- ifelse(p_data$sig == 0, p_data$val, NA)

p_data$label <- ifelse(is.na(p_data$p), paste0(round(p_data$coef,4)), 
       ifelse(p_data$p > .1,  round(p_data$coef*100,4), 
              ifelse(p_data$p > .05,  paste0( round(p_data$coef*100,4), "*"), 
                     ifelse(p_data$p > .01,  paste0( round(p_data$coef*100,4), "**"), paste0( round(p_data$coef*100,4), "***")))))
p_data$label <- ifelse(p_data$x == 5, round(p_data$coef,4),p_data$label)


#coef_mutohyo$topic <- factor(coef_mutohyo$topic, levels=(coef_mutohyo$topic)[order(coef_mutohyo$val)])

p_data <- p_data %>%
  mutate(topic =  factor(topic, levels = order)) %>%
  arrange(topic) 

#ggplot(coef_mutohyo, aes(x = x, y = topic, label = round(coef,4))) + geom_text(color = "black", size = 5)
{
  df <- data.frame(matrix(ncol = ncol(p_data), nrow = length(unique(p_data$topic))))
  colnames(df) <- colnames(p_data)
  df$x <- 4
  df$topic <- unique(p_data$topic)
  df$label <- "YES"
  p_data <- rbind(p_data, df)
  }


ggplot(p_data, aes(x = factor(x), y = topic, label = label)) + 
  annotate(
    geom="tile", x=factor(p_data$x), y=p_data$topic,
    fill = scales::colour_ramp(c("blue4", "white","darkred"))(p_data$V2)
  ) +
  annotate(
    geom="tile", x=factor(p_data$x), y=p_data$topic,
    fill = scales::colour_ramp(c("blue", "white","red"))(p_data$V3)
  ) + geom_text(color = "black", size = 5) + 
  scale_x_discrete(labels = c("無投票ダミー", "当選回数", "選挙区定数", "統制変数", "R^2 (adjusted)", "トピック平均(%)"), position = "top")  + 
  xlab("") + ylab("") +
  ggtitle("") + 
  theme(text = element_text(size=15), panel.background = element_blank()) 


### create the regression tables
coef_mutohyo <- coef_mutohyo %>% arrange(coef) 
order <- c(coef_mutohyo$topic)

screenreg(list(mods[[16]],mods[[24]],mods[[4]],mods[[9]],mods[[21]],
               mods[[13]],mods[[17]],mods[[18]],mods[[5]],mods[[23]],
               mods[[12]],mods[[10]],mods[[20]],mods[[8]],mods[[14]],
               mods[[7]],mods[[15]],mods[[11]],mods[[19]],mods[[6]],
               mods[[22]],mods[[2]],mods[[1]],mods[[3]]),
          custom.note = paste("%stars.",
                              "All models were estimated using OLS regression.  Standard errors are in parenthesis. Unit of analysis ",
                              " is how much legislators had spoken on the specific topics. "),
          caption = "Table of figure 1 full regression results: Effect political competition had on topics covered by legislator",
          head.tag = TRUE,
          custom.model.names=unlist(list(unlist(dvs[16]),unlist(dvs[24]),unlist(dvs[4]),unlist(dvs[9]),unlist(dvs[21]),
                                  unlist(dvs[13]),unlist(dvs[17]),unlist(dvs[18]),unlist(dvs[5]),unlist(dvs[23]),
                                  unlist(dvs[12]),unlist(dvs[10]),unlist(dvs[20]),unlist(dvs[8]),unlist(dvs[14]),
                                  unlist(dvs[7]),unlist(dvs[15]),unlist(dvs[11]),unlist(dvs[19]),unlist(dvs[6]),
                                  unlist(dvs[22]),unlist(dvs[2]),unlist(dvs[1]),unlist(dvs[3]))),
          include.loglik = T, digits = 4)

htmlreg(list(mods[[16]],mods[[24]],mods[[4]],mods[[9]],mods[[21]],
               mods[[13]],mods[[17]],mods[[18]],mods[[5]],mods[[23]],
               mods[[12]],mods[[10]],mods[[20]],mods[[8]],mods[[14]],
               mods[[7]],mods[[15]],mods[[11]],mods[[19]],mods[[6]],
               mods[[22]],mods[[2]],mods[[1]],mods[[3]]),
          custom.note = paste("%stars.",
                              "All models were estimated using OLS regression.  Standard errors are in parenthesis. Unit of analysis ",
                              " is how much legislators had spoken on the specific topics. Results columns are ordered by the magnitude ",
                              "of the uncontested (mutohyo) effect (positive first)."),
          caption = "Table of figure 1 full regression results: Effect political competition had on topics covered by legislator",
          head.tag = TRUE,
          custom.model.names=unlist(list(unlist(dvs[16]),unlist(dvs[24]),unlist(dvs[4]),unlist(dvs[9]),unlist(dvs[21]),
                                         unlist(dvs[13]),unlist(dvs[17]),unlist(dvs[18]),unlist(dvs[5]),unlist(dvs[23]),
                                         unlist(dvs[12]),unlist(dvs[10]),unlist(dvs[20]),unlist(dvs[8]),unlist(dvs[14]),
                                         unlist(dvs[7]),unlist(dvs[15]),unlist(dvs[11]),unlist(dvs[19]),unlist(dvs[6]),
                                         unlist(dvs[22]),unlist(dvs[2]),unlist(dvs[1]),unlist(dvs[3]))),
          include.loglik = T, digits = 4, file = "C:/Users/robba/Dropbox/uncontested_text_analysis_paper/regression_models/poster_fig1_fullreg_sep28.html")


### DiD models

mods <- c()
d_1 <- data.frame()
d_2 <- data.frame()
d_3 <- data.frame()
d_4 <- data.frame()
d_5 <- data.frame()
d_6 <- data.frame()
d_7 <- data.frame()
d_8 <- data.frame()
d_9 <- data.frame()
d_10 <- data.frame()
d_11 <- data.frame()
d_12 <- data.frame()
#a_data$N
#feols(ent ~ gov_ruling + i(time, treat, "2007_2") | id + time, a_data)
#summary(feglm.nb(N ~ gov_ruling + mutohyo | id + time, a_data))
unique(a_data$time)
did_data$time <- paste0(did_data$year, "_", did_data$period)
a_data <- did_data %>% group_by(id) %>% filter(min(period) == 2)
a_data <- a_data %>% group_by(id) %>% filter(max(year) <= min(year) + 8)

for (i in 1:length(dvs)){
  print(i)  
  dv <- dvs[i]

  #if(dv %in% c("length"))
  res <- lapply(dv, function(var) {
    
    #eq <- xpd(..lhs ~ mutohyo+ (1|year + id), ..lhs = var)
    #res <- glmmTMB(eq, data=a_data, ziformula=~1, family=gaussian)
    res <- fixest::feols(xpd(..lhs ~ mutohyo+ gov_ruling | id + time, ..lhs = var), data = a_data)
    return(res)
  })  

  
{
  
  ##save all coefficients into memory
  coef <- res[[1]]$coeftable[[1]][1]
  p_value <- res[[1]]$coeftable[[4]][1]
  #as.data.frame(coef, p_value)
  d <- data.frame(coef = coef, p = p_value, topic = dvs[i])
  d$sig <- int(d$p <.1)
  d$sig <-ifelse(is.na(d$sig), 0 ,d$sig)
  d_1 <- rbind(d,d_1)
  
  ##save all coefficients into memory
  coef <- res[[1]]$coeftable[[1]][2]
  p_value <- res[[1]]$coeftable[[4]][2]
  #as.data.frame(coef, p_value)
  d <- data.frame(coef = coef, p = p_value, topic = dvs[i])
  d$sig <- int(d$p <.1)
  d$sig <-ifelse(is.na(d$sig), 0 ,d$sig)
  d_2 <- rbind(d,d_2)
  
  ##save all coefficients into memory
  coef <- res[[1]]$coeftable[[1]][3]
  p_value <- res[[1]]$coeftable[[4]][3]
  #as.data.frame(coef, p_value)
  d <- data.frame(coef = coef, p = p_value, topic = dvs[i])
  d$sig <- int(d$p <.1)
  d$sig <-ifelse(is.na(d$sig), 0 ,d$sig)
  d_3 <- rbind(d,d_3)
  
  ##save all coefficients into memory
  coef <- res[[1]]$coeftable[[1]][4]
  p_value <- res[[1]]$coeftable[[4]][4]
  #as.data.frame(coef, p_value)
  d <- data.frame(coef = coef, p = p_value, topic = dvs[i])
  d$sig <- int(d$p <.1)
  d$sig <-ifelse(is.na(d$sig), 0 ,d$sig)
  d_4 <- rbind(d,d_4)

  ##save all coefficients into memory
  coef <- res[[1]]$coeftable[[1]][5]
  p_value <- res[[1]]$coeftable[[4]][5]
  #as.data.frame(coef, p_value)
  d <- data.frame(coef = coef, p = p_value, topic = dvs[i])
  d$sig <- int(d$p <.1)
  d$sig <-ifelse(is.na(d$sig), 0 ,d$sig)
  d_5 <- rbind(d,d_5)
  
  ##save all coefficients into memory
  coef <- res[[1]]$coeftable[[1]][6]
  p_value <- res[[1]]$coeftable[[4]][6]
  #as.data.frame(coef, p_value)
  d <- data.frame(coef = coef, p = p_value, topic = dvs[i])
  d$sig <- int(d$p <.1)
  d$sig <-ifelse(is.na(d$sig), 0 ,d$sig)
  d_6 <- rbind(d,d_6)
  
  ##save all coefficients into memory
  coef <- res[[1]]$coeftable[[1]][7]
  p_value <- res[[1]]$coeftable[[4]][7]
  #as.data.frame(coef, p_value)
  d <- data.frame(coef = coef, p = p_value, topic = dvs[i])
  d$sig <- int(d$p <.1)
  d$sig <-ifelse(is.na(d$sig), 0 ,d$sig)
  d_7 <- rbind(d,d_7)
  
  ##save all coefficients into memory
  coef <- res[[1]]$coeftable[[1]][8]
  p_value <- res[[1]]$coeftable[[4]][8]
  #as.data.frame(coef, p_value)
  d <- data.frame(coef = coef, p = p_value, topic = dvs[i])
  d$sig <- int(d$p <.1)
  d$sig <-ifelse(is.na(d$sig), 0 ,d$sig)
  d_8 <- rbind(d,d_8)
  
  ##save all coefficients into memory
  coef <- res[[1]]$coeftable[[1]][9]
  p_value <- res[[1]]$coeftable[[4]][9]
  #as.data.frame(coef, p_value)
  d <- data.frame(coef = coef, p = p_value, topic = dvs[i])
  d$sig <- int(d$p <.1)
  d$sig <-ifelse(is.na(d$sig), 0 ,d$sig)
  d_9 <- rbind(d,d_9)
  
  ##save all coefficients into memory
  coef <- res[[1]]$coeftable[[1]][10]
  p_value <- res[[1]]$coeftable[[4]][10]
  #as.data.frame(coef, p_value)
  d <- data.frame(coef = coef, p = p_value, topic = dvs[i])
  d$sig <- int(d$p <.1)
  d$sig <-ifelse(is.na(d$sig), 0 ,d$sig)
  d_10 <- rbind(d,d_10)
  
  ### get r2 
  coef <- r2(res[[1]])[3]
  p_value <- NA
  d <- data.frame(coef = coef, p = p_value, topic = dvs[i])
  d$sig <- NA
  d_11 <- rbind(d,d_11)
  
  ### get topic mean
  coef <- mean(did_data[,dvs[i]], na.rm = TRUE)
  p_value <- NA
  d <- data.frame(coef = coef, p = p_value, topic = dvs[i])
  d$sig <- NA
  d_12 <- rbind(d,d_12)
  }
  mods[[i]] <- res
}

d_1$x <- 1
d_2$x <- 2
d_3$x <- 3
d_4$x <- 4
d_5$x <- 5
d_6$x <- 6
d_7$x <- 7
d_8$x <- 8
d_9$x <- 9
d_10$x <- 10
d_11$x <- 11
d_12$x <- 12



d_1 <- d_1 %>%mutate( max = max(abs(d_1$coef)),
                                        val = (coef/max)/2+.5  )
d_2 <- d_2 %>%mutate( max = max(abs(d_2$coef)),
                      val = (coef/max)/2+.5  )
d_3 <- d_3 %>%mutate( max = max(abs(d_3$coef)),
                      val = (coef/max)/2+.5  )
d_4 <- d_4 %>%mutate( max = max(abs(d_4$coef)),
                      val = (coef/max)/2+.5  )
d_5 <- d_5 %>%mutate( max = max(abs(d_5$coef)),
                      val = (coef/max)/2+.5  )
d_6 <- d_6 %>%mutate( max = max(abs(d_6$coef)),
                      val = (coef/max)/2+.5  )
d_7 <- d_7 %>%mutate( max = max(abs(d_7$coef)),
                      val = (coef/max)/2+.5  )
d_8 <- d_8 %>%mutate( max = max(abs(d_8$coef)),
                      val = (coef/max)/2+.5  )
d_9 <- d_9 %>%mutate( max = max(abs(d_9$coef)),
                      val = (coef/max)/2+.5  )
d_10 <- d_10 %>%mutate( max = max(abs(d_10$coef)),
                      val = (coef/max)/2+.5  )
d_11 <- d_11 %>%mutate( max = max(abs(d_11$coef)),
                      val = (coef/max)/2+.5  )
d_12 <- d_12 %>%mutate( max = max(abs(d_12$coef)),
                      val = (coef/max)/2+.5  )

p_data <- rbind(d_1,d_2,d_3,d_4,d_5,d_6,d_7,d_8,d_9,d_10,d_11,d_12)

p_data$V2 <- ifelse(p_data$sig == 1, p_data$val, NA)
p_data$V3 <- ifelse(p_data$sig == 0, p_data$val, NA)

p_data$label <- ifelse(is.na(p_data$p), paste0(round(p_data$coef,4)), 
                       ifelse(p_data$p > .1,  round(p_data$coef*100,4), 
                              ifelse(p_data$p > .05,  paste0( round(p_data$coef*100,4), "*"), 
                                     ifelse(p_data$p > .01,  paste0( round(p_data$coef*100,4), "**"), paste0( round(p_data$coef*100,4), "***")))))
p_data$label <- ifelse(p_data$x == 5, round(p_data$coef,4),p_data$label)
p_data <- p_data %>%
  mutate(topic =  factor(topic, levels = order)) %>%
  arrange(topic) 

#ggplot(coef_mutohyo, aes(x = x, y = topic, label = round(coef,4))) + geom_text(color = "black", size = 5)
factor(p_data$x)

ggplot(p_data, aes(x = factor(x), y = topic, label = label)) +
  annotate(
    geom="tile", x=factor(p_data$x), y=p_data$topic,
    fill = scales::colour_ramp(c("blue4", "white","darkred"))(p_data$V2)
  ) +
  annotate(
    geom="tile", x=factor(p_data$x), y=p_data$topic,
    fill = scales::colour_ramp(c("blue", "white","red"))(p_data$V3)
  ) + geom_text(color = "black", size = 4) + geom_vline(xintercept = 5.5, size = 1.5) +
  annotate(x=5.5,y=26.2,label="Uncontested occurs",vjust=2,geom="label")  +  
  scale_x_discrete(labels = c("Gov. party", "2008", "2009","2010","2011_before", "2011_after", 
                              "2012", "2013", "2014", "2015", "R^2 (adjusted)", "Topic mean"), position = "top")  + 
  theme(text = element_text(size=13), panel.background = element_blank()) + xlab("") + ylab("") +
  ggtitle("")


### comment on p.19
## avg. reelect
reg_data %>% group_by(mutohyo) %>% summarise(mean(reelect), mean(totseat))
## 
see_1 <- reg_data %>% group_by(mutohyo) %>% summarise(N_mut = n())
see_2 <- reg_data %>% group_by(totseat) %>% summarise(N_tot = n())

see_3 <- reg_data %>% group_by(mutohyo, totseat) %>% summarise(N = n())
see <- merge(see_1, see_3, by = c("mutohyo"))
see <- merge(see, see_2, by = c("totseat"))
see$N
see %>% mutate(N/N_mut, N/N_tot)


### policy area analysis


reg_data <- reg_data %>% mutate(pork_barrel = agriculture + roads_infrastructure,
                    new_policy_area = forestry_fishing_environment_renewables + women_work_children,
                    healthcare = hospitals_patients + disabilities + cancer_health_patients)

dvs_2 <- c("pork_barrel", "new_policy_area", "healthcare")

###reg modellings
mods <- c()
coef_mutohyo <- data.frame()
d_2 <- data.frame()
d_3 <- data.frame()
d_4 <- data.frame()
d_5 <- data.frame()


for (i in 1:length(dvs_2)){
  dv <- dvs_2[i]
  #if(dv %in% c("length"))
  res  <- lm(get(dv) ~ mutohyo + reelect + totseat + secretary + city_assembly + party + gov_ruling + age + factor(pref_id) + factor(period)
             , data= reg_data)
  ### get mutohyo values
  coef <- res$coefficients[2]
  res_cl <- coeftest(res, vcov = vcovCL, cluster = ~pref_id + period)
  p_value <- res_cl[2,4]
  #as.data.frame(coef, p_value)
  d <- data.frame(coef = coef, p = p_value, topic = dvs_2[i])
  d$sig <- int(d$p <.1)
  d$sig <-ifelse(is.na(d$sig), 0 ,d$sig)
  
  coef_mutohyo <- rbind(d,coef_mutohyo)
  
  ### get reelect values
  summary(res)
  coef <- res$coefficients[3]
  res_cl <- coeftest(res, vcov = vcovCL, cluster = ~pref_id + period)
  p_value <- res_cl[3,4]
  #as.data.frame(coef, p_value)
  d <- data.frame(coef = coef, p = p_value, topic = dvs_2[i])
  d$sig <- int(d$p <.1)
  d$sig <-ifelse(is.na(d$sig), 0 ,d$sig)
  
  d_2 <- rbind(d,d_2)
  
  ### get totseat values
  summary(res)
  coef <- res$coefficients[4]
  res_cl <- coeftest(res, vcov = vcovCL, cluster = ~pref_id + period)
  p_value <- res_cl[4,4]
  #as.data.frame(coef, p_value)
  d <- data.frame(coef = coef, p = p_value, topic = dvs_2[i])
  d$sig <- int(d$p <.1)
  d$sig <-ifelse(is.na(d$sig), 0 ,d$sig)
  
  d_3 <- rbind(d,d_3)
  
  ### get r2 value
  coef <- summary(res)$adj.r.squared
  p_value <- NA
  d <- data.frame(coef = coef, p = p_value, topic = dvs_2[i])
  d$sig <- NA
  d_4 <- rbind(d,d_4)
  
  ### get topic mean
  coef <- mean(reg_data[,dvs[i]], na.rm = TRUE)
  p_value <- NA
  d <- data.frame(coef = coef, p = p_value, topic = dvs_2[i])
  d$sig <- NA
  d_5 <- rbind(d,d_5)
  
  mods[[i]] <- res
}
coef_mutohyo$x <- 1
d_2$x <- 2
d_3$x <- 3
### x = 4 is the one that says "YES"
d_4$x <- 5
d_5$x <- 6

coef_mutohyo <- coef_mutohyo %>% arrange(coef) 
order <- c(coef_mutohyo$topic)

coef_mutohyo <- coef_mutohyo %>%mutate( max = max(abs(coef_mutohyo$coef)),
                                        val = (coef/max)/2+.5  )
d_2 <- d_2 %>%mutate( max = max(abs(d_2$coef)),
                      val = (coef/max)/2+.5  )
d_3 <- d_3 %>%mutate( max = max(abs(d_3$coef)),
                      val = (coef/max)/2+.5  )
d_4 <- d_4 %>%mutate( max = max(abs(d_4$coef)),
                      val = (coef/max)/2+.5  )
d_5 <- d_5 %>%mutate( max = max(abs(d_5$coef)),
                      val = (coef/max)/2+.5  )

p_data <- rbind(coef_mutohyo, d_2,d_3,d_4,d_5)

p_data$V2 <- ifelse(p_data$sig == 1, p_data$val, NA)
p_data$V3 <- ifelse(p_data$sig == 0, p_data$val, NA)

p_data$label <- ifelse(is.na(p_data$p), paste0(round(p_data$coef,4)), 
                       ifelse(p_data$p > .1,  round(p_data$coef*100,4), 
                              ifelse(p_data$p > .05,  paste0( round(p_data$coef*100,4), "*"), 
                                     ifelse(p_data$p > .01,  paste0( round(p_data$coef*100,4), "**"), paste0( round(p_data$coef*100,4), "***")))))
p_data$label <- ifelse(p_data$x == 5, round(p_data$coef,4),p_data$label)

#coef_mutohyo$topic <- factor(coef_mutohyo$topic, levels=(coef_mutohyo$topic)[order(coef_mutohyo$val)])

p_data <- p_data %>%
  mutate(topic =  factor(topic, levels = order)) %>%
  arrange(topic) 

### create a column that says control variables are included
{
df <- data.frame(matrix(ncol = ncol(p_data), nrow = length(unique(p_data$topic))))
colnames(df) <- colnames(p_data)
df$x <- 4
df$topic <- unique(p_data$topic)
df$label <- "YES"
p_data <- rbind(p_data, df)
}

ggplot(p_data, aes(x = factor(x), y = topic, label = label)) + 
  annotate(
    geom="tile", x=factor(p_data$x), y=p_data$topic,
    fill = scales::colour_ramp(c("blue3", "white","darkred"))(p_data$V2)
  ) +
  annotate(
    geom="tile", x=factor(p_data$x), y=p_data$topic,
    fill = scales::colour_ramp(c("blue", "white","red"))(p_data$V3)
  ) + geom_text(color = "black", size = 5) + 
  scale_x_discrete(labels = c("Uncontested", "Tenure", "District magnitude", "Control", "R^2 (adjusted)", "Topic mean (%)"), position = "top")  + 
  scale_y_discrete(labels = c("Healthcare", "New policy areas", "Pork-barrel"))  + 
  xlab("") + ylab("") +
  ggtitle("") + 
  theme(text = element_text(size=15), panel.background = element_blank()) 



htmlreg(list(mods[[1]],mods[[2]],mods[[3]]),  type = "clustered",
        cluster = ~pref_id + period , stars = c(0.01, 0.05, 0.1, 0.15),
        custom.note = paste("%stars.",
                            "Time period was 2003 - 2019 with a total 4 distinct electoral periods included. N = 2682",
                            "in all models"),
        caption = "Table 2: Regular regression models for effect on quantitative measurements",
        head.tag = TRUE,
        custom.model.names=c(dvs_2[1], dvs_2[2], dvs_2[3]),
        include.loglik = T, file = "C:/Users/robba/Dropbox/uncontested_text_analysis_paper/regression_models/fig5_full_regression_sep1.html")








