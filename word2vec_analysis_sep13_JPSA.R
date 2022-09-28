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


### analysis 3 - word2vec_tf-idf model ###

reg_data <- read_csv("word2vec_data_period_aug31.csv")
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
  
}

colnames(reg_data)
dvs <-  c("global_warming","renewable_energy","gender_equality","lgbt_samesex")




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
  res  <- lm(get(dv) ~ mutohyo + reelect + totseat + secretary + city_assembly + party_2 + gov_ruling + age + factor(pref_id) + factor(period)
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
  scale_x_discrete(labels = c("無投票ダミー", "当選回数", "選挙区定数", "統制変数", "R^2 (adjusted)", "トピック平均 (%)"), position = "top")  + 
  scale_y_discrete(labels = c( "男女平等","地球温暖化","LGBT,性的マイノリティー", "再生可能エネルギー"))  + 
  xlab("") + ylab("") +
  ggtitle("") + 
  theme(text = element_text(size=15), panel.background = element_blank()) 


#### create the regression table for models to have at poster session
### create the regression tables
coef_mutohyo <- coef_mutohyo %>% arrange(coef) 
order <- c(coef_mutohyo$topic)

screenreg(list(mods[[2]],mods[[4]],mods[[1]],mods[[3]]),
          custom.note = paste("%stars.",
                              "All models were estimated using OLS regression.  Standard errors are in parenthesis. Unit of analysis ",
                              " is how much legislators had spoken on the specific topics. "),
          caption = "Table of figure 2 full regression results: Effect political competition had on new policy areas",
          head.tag = TRUE,
          custom.model.names=unlist(list(unlist(dvs[2]),unlist(dvs[4]),unlist(dvs[1]),unlist(dvs[3]))),
          include.loglik = T, digits = 4)

htmlreg(list(mods[[2]],mods[[4]],mods[[1]],mods[[3]]),
        custom.note = paste("%stars.",
                            "All models were estimated using OLS regression.  Standard errors are in parenthesis. Unit of analysis ",
                            " is how much legislators had spoken on the specific topics. "),
        caption = "Table of figure 2 full regression results: Effect political competition had on new policy areas",
        head.tag = TRUE,
        custom.model.names=unlist(list(unlist(dvs[2]),unlist(dvs[4]),unlist(dvs[1]),unlist(dvs[3]))),
        include.loglik = T, digits = 4, file = "C:/Users/robba/Dropbox/uncontested_text_analysis_paper/regression_models/poster_fig2_fullreg_sep28.html")



