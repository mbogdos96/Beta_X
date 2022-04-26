#start here
install.packages('tidyverse')
install.packages('sicegar')
install.packages('nlme')
install.packages('minpack.lm')
install.packages('propagate')
install.packages('investr')
install.packages('ggrepel')
install.packages("ragg")
#next here
library(investr)
library(sicegar)
library(minpack.lm)
library(nlme)
library(nls2)
library(tidyverse)
library(propagate)
library(ggrepel)

#next here

tetrakis_data <- read_csv('beta_X_R_import_tetrakis.csv')
kinetics_data <- read_csv('beta_X_R_import_kinetics.csv')
tBu_data <- read_csv('beta_X_R_import_tBu.csv')
kinetics <- kinetics_data %>%
  select(LG,Time_hours,Selectivity) %>%
  rename(Leaving_Group = 'LG', Time_hours = 'Time_hours', Selectivity_Scale = 'Selectivity')
PtBu3 <- tBu_data %>%
  rename(Leaving_Group = 'LG') %>%
  select_if(~ !any(is.na(.))) 
PPh3 <- tetrakis_data %>%
  select_if(~ !any(is.na(.))) %>%
  rename(Sel_scale = 'Selectivity', LGpKa = 'pKa', X_group = 'LG')
TEPdata <- read_csv('TEPOAc.csv')

#PtBu3

attach(PtBu3) 
gompertz.ptbu3=nlsLM(Selectivity~A*exp(-exp(B-C*pKa))+D, start = list(A=-2,B=-20,C=10.5,D=1))
coef(gompertz.ptbu3) 
detach(PtBu3)

#investr method; works fit + CI

new.data <- data.frame(pKa=seq(-5,2,by=0.1))
interval <- as_tibble(predFit(gompertz.ptbu3, newdata = new.data, interval = "confidence",level = 0.67)) %>%
  mutate(pKa = new.data$pKa)
new.dataPtbu3 <- data.frame(pKa=seq(-5,2,by=0.1))
pred.interv <- as_tibble(predFit(gompertz.ptbu3, newdata = new.dataPtbu3, interval = "prediction",level = 0.67)) %>%
  mutate(pKa= new.dataPtbu3$pKa)

#PPh3 time

attach(PPh3)
gompertz.pph3=nlsLM(Sel_scale~A*exp(-exp(B-C*LGpKa))+D, start = list(A=-2,B=7,C=2,D=1))
coef(gompertz.pph3)
detach(PPh3)

new.dataPPh3 <- data.frame(LGpKa=seq(2.5,6,by=0.1))
interval_PPh3 <- as_tibble(predFit(gompertz.pph3,newdata = new.dataPPh3,interval = "confidence", level = 0.67)) %>%
  mutate(pKa = new.dataPPh3$LGpKa)
new.datatetrakis <- data.frame(LGpKa=seq(2.5,6,by=0.1))
pred.intervPPh3 <- as_tibble(predFit(gompertz.pph3, newdata = new.datatetrakis, interval = "prediction", level = 0.67)) %>%
  mutate(pKa= new.datatetrakis$LGpKa)

#plot TEP stuff

TEP_lr <- lm(Selectivity ~ TEP,(data = TEPdata))

PtBu3_ee <- tBu_data %>% mutate(pct_Hyd = case_when(H_X > 0 ~ (100*H_X)/(1+H_X), Selectivity >= 1 ~ 0, Selectivity <= -1 ~ 100), pct_X = case_when(X_H > 0 ~ (100*X_H)/(1+X_H), Selectivity >= 1 ~ 100, Selectivity <= -1 ~ 0), ee_t_Scale = pct_X - pct_Hyd)

gompertz.ptbu3ee <- nlsLM(ee_t_Scale~A*exp(-exp(B-C*pKa))+D, start = list(A=-200,B=-500,C=200,D=100), data = PtBu3_ee)
coef(gompertz.ptbu3ee)

ptbu3ee_new <- data.frame(pKa=seq(-5,2,by=0.1))
interval_ptbuee <- as_tibble(predFit(gompertz.ptbu3ee, newdata = new.data, interval = "confidence",level = 0.67)) %>%
  mutate(pKa = new.data$pKa)
new.dataPtbu3ee <- data.frame(pKa=seq(-5,2,by=0.1))
pred.intervee <- as_tibble(predFit(gompertz.ptbu3ee, newdata = new.dataPtbu3, interval = "prediction",level = 0.67)) %>%
  mutate(pKa= new.dataPtbu3$pKa)

tBu3pl <- nlsLM(Selectivity~A + B/(1+exp(C*pKa)), start = list(A=-1,B=2,C=-1), data = PtBu3_ee)

tBu4pl <- nlsLM(Selectivity~(A - D)/(1 + B*exp(-C*pKa))+D, start = list(A=-1,B=0,C=8,D=1), data = PtBu3_ee)

tBurichards <- nlsLM(Selectivity~A*(1+(B-1)*exp(-C*(pKa-D))^(1/(1-B))+E), start = list(A=-2,B=1.5,C=7,D=-3,E=1), data = PtBu3_ee)

comp_funct_tBu <- ggplot(data=NULL) +
  geom_function(fun = ~predict(tBu3pl,data.frame(pKa=.x)), size = 1, colour = "red") +
  geom_function(fun = ~predict(tBu4pl,data.frame(pKa=.x)), size = 1, colour = "gray") +
  geom_function(fun = ~predict(gompertz.ptbu3,data.frame(pKa = .x)), size = 1, colour = "blue") +
  geom_point(data=PtBu3_ee, mapping = aes(x=pKa, y=Selectivity))

#plots for publication

themething <- theme(axis.line = element_line(size=1), panel.background=element_rect(fill="white"), axis.ticks = element_line(size=1), axis.text = element_text(size = rel(1.5)), axis.title = element_text(size = rel(1.5)),plot.title = element_text(size=rel(1.75),face="bold"),plot.title.position = "panel",legend.key = element_blank(),legend.background = element_blank())

themething_big_labs <- theme(axis.line = element_line(size=1), text = element_text(size=20), panel.background=element_rect(fill="white"), axis.ticks = element_line(size=1), axis.text = element_text(size = rel(1.5)), axis.title = element_text(size = rel(1.5)),plot.title = element_text(size=rel(1.75),face="bold"),plot.title.position = "panel")

sel_PtBu3_imgs_pub <- ggplot(data=NULL) +
  geom_ribbon(data = interval, aes(x= pKa, ymin=lwr, ymax=upr), alpha = 0.7, inherit.aes = FALSE, fill="gray") + 
  geom_ribbon(data = pred.interv, aes(x=pKa,ymin=lwr, ymax=upr), alpha = 0.3, inherit.aes = FALSE, fill = "gray") +
  geom_function(fun = ~predict(gompertz.ptbu3,data.frame(pKa = .x)),size=1,colour="black") +
  geom_point(data = filter(PtBu3, Selectivity > 0), mapping = aes(x=pKa,y=Selectivity), colour = "#A943E2",size=6) +
  geom_point(data = filter(PtBu3, Selectivity < 0), mapping = aes(x=pKa,y=Selectivity), colour = "#377337",size=6) +
  geom_label_repel(data = filter(PtBu3, Selectivity >= 0.96), mapping = aes(x=pKa,y=Selectivity, label = paste0(Leaving_Group)), nudge_y = -0.1, nudge_x = 0, size=7) +
  geom_label_repel(data = filter(PtBu3, Selectivity < 0.96 & Selectivity > -0.5), mapping = aes(x=pKa,y=Selectivity, label = paste0(Leaving_Group)), nudge_y = 0, nudge_x = 3, size=7) +
  geom_label_repel(data = filter(PtBu3, Selectivity < -0.5 & Selectivity > -1), mapping = aes(x=pKa,y=Selectivity, label = paste0(Leaving_Group)), nudge_y = -0.1, nudge_x = -5, size=7) +
  geom_label_repel(data = filter(PtBu3, Selectivity <= -1 & pKa < 15), mapping = aes(x=pKa,y=Selectivity, label = paste0(Leaving_Group)), nudge_y = 0.4, nudge_x = 0, size=7) +
  labs(y=expression(paste(Selectivity[X/H])), x=expression(paste(pK[aH]))) +
  scale_x_continuous(breaks = c(-10,-5,0,5,10,15),limits = c(-10,15)) +
  themething_big_labs

sel_PPh3_pub_img <- ggplot(data=NULL) +
  geom_ribbon(data=interval_PPh3, aes(x=pKa, ymin=lwr, ymax=upr),alpha=0.2) +
  geom_ribbon(data = pred.intervPPh3, aes(x=pKa,ymin=lwr, ymax=upr),alpha=0.08) +
  geom_function(fun = ~predict(gompertz.pph3,data.frame(LGpKa=.x)),size=1,colour="black") +
  geom_point(data = filter(PPh3, Sel_scale > 0), mapping = aes(x=LGpKa,y=Sel_scale), colour = "#A943E2",size=6) +
  geom_point(data = filter(PPh3, Sel_scale < 0), mapping = aes(x=LGpKa,y=Sel_scale), colour = "#377337",size=6) +
  geom_label_repel(data = filter(PPh3, Sel_scale >= 0.96), mapping = aes(x=LGpKa,y=Sel_scale, label = paste0(X_group)), nudge_y = -0.3, nudge_x = -3, size=7) +
  geom_label_repel(data = filter(PPh3, Sel_scale < 0.96 & Sel_scale > -0.5), mapping = aes(x=LGpKa,y=Sel_scale, label = paste0(X_group)), nudge_y = 0, nudge_x = 2, size=7) +
  geom_label_repel(data = filter(PPh3, Sel_scale < -0.5 & Sel_scale >= -1), mapping = aes(x=LGpKa,y=Sel_scale, label = paste0(X_group)), nudge_y = -0.1, nudge_x = -1, size=7) +
  labs(y=expression(paste(Selectivity[X/H])), x=expression(paste(pK[aH]))) +
  scale_x_continuous(breaks = c(-10,-5,0,5,10,15),limits = c(-10,15)) +
  themething_big_labs
  
both_L_labels_pub <- ggplot(data=NULL) +
  geom_ribbon(data = interval, aes(x= pKa, ymin=lwr, ymax=upr), alpha = 0.7, inherit.aes = FALSE, fill="gray") + 
  geom_ribbon(data = pred.interv, aes(x=pKa,ymin=lwr, ymax=upr), alpha = 0.3, inherit.aes = FALSE, fill = "gray") +
  geom_ribbon(data=interval_PPh3, aes(x=pKa, ymin=lwr, ymax=upr),alpha=0.2) +
  geom_ribbon(data = pred.intervPPh3, aes(x=pKa,ymin=lwr, ymax=upr),alpha=0.08) +
  geom_function(fun = ~predict(gompertz.ptbu3,data.frame(pKa = .x)),size=1,colour="black") +
  geom_function(fun = ~predict(gompertz.pph3,data.frame(LGpKa=.x)),size=1,colour="black") +
  geom_point(data = filter(PtBu3, Selectivity > 0), mapping = aes(x=pKa,y=Selectivity), colour = "#A943E2",size=6) +
  geom_point(data = filter(PtBu3, Selectivity < 0), mapping = aes(x=pKa,y=Selectivity), colour = "#377337",size=6) +
  geom_point(data = filter(PPh3, Sel_scale > 0), mapping = aes(x=LGpKa,y=Sel_scale), colour = "#A943E2",size=6) +
  geom_point(data = filter(PPh3, Sel_scale < 0), mapping = aes(x=LGpKa,y=Sel_scale), colour = "#377337",size=6) +
  labs(y=expression(paste(Selectivity[X/H])), x=expression(paste(pK[aH]))) +
  scale_x_continuous(breaks = c(-10,-5,0,5,10,15),limits = c(-10,15)) +
  themething_big_labs
  
TEP_lab_pub_alt <- ggplot(data = NULL) +
  geom_point(data = filter(TEPdata, Selectivity > 0), mapping = aes(x=TEP,y=Selectivity), colour = "#A943E2",size=10) +
  geom_point(data = filter(TEPdata, Selectivity < 0), mapping = aes(x=TEP,y=Selectivity), colour = "#377337",size=10) +
  geom_label_repel(data = filter(TEPdata, TEP < 2067), mapping = aes(x=TEP,y=Selectivity, label = paste0(Substituent)), nudge_y = -0.375, nudge_x = 0, size=12) +
  geom_label_repel(data = filter(TEPdata, TEP > 2067), mapping = aes(x=TEP,y=Selectivity, label = paste0(Substituent)), nudge_y = 0.1, nudge_x = -1.5, size=12) +
  labs(y=expression(paste(Selectivity[X/H])), x=expression(paste(TEP (cm^-1)))) +
  scale_y_continuous(breaks = c(-1,-0.5,0,0.5,1),limits = c(-1.2,1.2)) +
  scale_x_continuous(breaks = c(2055,2060,2065,2070,2075),limits = c(2054,2076)) +
  themething_big_labs

kinetics_plot_pub <- ggplot(data = kinetics_data, aes(x = Time_hours, y=Selectivity, group = LG, colour = LG)) +
  geom_point(mapping = aes(x=Time_hours, y=Selectivity), data = . %>% filter(LG %in% c("OMs")), size=4) +
  geom_point(mapping = aes(x=Time_hours, y=Selectivity), data = . %>% filter(LG %in% c("OTFA")), size=4) +
  geom_point(mapping = aes(x=Time_hours, y=Selectivity), data = . %>% filter(LG %in% c("OAc")), size=4) +
  geom_point(mapping = aes(x=Time_hours, y=Selectivity), data = . %>% filter(LG %in% c("OPO(OPh)2")), size=4) + 
  labs(y=expression(paste(Selectivity[H/X])), x='Time (hours)') +
  themething

b_acyl_pub <- ggplot(data = kinetics_data) +
  geom_point(mapping = aes(x = Time_hours, y = B_acyl_pct), data = . %>% filter(LG %in% c("OTFA")), size = 5) +
  labs(y=expression(paste(beta,"-Acyl product (%)")),x='Time (hours)') +
  themething

#E/Z ratios and barriers

ez_data <- read.csv("E_Z_ratios.csv", sep = "")

EZ_ratio_plot <- ggplot(data = ez_data, aes(x = time_min, y=Ezratio, colour = X)) +
  geom_point(mapping = aes(x=time_min, y=Ezratio), data = . %>% filter(X %in% c("OMs")), size=4) +
  geom_point(mapping = aes(x=time_min, y=Ezratio), data = . %>% filter(X %in% c("OTFA")), size=4) +
  geom_point(mapping = aes(x=time_min, y=Ezratio), data = . %>% filter(X %in% c("OAc")), size=4) +
  labs(y="E/Z Ratio", x='Time (min)') +
  themething

ddg_pka <- function(pKa) {
  ddg <- (-8.3145 * 298 * log(-(-1.9693456 + (0.9946751 - 1)*exp(exp(-3.8275428 -1.5370155*pKa)))/(-1.9693456 + (0.9946751 + 1)*exp(exp(-3.8275428 -1.5370155*pKa)))))/(4184)
  return(ddg)
}

pka_for_ddg_test <- seq(-5,5,by=0.1)

ddg_gen <- data.frame(pka_for_ddg_test) %>%
  rename(pKa = "pka_for_ddg_test") %>%
  rowwise() %>%
  mutate(ddg = ddg_pka(pKa))

plot_ddg_pka <- ggplot(data = ddg_gen, mapping = aes(x=pKa,y=ddg)) +
  geom_point() + 
  labs(y="(kcal/mol)", x=expression(paste(pK[aH]))) +
  themething_big_labs

ddg_on_ee <- PtBu3_ee %>%
  mutate(ddg_pred = ddg_pka(pKa))

plot_ddg_ee <- ggplot(data = ddg_on_ee, mapping = aes(y=ee_Scale,x=ddg_pred)) +
  geom_point() + 
  themething_big_labs

#in order to reproduce any of the graphs, simply type the name of any of them into the console e.g. kinetics_plot_pub
