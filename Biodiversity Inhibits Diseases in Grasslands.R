
setwd("C:\\Users\\Dell\\Desktop\\zp\\Others\\Jianghongying")

library(tidyverse)
library(readxl)
library(lme4)
library(lmerTest)
library(ggpubr)
library(viridis)
library(ggsci)
library(performance)
library(metafor)
library(eoffice)
library(vegan)
library(glmmTMB)
library(piecewiseSEM)
library(tidybayes)

mytheme <- theme_test()+theme(axis.title.x = element_text(size = 24, face = "bold", colour = "black"),
                              axis.title.y = element_text(size = 24, face = "bold", colour = "black"),
                              axis.text.x = element_text(size = 20, colour = "black",hjust = 0.5),
                              axis.text.y = element_text(size = 20, colour = "black"))+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth = 1, linetype="solid"))+
  #theme(panel.grid = element_blank())+
  theme(plot.margin = unit(x = c(0.3,0.3,0.2,0.2),units = "cm"))+
  theme(legend.text = element_text(size = 14, colour = "black"))+
  theme(plot.title = element_text(size = 25,          #字体大小
                                  hjust = 0.5,          #字体左右的位置
                                  angle = 0))

trait.flex.anova <-
  function(formula, specif.avg, const.avg, ...) 
  {
    # Formula into string form
    form.string<-deparse(substitute(formula))
    
    # Check formula parameter
    form.parts<-unlist(strsplit(form.string,"~"))
    if(length(form.parts) != 2)
      stop("First parameter must be valid one-sided formula, like ~A*B");
    if(nchar(form.parts[1])>0)
      warning("Left side of the formula was ignored!");
    
    # The two average variables into string form
    spec.av <- deparse(substitute(specif.avg))
    cons.av <- deparse(substitute(const.avg))
    
    test.has.onelevel<-function(aov.summ)
    {
      (length(aov.summ) == 1);
    }
    
    test.has.resid<-function(aov.one)
    {
      if(class(aov.one)[1] != "anova")
        warning("specified object is not aov result!");
      nrows <- dim(aov.one)[1]
      if(nrows == 0)
        return( FALSE);
      last.row.lbl <- dimnames(aov.one)[[1]][nrows]
      if(unlist(strsplit(last.row.lbl, " "))[1] != "Residuals")
        return( FALSE);
      if(dim(aov.one)[2] < 5)  # P and F are missing
        return( FALSE);
      TRUE;
    }
    
    # Specific averages ANOVA
    form.1 <- as.formula(
      paste(spec.av,form.parts[2],sep="~"))
    res.1 <- summary( aov(form.1,...))
    if(test.has.onelevel( res.1) == FALSE)
      stop("Cannot evaluate ANOVAs with multiple error levels!")
    res.1 <- res.1[[1]]
    if(test.has.resid( res.1) == FALSE)
      stop("No residual DFs left, cannot continue");
    
    # Constant averages ANOVA
    form.2 <- as.formula(
      paste(cons.av,form.parts[2],sep="~"))
    # no need to test for multilevels by now
    res.2 <- summary( aov(form.2,...))[[1]]
    if(test.has.resid( res.2) == FALSE)
      stop("No residual DFs left in constant averages ANOVA, cannot continue");
    
    
    # Now the differences:
    spec.const.diff <- paste("I(", spec.av, "-", cons.av, ")", sep="")
    
    form.3 <- as.formula(
      paste(spec.const.diff,form.parts[2],sep="~"))
    # no need to test for multilevels by now
    res.3 <- summary( aov(form.3,...))[[1]]
    if(test.has.resid( res.3) == FALSE)
      stop("No residual DFs left in (specif-const) ANOVA, cannot continue");
    
    
    if((dim(res.1) != dim(res.2)) || (dim(res.1) != dim(res.3)))
      stop("Tables from the three ANOVAs have incompatible sizes");
    
    # Create sum of squares table: add SS(Tot) except for null models    
    nrows <- dim(res.1)[1]
    ss.turn <- res.2[,2]
    ss.var  <- res.3[,2]
    ss.tot  <- res.1[,2]
    ss.covar<- ss.tot - ss.turn - ss.var
    ss.row.names <- dimnames(res.1)[[1]]
    if(nrows > 1)
    {
      ss.turn <- c(ss.turn, sum(ss.turn))
      ss.var  <- c(ss.var,  sum(ss.var))
      ss.tot  <- c(ss.tot,  sum(ss.tot))
      ss.covar<- c(ss.covar,sum(ss.covar))
      ss.row.names <- c(ss.row.names, "Total")
      nrows   <- nrows + 1
    }
    else
    {
      # replace row title
      ss.row.names[1] <- "Total"
    }
    SS.tab <- data.frame( Turnover = ss.turn, Intraspec. = ss.var,
                          Covariation = ss.covar, Total = ss.tot,
                          row.names = ss.row.names)
    # Calculate relative fractions
    TotalSS <- SS.tab[nrows, 4] # lower right corner
    SS.tab.rel <- SS.tab / TotalSS
    
    # Collect significances
    if(nrows > 1)  # get rid of the "Total" label again
      ss.row.names <- ss.row.names[-nrows]
    P.tab <- data.frame( Turnover = res.2[,5], Intraspec. = res.3[,5],
                         Total = res.1[,5], row.names = ss.row.names)
    
    res <- list( SumSq=SS.tab, RelSumSq=SS.tab.rel, Pvals=P.tab, 
                 anova.turnover=res.2, anova.total=res.1, anova.diff=res.3)
    class(res)<- "trait.flex"
    res
  }
##  Read data   #### 
env <- read.csv('Supplementary Data 1.csv')
env

ggplot(env)+
  geom_point(aes(x = MAT, y = MAP),shape = 1,size = 5,alpha = 2)+
  mytheme+theme(legend.position = 'none')+
  scale_x_continuous(breaks = seq(-7,9,3))+
  labs(x = "MAT",y = "MAP")

data_all <- read.csv('Supplementary Data 2.csv')
data_all <- data_all%>% mutate(FishZ_ITV_PL=0.5*log((1+ITV_PL+0.0001)/(1-ITV_PL-0.0001)),  ## Fisher Z transformation 
                               FishZ_PL=0.5*log((1+PL+0.0001)/(1-PL-0.0001)),
                               FishZ_Turn=0.5*log((1+Proness+0.0001)/(1-Proness-0.0001)))

data_all_tra <- data_all %>% 
  gather(key = 'Type',value = 'PL',-c(Site,Plot,SR))


####  Linear mixed-effects models   ####
model <- lmer(FishZ_PL~SR+(1|Site),data_all)
model2 <- lmer(FishZ_Turn~SR+(1|Site),data_all)
model3 <- lmer(FishZ_ITV_PL~SR+(1|Site),data_all)
anova(model)
anova(model2)
anova(model3)
fm1<- model.matrix(model) %*% fixef(model)
var(fm1)
fm2<- model.matrix(model2) %*% fixef(model2)
var(fm2)
fm3<- model.matrix(model3) %*% fixef(model3)
var(fm3)
rvd1<-as.data.frame(VarCorr(model))
Sum_var = var(fm1)+rvd1$vcov[1]+ rvd1$vcov[2]

expla_turn=var(fm2)/Sum_var
expla_itv=var(fm3)/Sum_var

## calculate slope and explanation for each general linear model each site ####
slope <- matrix(ncol = 6,nrow = 63)
i=1
for (i in 1:63) {
  data_lm <- filter(data_all,Site == env$Site[i])
  mod1 <- lm(FishZ_PL~SR,data_lm)
  mod2 <- lm(FishZ_Turn~SR,data_lm)
  mod3 <- lm(FishZ_ITV_PL~SR,data_lm)
  anova_decomp <- trait.flex.anova(~SR, FishZ_PL, FishZ_Turn, data=data_lm)
  
  slope[i,1] <- env$Site[i]
  slope[i,2] <- mod1$coefficients[2]
  slope[i,3] <- mod2$coefficients[2]
  slope[i,4] <- mod3$coefficients[2]
  slope[i,5] <- anova_decomp$RelSumSq[1,1:2][1]$Turnover
  slope[i,6] <- anova_decomp$RelSumSq[1,1:2][2]$Intraspec.
  
}
slope <- slope %>% as.data.frame()
colnames(slope)[1] <- "Slope"
slope[is.na(slope)] <- 0
colnames(slope) <- c('Slope','PL','Proness','ITV','ex_Turnover','ex_ITV')

slope


#### Linear mixed-effects model and General linear model for Fig2  ####
Site_data <- env%>% 
  gather(key = 'Type',value = 'Slope',-c(Site,MAT,MAP,LON,LAT,ALT,SWC,SWC_2mm,
                                         pH,Potential,STC,STN,STP,Soil_AP,NH4N,NO3N,
                                         TIC,ex_Turnover,ex_ITV,rate,NestedT))


model1 <- lm(Slope ~ LON,data = Site_data %>% filter(Type=='PL')) ## Not Significant
summary(model1)
a1 <- anova(model1)
model2 <- lm(Slope ~ LAT,data = Site_data %>% filter(Type=='PL')) ## Not Significant
summary(model2)
a2 <- anova(model2)
model3 <- lm(Slope ~ ALT,data = Site_data %>% filter(Type=='PL')) ## Marginally Significant
summary(model3)
a3 <- anova(model3)
model4 <- lm(Slope ~ MAT,data = Site_data %>% filter(Type=='PL')) ## Significant
summary(model4)
a4 <- anova(model4)
model5 <- lm(Slope ~ MAP,data = Site_data %>% filter(Type=='PL')) ## Not Significant
summary(model5)
a5 <- anova(model5)
model6 <- lm(Slope ~ SWC,data = Site_data %>% filter(Type=='PL')) ## Not Significant
summary(model6)
a6 <- anova(model6)
model7 <- lm(Slope ~ pH,data = Site_data %>% filter(Type=='PL')) ## Marginally Significant
summary(model7)
a7 <- anova(model7)
model8 <- lm(Slope ~ STC,data = Site_data %>% filter(Type=='PL')) ## Significant
summary(model8)
a8 <- anova(model8)
model9 <- lm(Slope ~ STN,data = Site_data %>% filter(Type=='PL')) ## Not Significant
summary(model9)
a9 <- anova(model9)
model10 <- lm(Slope ~ STP,data = Site_data %>% filter(Type=='PL')) ## Not Significant
summary(model10)
a10 <- anova(model10)
model11 <- lm(Slope ~ Soil_AP,data = Site_data %>% filter(Type=='PL')) ## Not Significant
summary(model11)
a11 <- anova(model11)
model12 <- lm(Slope ~ NH4N,data = Site_data %>% filter(Type=='PL')) ## Significant
summary(model12)
a12 <- anova(model12)
model13 <- lm(Slope ~ NO3N,data = Site_data %>% filter(Type=='PL')) ## Not Significant
summary(model13)
a13 <- anova(model13)
model14 <- lm(Slope ~ TIC,data = Site_data %>% filter(Type=='PL')) ## Not Significant
summary(model14)
a14 <- anova(model14)


model1 <- lm(Slope ~ LON,data = Site_data %>% filter(Type=='Proness')) ## Not Significant
summary(model1)
a1 <- anova(model1)
model2 <- lm(Slope ~ LAT,data = Site_data %>% filter(Type=='Proness')) ##  Not Significant
summary(model2)
a2 <- anova(model2)
model3 <- lm(Slope ~ ALT,data = Site_data %>% filter(Type=='Proness')) ## Marginally Significant
summary(model3)
a3 <- anova(model3)
model4 <- lm(Slope ~ MAT,data = Site_data %>% filter(Type=='Proness')) ## Significant
summary(model4)
a4 <- anova(model4)
model5 <- lm(Slope ~ MAP,data = Site_data %>% filter(Type=='Proness')) ## Marginally Significant
summary(model5)
a5 <- anova(model5)
model6 <- lm(Slope ~ SWC,data = Site_data %>% filter(Type=='Proness')) ## Not Significant
summary(model6)
a6 <- anova(model6)
model7 <- lm(Slope ~ pH,data = Site_data %>% filter(Type=='Proness')) ## Not Significant
summary(model7)
a7 <- anova(model7)
model8 <- lm(Slope ~ STC,data = Site_data %>% filter(Type=='Proness')) ## Significant
summary(model8)
a8 <- anova(model8)
model9 <- lm(Slope ~ STN,data = Site_data %>% filter(Type=='Proness')) ## Not Significant
summary(model9)
a9 <- anova(model9)
model10 <- lm(Slope ~ STP,data = Site_data %>% filter(Type=='Proness')) ## Not Significant
summary(model10)
a10 <- anova(model10)
model11 <- lm(Slope ~ Soil_AP,data = Site_data %>% filter(Type=='Proness')) ## Not Significant
summary(model11)
a11 <- anova(model11)
model12 <- lm(Slope ~ NH4N,data = Site_data %>% filter(Type=='Proness')) ## Marginally Significant
summary(model12)
a12 <- anova(model12)
model13 <- lm(Slope ~ NO3N,data = Site_data %>% filter(Type=='Proness')) ## Not Significant
summary(model13)
a13 <- anova(model13)
model14 <- lm(Slope ~ TIC,data = Site_data %>% filter(Type=='Proness')) ## Not Significant
summary(model14)
a14 <- anova(model14)


model1 <- lm(Slope ~ LON,data = Site_data %>% filter(Type=='ITV_PL')) ## Not Significant
summary(model1)
a1 <- anova(model1)
model2 <- lm(Slope ~ LAT,data = Site_data %>% filter(Type=='ITV_PL')) ## Not Significant
summary(model2)
a2 <- anova(model2)
model3 <- lm(Slope ~ ALT,data = Site_data %>% filter(Type=='ITV_PL')) ## Marginally Significant
summary(model3)
a3 <- anova(model3)
model4 <- lm(Slope ~ MAT,data = Site_data %>% filter(Type=='ITV_PL')) ## Significant
summary(model4)
a4 <- anova(model4)
model5 <- lm(Slope ~ MAP,data = Site_data %>% filter(Type=='ITV_PL')) ## Not Significant
summary(model5)
a5 <- anova(model5)
model6 <- lm(Slope ~ SWC,data = Site_data %>% filter(Type=='ITV_PL')) ## Not Significant
summary(model6)
a6 <- anova(model6)
model7 <- lm(Slope ~ pH,data = Site_data %>% filter(Type=='ITV_PL')) ## Significant
summary(model7)
a7 <- anova(model7)
model8 <- lm(Slope ~ STC,data = Site_data %>% filter(Type=='ITV_PL')) ## Not Significant
summary(model8)
a8 <- anova(model8)
model9 <- lm(Slope ~ STN,data = Site_data %>% filter(Type=='ITV_PL')) ## Not Significant
summary(model9)
a9 <- anova(model9)
model10 <- lm(Slope ~ STP,data = Site_data %>% filter(Type=='ITV_PL')) ## Not Significant
summary(model10)
a10 <- anova(model10)
model11 <- lm(Slope ~ Soil_AP,data = Site_data %>% filter(Type=='ITV_PL')) ## Marginally Significant
summary(model11)
a11 <- anova(model11)
model12 <- lm(Slope ~ NH4N,data = Site_data %>% filter(Type=='ITV_PL')) ## Significant
summary(model12)
a12 <- anova(model12)
model13 <- lm(Slope ~ NO3N,data = Site_data %>% filter(Type=='ITV_PL')) ## Not Significant
summary(model13)
a13 <- anova(model13)
model14 <- lm(Slope ~ TIC,data = Site_data %>% filter(Type=='ITV_PL')) ## Not Significant
summary(model14)
a14 <- anova(model14)


model1 <- lm(log(rate) ~ LON,data = Site_data) ## Significant
summary(model1)
a1 <- anova(model1)
model2 <- lm(log(rate) ~ LAT,data = Site_data) ## Significant
summary(model2)
a2 <- anova(model2)
model3 <- lm(log(rate) ~ ALT,data = Site_data) ## Significant
summary(model3)
a3 <- anova(model3)
model4 <- lm(log(rate) ~ MAT,data = Site_data) ## Significant
summary(model4)
a4 <- anova(model4)
model5 <- lm(log(rate) ~ MAP,data = Site_data) ## Not Significant
summary(model5)
a5 <- anova(model5)
model6 <- lm(log(rate) ~ SWC,data = Site_data) ## Not Significant
summary(model6)
a6 <- anova(model6)
model7 <- lm(log(rate) ~ pH,data = Site_data) ## Not Significant
summary(model7)
a7 <- anova(model7)
model8 <- lm(log(rate) ~ STC,data = Site_data) ## Significant
summary(model8)
a8 <- anova(model8)
model9 <- lm(log(rate) ~ STN,data = Site_data) ## Not Significant
summary(model9)
a9 <- anova(model9)
model10 <- lm(log(rate) ~ STP,data = Site_data) ## Not Significant
summary(model10)
a10 <- anova(model10)
model11 <- lm(log(rate) ~ Soil_AP,data = Site_data) ## Significant
summary(model11)
a11 <- anova(model11)
model12 <- lm(log(rate) ~ NH4N,data = Site_data) ## Not Significant
summary(model12)
a12 <- anova(model12)
model13 <- lm(log(rate) ~ NO3N,data = Site_data) ## Significant
summary(model13)
a13 <- anova(model13)
model14 <- lm(log(rate) ~ TIC,data = Site_data) ## Not Significant
summary(model14)
a14 <- anova(model14)

#### Linear mixed-effects model and General linear model for Fig 3 #### 

model1 <- lm(NestedT ~ LON,data = env) ## Significant
summary(model1)
a1 <- anova(model1)
model2 <- lm(NestedT ~ LAT,data = env) ##  Significant
summary(model2)
a2 <- anova(model2)
model3 <- lm(NestedT ~ ALT,data = env) ##  Significant
summary(model3)
a3 <- anova(model3)
model4 <- lm(NestedT ~ MAT,data = env) ## Significant
summary(model4)
a4 <- anova(model4)
model5 <- lm(NestedT ~ MAP,data = env) ## Not Significant
summary(model5)
a5 <- anova(model5)
model6 <- lm(NestedT ~ SWC,data = env) ## Not Significant
summary(model6)
a6 <- anova(model6)
model7 <- lm(NestedT ~ pH,data = env) ## Not Significant
summary(model7)
a7 <- anova(model7)
model8 <- lm(NestedT ~ STC,data = env) ## Not Significant
summary(model8)
a8 <- anova(model8)
model9 <- lm(NestedT ~ STN,data = env) ## Significant
summary(model9)
a9 <- anova(model9)
model10 <- lm(NestedT ~ STP,data = env) ## Not Significant
summary(model10)
a10 <- anova(model10)
model11 <- lm(NestedT ~ Soil_AP,data = env) ## Not Significant
summary(model11)
a11 <- anova(model11)
model12 <- lm(NestedT ~ NH4N,data = env) ## Not Significant
summary(model12)
a12 <- anova(model12)
model13 <- lm(NestedT ~ NO3N,data = env) ## Not Significant
summary(model13)
a13 <- anova(model13)
model14 <- lm(NestedT ~ TIC,data = env) ## Not Significant
summary(model14)
a14 <- anova(model14)



#### SEM for Fig 4 a  ####
sem1 <- env %>% dplyr::select(MAT,MAP,Site) %>% left_join(data_all) %>% na.omit() %>% mutate(SR_MAT=SR*MAT,SR_MAP=SR*MAP)

model <- psem(lmer(FishZ_Turn~SR+MAT+MAP+SR:MAP+SR:MAT+(1|Site),sem1),
              lmer(FishZ_ITV_PL~SR+MAT+MAP+SR:MAP+SR:MAT+(1|Site),sem1),
              FishZ_Turn%~~%FishZ_ITV_PL)
summary(model)

model <- psem(lmer(FishZ_Turn~SR+MAT+MAP+SR_MAP+SR_MAT+(1|Site),sem1),
              lmer(FishZ_ITV_PL~SR+MAT+MAP+SR_MAP+SR_MAT+(1|Site),sem1),
              FishZ_Turn%~~%FishZ_ITV_PL)
summary(model)

####   Partial regression plot for Fig 4 b  ####
## Species turnover ~ SR:MAT
yresid <- resid(lmer(FishZ_Turn~SR+MAT+MAP+SR_MAP+(1|Site),sem1))
xresid <- resid(lmer(SR_MAT~SR+MAP+SR_MAP+MAT+(1|Site),sem1))
lm(yresid~xresid) %>% summary()
p1 <- data.frame(yresid,xresid) %>% ggplot(aes(xresid,yresid))+
  geom_point(size=3,shape=1,color="#F39B7FFF")+
  geom_smooth(method = 'lm',color='#0225F8',linewidth=0.5)+
  mytheme+
  labs(x='',y='Species turnover effect | Others')

## Species turnover ~ SR:MAP
yresid <- resid(lmer(FishZ_Turn~SR+MAT+MAP+SR_MAT+(1|Site),sem1))
xresid <- resid(lmer(SR_MAP~SR+MAP+SR_MAT+MAT+(1|Site),sem1))
lm(yresid~xresid) %>% anova()
p2 <- data.frame(yresid,xresid) %>% ggplot(aes(xresid,yresid))+
  geom_point(size=3,shape=1,color="#F39B7FFF")+
  geom_smooth(method = 'lm',color='#0CA133',linewidth=0.5)+
  mytheme+
  labs(x='',y='')

yresid <- resid(lmer(FishZ_ITV_PL~SR+MAT+MAP+SR_MAP+(1|Site),sem1))
xresid <- resid(lmer(SR_MAT~SR+MAP+SR_MAP+MAT+(1|Site),sem1))
p3 <- data.frame(yresid,xresid) %>% ggplot(aes(xresid,yresid))+
  geom_point(size=3,shape=1,color="#8491B4FF")+
  geom_smooth(method = 'lm',color='#0225F8',linewidth=0.5)+
  mytheme+
  labs(x='SR:MAT | Others',y='Intraspecific variation | Others')

yresid <- resid(lmer(FishZ_ITV_PL~SR+MAT+MAP+SR_MAT+(1|Site),sem1))
xresid <- resid(lmer(SR_MAP~SR+MAP+MAT+SR_MAT+(1|Site),sem1))
p4 <- data.frame(yresid,xresid) %>% ggplot(aes(xresid,yresid))+
  geom_point(size=3,shape=1,color="#8491B4FF")+
  #geom_smooth(method = 'lm',color='#0225F8')+
  mytheme+
  labs(x='SR:MAP | Others',y='')

p0 <- ggarrange(p1,p2,p3,p4,nrow = 2,ncol = 2,labels = c('a','b','c','d'),align = 'hv',
                label.x = 0.23,label.y = 0.95,font.label = list(size = 12, face = "bold"))
p0
ggsave('SEM_partial.tiff',height = 6,width = 6)

