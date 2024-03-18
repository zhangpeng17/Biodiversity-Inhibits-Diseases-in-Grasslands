
setwd("C:\\Users\\Dell\\Desktop\\zp\\Zhangpeng")

library(tidyverse)
library(readxl)
library(lme4)
library(lmerTest)
library(nlme)
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

mytheme <- theme_test()+theme(axis.title.x = element_text(size = 12, face = "bold", colour = "black"),
                              axis.title.y = element_text(size = 12, face = "bold", colour = "black"),
                              axis.text.x = element_text(size = 10, colour = "black",hjust = 0.5),
                              axis.text.y = element_text(size = 10, colour = "black"))+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth = 1, linetype="solid"))+
  #theme(panel.grid = element_blank())+
  theme(plot.margin = unit(x = c(0.3,0.3,0.2,0.2),units = "cm"))+
  theme(legend.text = element_text(size = 14, colour = "black"))+
  theme(plot.title = element_text(size = 25, hjust = 0.5, angle = 0))

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


####  Linear mixed-effects models  for fig 1b ####
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

expla_turn/(expla_turn+expla_itv)
expla_itv/(expla_turn+expla_itv)

## calculate slope and explanation for each general linear model each site ####
slope <- matrix(ncol = 9,nrow = 63)
i=1
for (i in 1:63) {
  data_lm <- filter(data_all,Site == env$Site[i])
  mod1 <- lm(FishZ_PL~SR,data_lm)%>% summary()
  mod2 <- lm(FishZ_Turn~SR,data_lm)%>% summary()
  mod3 <- lm(FishZ_ITV_PL~SR,data_lm)%>% summary()
  anova_decomp <- trait.flex.anova(~SR, FishZ_PL, FishZ_Turn, data=data_lm)
  
  slope[i,1] <- env$Site[i]
  slope[i,2] <- mod1$coefficients[2,1]
  slope[i,3] <- mod2$coefficients[2,1]
  slope[i,4] <- mod3$coefficients[2,1]
  slope[i,5] <- anova_decomp$RelSumSq[1,1:2][1]$Turnover
  slope[i,6] <- anova_decomp$RelSumSq[1,1:2][2]$Intraspec.
  slope[i,7] <- mod1$coefficients[2,4]
  slope[i,8] <- mod2$coefficients[2,4]
  slope[i,9] <- mod3$coefficients[2,4]
}
slope <- slope %>% as.data.frame()
colnames(slope)[1] <- "Slope"
slope[is.na(slope)] <- 0
colnames(slope) <- c('Slope','PL','STP','ITV','ex_Turnover','ex_ITV','P-PL','P-STP','P-ITV')

#slope %>% write.csv('Slpoe.csv')


#### Linear mixed-effects model and General linear model for fig 2  ####
Site_data <- env %>% 
  gather(key = 'Type',value = 'Slope',-c(Site,MAT,MAP,LON,LAT,ALT,ex_Turnover,ex_ITV,rate,NestedT,P.PL,P.STP,P.ITV,GrasslandType))

model <- lm(Slope~NestedT,data = Site_data)
summary(model)

model <- lm(Slope~GrasslandType,data = Site_data %>% filter(Type=='Slope_PL') )
summary(model);anova(model)
model <- lm(Slope~GrasslandType,data = Site_data %>% filter(Type=='Slope_STP') )
summary(model);anova(model)
model <- lm(Slope~GrasslandType,data = Site_data %>% filter(Type=='Slope_ITV') )
summary(model);anova(model)

model1 <- lm(Slope ~ LON,data = Site_data %>% filter(Type=='Slope_PL')) ## Not Significant
summary(model1)
anova(model1)
model2 <- lm(Slope ~ LAT,data = Site_data %>% filter(Type=='Slope_PL')) ## Not Significant
summary(model2)
anova(model2)
model3 <- lm(Slope ~ ALT,data = Site_data %>% filter(Type=='Slope_PL')) ## Marginally Significant
summary(model3)
anova(model3)
model4 <- lm(Slope ~ MAT,data = Site_data %>% filter(Type=='Slope_PL')) ## Significant
summary(model4)
anova(model4)
model5 <- lm(Slope ~ MAP,data = Site_data %>% filter(Type=='Slope_PL')) ## Not Significant
summary(model5)
anova(model5)



model1 <- lm(Slope ~ LON,data = Site_data %>% filter(Type=='Slope_STP')) ## Not Significant
summary(model1)
anova(model1)
model2 <- lm(Slope ~ LAT,data = Site_data %>% filter(Type=='Slope_STP')) ##  Not Significant
summary(model2)
anova(model2)
model3 <- lm(Slope ~ ALT,data = Site_data %>% filter(Type=='Slope_STP')) ## Marginally Significant
summary(model3)
anova(model3)
model4 <- lm(Slope ~ MAT,data = Site_data %>% filter(Type=='Slope_STP')) ## Significant
summary(model4)
anova(model4)
model5 <- lm(Slope ~ MAP,data = Site_data %>% filter(Type=='Slope_STP')) ## Marginally Significant
summary(model5)
anova(model5)



model1 <- lm(Slope ~ LON,data = Site_data %>% filter(Type=='Slope_ITV')) ## Not Significant
summary(model1)
anova(model1)
model2 <- lm(Slope ~ LAT,data = Site_data %>% filter(Type=='Slope_ITV')) ## Not Significant
summary(model2)
anova(model2)
model3 <- lm(Slope ~ ALT,data = Site_data %>% filter(Type=='Slope_ITV')) ## Marginally Significant
summary(model3)
anova(model3)
model4 <- lm(Slope ~ MAT,data = Site_data %>% filter(Type=='Slope_ITV')) ## Significant
summary(model4)
anova(model4)
model5 <- lm(Slope ~ MAP,data = Site_data %>% filter(Type=='Slope_ITV')) ## Not Significant
summary(model5)
anova(model5)



model1 <- lm(log(rate) ~ LON,data = env) ## Significant
summary(model1)
anova(model1)
model2 <- lm(log(rate) ~ LAT,data = env) ## Marginally significant
summary(model2)
anova(model2)
model3 <- lm(log(rate) ~ ALT,data = env) ## Significant
summary(model3)
anova(model3)
model4 <- lm(log(rate) ~ MAT,data = env) ## Significant
summary(model4)
anova(model4)
model5 <- lm(log(rate) ~ MAP,data = env) ## Not Significant
summary(model5)
anova(model5)


#### Linear mixed-effects model and General linear model for fig 3 #### 

model1 <- lm(NestedT ~ LON,data = env) ## Significant
summary(model1)
anova(model1)
model2 <- lm(NestedT ~ LAT,data = env) ##  Significant
summary(model2)
anova(model2)
model3 <- lm(NestedT ~ ALT,data = env) ##  Significant
summary(model3)
anova(model3)
model4 <- lm(NestedT ~ MAT,data = env) ## Significant
summary(model4)
anova(model4)
model5 <- lm(NestedT ~ MAP,data = env) ## Not Significant
summary(model5)
anova(model5)


data_PL_env <- env %>% select(Site,LON,LAT,ALT,MAT,MAP) %>% left_join(data_all %>% select(Site,PL,FishZ_PL),multiple = "all" )

model1 <- lmer(FishZ_PL ~ LON+(1|Site),data = data_PL_env) ## Significant
summary(model1)
anova(model1)
model2 <- lmer(FishZ_PL ~ LAT+(1|Site),data = data_PL_env) ##  Marginally Significant
summary(model2)
anova(model2)
model3 <- lmer(FishZ_PL ~ log(ALT)+(1|Site),data = data_PL_env) ## Significant
summary(model3)
anova(model3)
model4 <- lmer(FishZ_PL ~ MAT+(1|Site),data = data_PL_env) ## Significant
summary(model4)
anova(model4)
model5 <- lmer(FishZ_PL ~ MAP+(1|Site),data = data_PL_env) ## Significant
summary(model5)
anova(model5)


#### SEM for Fig 4 a  ####

sem1 <- env %>% dplyr::select(MAT,MAP,Site) %>% left_join(data_all,multiple = "all") %>% mutate(SR = scale(SR),
                                                                                             MAT = scale(MAT),
                                                                                             MAP = scale(MAP),
                                                                                             SR_MAT=scale(SR)*scale(MAT),
                                                                                             SR_MAP=scale(SR)*scale(MAP))

model <- psem(lme(FishZ_Turn~SR+MAT+MAP+SR_MAP+SR_MAT,random= ~ 1|Site,sem1),
              lme(FishZ_ITV_PL~SR+MAT+MAP+SR_MAP+SR_MAT,random= ~ 1|Site,sem1),
              FishZ_Turn%~~%FishZ_ITV_PL)
summary(model)
piecewiseSEM::rsquared(model,method='theoretical')

lme(FishZ_Turn~SR+MAT+MAP+SR_MAP+SR_MAT,random= ~ 1|Site,sem1) %>% car::vif()
lme(FishZ_ITV_PL~SR+MAT+MAP+SR_MAP+SR_MAT,random= ~ 1|Site,sem1) %>% car::vif()



####   Partial regression plot for Fig 4 b  ####
## Species turnover ~ SR:MAT
# yresid <- resid(lmer(FishZ_Turn~SR+MAT+MAP+SR_MAP+(1|Site),sem1))
# xresid <- resid(lmer(SR_MAT~SR+MAP+SR_MAP+MAT+(1|Site),sem1))
resid1 <- partialResid(FishZ_Turn ~ SR_MAT,model)
lm(yresid~xresid,resid1) %>% summary()

p1 <- resid1 %>% ggplot(aes(xresid,yresid))+
  geom_point(size=3,shape=1,color="#F39B7FFF")+
  geom_smooth(method = 'lm',color='#0225F8',linewidth=0.5)+
  mytheme+
  labs(x='',y='Species turnover effect | Others')

## Species turnover ~ SR:MAP
# yresid <- resid(lmer(FishZ_Turn~SR+MAT+MAP+SR_MAT+(1|Site),sem1))
# xresid <- resid(lmer(SR_MAP~SR+MAP+SR_MAT+MAT+(1|Site),sem1))
resid1 <- partialResid(FishZ_Turn ~ SR_MAP,model)
lm(yresid~xresid,resid1) %>% summary()

p2 <- resid1 %>% ggplot(aes(xresid,yresid))+
  geom_point(size=3,shape=1,color="#F39B7FFF")+
  geom_smooth(method = 'lm',color='#0CA133',linewidth=0.5)+
  mytheme+
  labs(x='',y='')

# yresid <- resid(lmer(FishZ_ITV_PL~SR+MAT+MAP+SR_MAP+(1|Site),sem1))
# xresid <- resid(lmer(SR_MAT~SR+MAP+SR_MAP+MAT+(1|Site),sem1))
resid1 <- partialResid(FishZ_ITV_PL ~ SR_MAT,model)
lm(yresid~xresid,resid1) %>% summary()

p3 <- resid1 %>% ggplot(aes(xresid,yresid))+
  geom_point(size=3,shape=1,color="#8491B4FF")+
  geom_smooth(method = 'lm',color='#0225F8',linewidth=0.5)+
  mytheme+
  labs(x='SR:MAT | Others',y='Intraspecific variation | Others')

# yresid <- resid(lmer(FishZ_ITV_PL~SR+MAT+MAP+SR_MAT+(1|Site),sem1))
# xresid <- resid(lmer(SR_MAP~SR+MAP+MAT+SR_MAT+(1|Site),sem1))
resid1 <- partialResid(FishZ_ITV_PL ~ SR_MAP,model)
lm(yresid~xresid,resid1) %>% summary()

p4 <- resid1 %>% ggplot(aes(xresid,yresid))+
  geom_point(size=3,shape=1,color="#8491B4FF")+
  #geom_smooth(method = 'lm',color='#0225F8')+
  mytheme+
  labs(x='SR:MAP | Others',y='')

p0 <- ggarrange(p1,p2,p3,p4,nrow = 2,ncol = 2,align = 'hv',
                label.x = 0.23,label.y = 0.95,font.label = list(size = 12, face = "bold"))
p0
ggsave('SEM_partial2.tiff',height = 6,width = 6)

