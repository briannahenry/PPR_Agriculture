library(brms)
library(ggplot2)



se_adult$trt0f<-as.factor(se_adult$trt)
se_adult$mon_yr2<-as.factor(se_adult$mon_yr)
se_adult$yr<-as.factor(se_adult$year)

#selenium adults by treatment
get_prior(se_ugg~trt*yr+(1|site/rep), data=se_adult, family=Gamma(link="log"))

seadu_trt<-brm(se_ugg~trt0f*yr+(1|site/rep), data=se_adult, family=Gamma(link="log"),
               prior=c(prior(normal(0,3),class="b"),
                       prior(normal(0,2),class="Intercept"),
                       prior(cauchy(0,1), class="sd")),
               chains=4,iter=2000)
pp_check(seadu_trt, type="boxplot")
print(seadu_trt)
marginal_effects(seadu_trt)

#plot
seatrtplot<-marginal_effects(seadu_trt, robust=FALSE)
seatrtplot<-as.data.frame(seatrtplot$`trt0f:yr`)
seatrtplot$trt2<-ifelse(seatrtplot$trt0f=="control", "Control",
                       ifelse(seatrtplot$trt0f=="surface", "Surface", "Tile"))
#Add lines for risk and tox
se_trtp<-ggplot(seatrtplot, aes(x=trt2, y=estimate__, ymin=lower__, ymax=upper__, fill=yr))+
  geom_point(size=4, position=position_dodge(width=0.4), shape=21)+
  geom_errorbar(width=0.1, position=position_dodge(width=0.4))+
  xlab("Treatment")+
  ylab(expression(paste("Selenium Body Burden (",mu,"g/g)")))+
  theme_classic()+
  scale_fill_manual(name="Sample Year", labels = c("2015","2016"),  values = c("grey0","grey60"))+
  theme(axis.text.x = element_text(color="black", size="10"))+
  theme(axis.text.y = element_text(color="black", size="10"))+
  theme(legend.title = element_text(face="bold", size="10"))

se_trtp

ggsave("Se_adultbytrt.tiff", se_trtp, dpi=400, width=5, height=3, units="in")


#Differences by treatment and year
seapost<-posterior_samples(seadu_trt)
str(seapost)

seac15<-seapost$b_Intercept
seas15<-seapost$b_Intercept+seapost$b_trt0fsurface
seat15<-seapost$b_Intercept+seapost$b_trt0ftile

#Control and surface
seacs_15<-exp(seac15)-exp(seas15)
mean(seacs_15)
quantile(seacs_15, probs=c(0.025,0.975))
sum(seacs_15>0)/4000 #78.4%

#surface and tile
seast_15<-exp(seat15)-exp(seas15)
mean(seast_15)
quantile(seast_15, probs=c(0.025,0.975))
sum(seast_15>0)/4000 #99.4%

#control and tile
seact_15<-exp(seat15)-exp(seac15)
mean(seact_15)
quantile(seact_15, probs=c(0.025,0.975))
sum(seact_15>0)/4000 #97.4%

exp(mean(seac15))#1.5
exp(mean(seas15))#1.19
exp(mean(seat15))#2.9
#Tile is 2 times higher than control
#Almost 2.5 times higher than surface

seac16<-seapost$b_Intercept+seapost$b_yr2016
seas16<-seapost$b_Intercept+seapost$b_trt0fsurface+seapost$b_yr2016+seapost$`b_trt0fsurface:yr2016`
seat16<-seapost$b_Intercept+seapost$b_trt0ftile+seapost$b_yr2016+seapost$`b_trt0ftile:yr2016`

#Control and surface
seacs_16<-exp(seas16)-exp(seac16)
mean(seacs_16)
quantile(seacs_16, probs=c(0.025,0.975))
sum(seacs_16>0)/4000 #88.8%

#surface and tile
seast_16<-exp(seat16)-exp(seas16)
mean(seast_16)
quantile(seast_16, probs=c(0.025,0.975))
sum(seast_16>0)/4000 #98.7%

#control and tile
seact_16<-exp(seat16)-exp(seac16)
mean(seact_16)
quantile(seact_16, probs=c(0.025,0.975))
sum(seact_16>0)/4000 #99.9%

exp(mean(seac16))#1.55
exp(mean(seas16))#2.34
exp(mean(seat16))#5.04

##########Adults vs water########

seadu_water1<-brm(se_ugg~watermay*yr+(1|site/rep), data=se_adult, family=Gamma(link="log"),
               prior=c(prior(normal(0,3),class="b"),
                       prior(normal(0,2),class="Intercept"),
                       prior(cauchy(0,1), class="sd")),
               chains=1,iter=2000)

print(seadu_water1)
pp_check(seadu_water1, type="boxplot")

seadu_water2<-brm(se_ugg~waterjune*yr+(1|site/rep), data=se_adult, family=Gamma(link="log"),
                  prior=c(prior(normal(0,3),class="b"),
                          prior(normal(0,2),class="Intercept"),
                          prior(cauchy(0,1), class="sd")),
                  chains=4,iter=2000)

print(seadu_water2)
pp_check(seadu_water2, type="boxplot")

WAIC(seadu_water1,seadu_water2) #no difference

seawatplot<-marginal_effects(seadu_water2, robust=FALSE)
seawatplot<-as.data.frame(seawatplot$`waterjune:yr`)

ggplot(seawatplot, aes(x=waterjune, y=estimate__, ymin=lower__, ymax=upper__))+
  geom_ribbon(alpha=0.3)+
  geom_line(alpha=0.3)+
  scale_y_log10()+
  theme_classic()+
  facet_wrap(~yr)
