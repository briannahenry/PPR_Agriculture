larvbio<-larvbiomass_transhab
larvbio$yr0f<-as.factor(larvbio$yr)
larvbio$trt0f<-as.factor(larvbio$trt)

mlbio<-brm(tot_transhab~trt0f*yr0f+(1|site/trans/hab),data=larvbio,family=Gamma(link="log"),
          prior=c(prior(normal(0,4),class="b"),
                  prior(normal(1,6),class="Intercept"),
                  prior(cauchy(0,1),class="sd")),
          chains=4,iter=2000)

print(mlbio)
pp_check(mlbio, type="boxplot")
marginal_effects(mlbio)

tdlbplot<-marginal_effects(mlbio, robust=FALSE)
tdlbplot2<-as.data.frame(tdlbplot$`trt0f:yr0f`)
tdlbplot2$est<-(tdlbplot2$estimate__)/1.5
tdlbplot2$lower<-(tdlbplot2$lower__)/1.5
tdlbplot2$upper<-(tdlbplot2$upper__)/1.5

larvbioplot<-ggplot(tdlbplot2, aes(x=trt0f, y=est, ymin=lower, ymax=upper, fill=yr0f))+
  geom_point(size=4,position=position_dodge(width=0.4), shape=21)+
  geom_errorbar(width=0.1, position=position_dodge(width=0.4))+
  xlab("Treatment")+
  ylab(expression(paste("Larval Biomass (mgDM/m"^2,")")))+
  scale_fill_manual(labels = c("2015", "2016"),  values = c("grey0","grey60"), name="Sample Year")+
  theme_classic()+
  theme(axis.text.x = element_text(color="black", size="10"))+
  theme(axis.text.y = element_text(color="black", size="10"))+
  theme(legend.title = element_text(face="bold", size="10"))

larvbioplot

ggsave("LarvalBiomass.tiff", larvbioplot, dpi=400, width=5, height=3, units="in")

library(brms)

#Differences among treatments
#2016 
mlbpost<-posterior_samples(mlbio)
str(mlbpost)
bc15<-mlbpost$b_Intercept
bt15<-mlbpost$b_Intercept+mlbpost$b_trt0fTile
bs15<-mlbpost$b_Intercept+mlbpost$b_trt0fSurface
  
bc16<-mlbpost$b_Intercept+mlbpost$b_yr0f2016
bt16<-mlbpost$b_Intercept+mlbpost$b_trt0fTile+mlbpost$`b_trt0fTile:yr0f2016`
bs16<-mlbpost$b_Intercept+mlbpost$b_trt0fSurface+mlbpost$`b_trt0fSurface:yr0f2016`

#2015
#surface and control
bsc15<-(exp(bs15)/1.5)-(exp(bc15)/1.5)
mean(bsc15)
quantile(bsc15, probs=c(0.025,0.975))
sum(bsc15>0)/4000 #87%

#surface and tile
bst15<-(exp(bs15)/1.5)-(exp(bt15)/1.5)
mean(bst15)
quantile(bst15, probs=c(0.025,0.975))
sum(bst15>0)/4000 #94%

#control and tile
bct15<-(exp(bc15)/1.5)-(exp(bt15)/1.5)
mean(bct15)
quantile(bct15, probs=c(0.025,0.975))
sum(bct15>0)/4000 #69%

(exp(mean(bs15))/1.5)#120
(exp(mean(bc15))/1.5)#81
(exp(mean(bt15))/1.5)#67

#2016
#surface and control
bsc16<-(exp(bs16)/1.5)-(exp(bc16)/1.5)
mean(bsc16)
quantile(bsc16, probs=c(0.025,0.975))
sum(bsc16>0)/4000 #99.8%

#surface and tile
bst16<-(exp(bs16)/1.5)-(exp(bt16)/1.5)
mean(bst16)
quantile(bst16, probs=c(0.025,0.975))
sum(bst16>0)/4000 #56%

#control and tile
bct16<-(exp(bt16)/1.5)-(exp(bc16)/1.5)
mean(bct16)
quantile(bct16, probs=c(0.025,0.975))
sum(bct16>0)/4000 #99.6%

(exp(mean(bs16))/1.5)#148
(exp(mean(bc16))/1.5)#35
(exp(mean(bt16))/1.5)#139
