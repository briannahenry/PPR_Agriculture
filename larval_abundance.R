library(brms)
library(ggplot2)

inv<-transhab_inverts
inv$trt0f<-as.factor(inv$trt)
inv$yr0f<-as.factor(inv$year)
######Model/Plot######
tdlam<-brm(totno~trt0f*yr0f+(1|site/trans/hab), data=inv, family=poisson(link="log"),
      prior=c(prior(normal(0,8),class="b"),
          prior(normal(4,6),class="Intercept"),
          prior(cauchy(0,1),class="sd")),
      chains=4,iter=2000)


pp_check(tdlam,type="boxplot") 

print(tdlam)

marginal_effects(tdlam, robust=FALSE)

#units: each is # bugs/ 10 sweeps and one sweep is 0.15m2 so each is # bugs/1.5 m2
##so divide by 1.5

tdlaplot<-marginal_effects(tdlam, robust=FALSE)
tdlaplot2<-as.data.frame(tdlaplot$`trt0f:yr0f`)
tdlaplot2$est<-(tdlaplot2$estimate__)/1.5
tdlaplot2$lower<-(tdlaplot2$lower__)/1.5
tdlaplot2$upper<-(tdlaplot2$upper__)/1.5


junelarvabund2<-ggplot(tdlaplot2, aes(x=trt0f, y=est, ymin=lower, ymax=upper, fill=yr0f))+
  geom_point(size=4,position=position_dodge(width=0.4), shape=21)+
  geom_errorbar(width=0.1, position=position_dodge(width=0.4))+
  xlab("Treatment")+
  ylab(expression(paste("Larval Abundance (#/m"^2,")")))+
  scale_fill_manual(labels = c("2015", "2016"),  values = c("grey0","grey60"), name="Sample Year")+
  theme_classic()+
  theme(axis.text.x = element_text(color="black", size="10"))+
  theme(axis.text.y = element_text(color="black", size="10"))+
  theme(legend.title = element_text(face="bold", size="10"))

junelarvabund2


ggsave("LarvalAbundance3.tiff", junelarvabund2, dpi=400, width=5, height=3, units="in")


#####Differences in treatments#######
larvpost<-posterior_samples(tdlam)
str(larvpost)

#2015
c_15<-larvpost$b_Intercept
s_15<-larvpost$b_Intercept+larvpost$b_trt0fSurface
t_15<-larvpost$b_Intercept+larvpost$b_trt0fTile

#Control and surface
lcs_15<-exp(s_15)-exp(c_15)
mean(lcs_15)
quantile(lcs_15, probs=c(0.025,0.975))
sum(lcs_15>0)/4000 #75.3%

#control and tile
lct_15<-exp(t_15)-exp(c_15)
mean(lct_15)
quantile(lct_15, probs=c(0.025,0.975))
sum(lct_15>0)/4000 #65%

#surface and tile
lst_15<-exp(s_15)-exp(t_15)
mean(lst_15)
quantile(lst_15, probs=c(0.025,0.975))
sum(lst_15>0)/4000 #61.6

(exp(mean(c_15)))/1.5#48.5
(exp(mean(s_15)))/1.5#61.94
(exp(mean(t_15)))/1.5#55.8

c_16<-larvpost$b_Intercept+larvpost$b_yr0f2016
s_16<-larvpost$b_Intercept+larvpost$b_yr0f2016+larvpost$b_trt0fSurface+larvpost$`b_trt0fSurface:yr0f2016`
t_16<-larvpost$b_Intercept+larvpost$b_yr0f2016+larvpost$b_trt0fTile+larvpost$`b_trt0fTile:yr0f2016`

#Control and surface
lcs_16<-exp(s_16)-exp(c_16)
mean(lcs_16)
quantile(lcs_16, probs=c(0.025,0.975))
sum(lcs_16>0)/4000 #98.3

#control and tile
lct_16<-exp(t_16)-exp(c_16)
mean(lct_16)
quantile(lct_16, probs=c(0.025,0.975))
sum(lct_16>0)/4000 #99.3%

#surface and tile
lst_16<-exp(t_16)-exp(s_16)
mean(lst_16)
quantile(lst_16, probs=c(0.025,0.975))
sum(lst_16>0)/4000 #63.4

(exp(mean(c_16)))/1.5#15.9
(exp(mean(s_16)))/1.5#37.04
(exp(mean(t_16)))/1.5#41.9
#surface and tile treatments are both more than twice the control trt
