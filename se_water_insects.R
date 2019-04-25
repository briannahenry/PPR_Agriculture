library(brms)

library(ggplot2)



#####THIS IS OUR FINAL MODEL#######

sellmmay<-brm(se_ugg~log(watermay)+(1|samp), data=se_bugs_water, family = Gamma(link="log"),
              
              prior=c(prior(normal(0,2),class="b"),
                      
                      prior(normal(0,3),class="Intercept")),
              
              chains=4, iter=4000, control = list(adapt_delta=0.99))



pp_check(sellmmay, type="boxplot") #posterior predictive distribution. New data that are simulated from the posterior looks like original data, which is good

print(sellmmay) #rhats are good. All below 1.1

postsellmm<-posterior_samples(sellmmay) #extract posteriors

sellmm_sl<-exp(postsellmm$b_logwatermay) #slope of the model

sum(sellmm_sl>0)/8000 #Probability that slope is greater than zero. That probability is >0.99.

quantile(sellmm_sl,probs=c(0.025,0.5,0.975)) #median and credible intervals of the slope.



#Plot it.

sellmmayplot<-marginal_effects(sellmmay, robust=FALSE)



sellmplot<-marginal_effects(sellmmay, robust=FALSE)

sellmplot<-as.data.frame(sellmplot$watermay)



se_bugs_water$stage2<-ifelse(se_bugs_water$stage=="larvae", "Larva", "Adult")



sebugwaterp<-ggplot()+
  
  geom_line(data=sellmplot, aes(x=watermay, y=estimate__))+
  
  geom_ribbon(data=sellmplot, aes(x=watermay, ymin=lower__, ymax=upper__), alpha=0.2)+
  
  geom_point(data=se_bugs_water,aes(x=waterjune,y=se_ugg,fill=stage2), shape=21, alpha=0.6)+
  
  scale_fill_manual(name="Stage", labels=c("Larva", "Adult"), values=c("white", "grey0"))+
  
  theme_classic()+
  
  scale_x_continuous(limits=c(0,4))+
  
  #scale_x_log10()+
  
  #scale_y_log10()+
  
  xlab(expression(paste("May Aqueous Selenium Concentrations (", mu, "g/l)")))+
  
  ylab(expression(paste("Insect Tissue Selenium Concentrations (",mu,"g/g)")))



sebugwaterp



ggsave("SeBugsWaterRegress_scaleLimit.tiff", sebugwaterp, dpi=400, width=5, height=4, units="in")

#Model comparison with water in June#
sellm<-brm(se_ugg~log(waterjune)+(1|samp), data=se_bugs_water, family = Gamma(link="log"),
           prior=c(prior(normal(0,2),class="b"),
                   prior(normal(0,4),class="Intercept")),
           chains=4, iter=4000,  control = list(adapt_delta=0.99))

print(sellm)

WAIC(sellm, sellmmay)
