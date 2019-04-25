library(brms)
library(ggplot2)

eee2$trt0f<-as.factor(eee2$trt)

eee2$mon_year2<-paste(eee2$month,"",eee2$year)


####Adult Abundance####
tdam<-brm(abund~trt0f*mon_year+(1|site/trapid),data=eee2,family=poisson(link="log"),
          prior=c(prior(normal(0,2),class="b"),
                  prior(normal(10,3),class="Intercept"),
                  prior(cauchy(0,1),class="sd")))

pp_check(tdam,type="boxplot")

print(tdam)

marginal_effects(tdam, robust=FALSE)

tdaplot<-marginal_effects(tdam, robust=FALSE)
tdaplot<-as.data.frame(tdaplot$`trt0f:mon_year`)
tdaplot$date<-factor(tdaplot$mon_year, levels = c("2015 May","2015 June", "2015 July", "2016 May", "2016 June", "2016 July"))
tdaplot$est2<-(tdaplot$estimate__)*2.78/4
tdaplot$upper2<-(tdaplot$upper__)*2.78/4
tdaplot$lower2<-(tdaplot$lower__)*2.78/4
tdaplot$year<-ifelse(tdaplot$mon_year=="2015 May", "2015",
                     ifelse(tdaplot$mon_year=="2015 June", "2015",
                            ifelse(tdaplot$mon_year=="2015 July", "2015", "2016")))
tdaplot$mth<-ifelse(tdaplot$mon_year=="2015 May", "May",
                    ifelse(tdaplot$mon_year=="2016 May", "May",
                           ifelse(tdaplot$mon_year=="2015 June", "June",
                                  ifelse(tdaplot$mon_year=="2016 June", "June", "July"))))
tdaplot$mth2<-factor(tdaplot$mth, levels=c("May", "June", "July"))

tdadultabunplot<-ggplot(tdaplot, aes(x=mth2, y=est2, ymin=lower2, ymax=upper2, color=trt0f, group=trt0f))+
  geom_point(position=pd, size=5)+
  geom_errorbar(width=0.1, position=pd)+
  geom_line(position=pd, size=5, alpha=0.2)+
  xlab("Date Sampled")+
  ylab(expression(paste("Adult Abundance (#/m"^2,"/day)")))+
  scale_color_manual(name="Treatment", labels = c("Control", "Surface", "Tile"),  values = c("grey0","dodgerblue4","tomato4"))+
  theme_classic()+
  facet_wrap(~year)+
  theme(axis.text.x = element_text(color="black"))+
  theme(axis.text.y = element_text(color="black"))+
  theme(legend.title = element_text(face="bold", size="11"))

tdadultabunplot      

ggsave("TileDrain_AdultAbund.tiff", tdadultabunplot, dpi=400, width=7, height=4.5, units="in")

##########Trt_Date Comparisions#######
aduabundpost<-posterior_samples(tdam)
str(aduabundpost)

#may 2015
con_515<-aduabundpost$b_Intercept+aduabundpost$b_mon_year2015May
sur_515<-aduabundpost$b_Intercept+aduabundpost$b_mon_year2015May+aduabundpost$b_trt0fSurface+aduabundpost$`b_trt0fSurface:mon_year2015May`
tile_515<-aduabundpost$b_Intercept+aduabundpost$b_mon_year2015May+aduabundpost$b_trt0fTile+aduabundpost$`b_trt0fTile:mon_year2015May`
  
#control and surface
con_sur_515<-((exp(con_515)*2.78)/4)-((exp(sur_515)*2.78)/4)
mean(con_sur_515)
quantile(con_sur_515, probs=c(0.025,0.975))
sum(con_sur_515>0)/4000 #77.3% prob

#control and tile
con_tile_515<-((exp(tile_515)*2.78)/4)-((exp(con_515)*2.78)/4)
mean(con_tile_515)
quantile(con_tile_515, probs=c(0.025,0.975))
sum(con_tile_515>0)/4000 #49.6% prob

#Surface and Tile
sur_tile_515<-((exp(tile_515)*2.78)/4)-((exp(sur_515)*2.78)/4)
mean(sur_tile_515)
quantile(sur_tile_515, probs=c(0.025,0.975))
sum(sur_tile_515>0)/4000 #23.4% prob

#June 2015
con_615<-aduabundpost$b_Intercept+aduabundpost$b_mon_year2015June
sur_615<-aduabundpost$b_Intercept+aduabundpost$b_trt0fSurface+aduabundpost$b_mon_year2015June+aduabundpost$`b_trt0fSurface:mon_year2015June`
tile_615<-aduabundpost$b_Intercept+aduabundpost$b_trt0fTile+aduabundpost$b_mon_year2015June+aduabundpost$`b_trt0fTile:mon_year2015June`

#control and surface
con_sur_615<-((exp(sur_615)*2.78)/4)-((exp(con_615)*2.78)/4)
mean(con_sur_615)
quantile(con_sur_615, probs=c(0.025,0.975))
sum(con_sur_615>0)/4000 #81.9% prob

#control and tile
con_tile_615<-((exp(con_615)*2.78)/4)-((exp(tile_615)*2.78)/4)
mean(con_tile_615)
quantile(con_tile_615, probs=c(0.025,0.975))
sum(con_tile_615>0)/4000 #97.55% prob

#Surface and Tile
sur_tile_615<-((exp(sur_615)*2.78)/4)-((exp(tile_615)*2.78)/4)
mean(sur_tile_615)
quantile(sur_tile_615, probs=c(0.025,0.975))
sum(sur_tile_615>0)/4000 #99.7% prob

  #Percent differences
  mean(exp(con_615))*2.78/4#17.5
  mean(exp(sur_615))*2.78/4#36.7
  mean(exp(tile_615))*2.78/4#7.4
  #Con is 2.3 times higher than tile
  #Sur is almost 5 times higher than tile

#July 2015
con_715<-aduabundpost$b_Intercept
sur_715<-aduabundpost$b_Intercept+aduabundpost$b_trt0fSurface
tile_715<-aduabundpost$b_Intercept+aduabundpost$b_trt0fTile

#control and surface
con_sur_715<-((exp(sur_715)*2.78)/4)-((exp(con_715)*2.78)/4)
mean(con_sur_715)
quantile(con_sur_715, probs=c(0.025,0.975))
sum(con_sur_715>0)/4000 #92% prob

#control and tile
con_tile_715<-((exp(con_715)*2.78)/4)-((exp(tile_715)*2.78)/4)
mean(con_tile_715)
quantile(con_tile_715, probs=c(0.025,0.975))
sum(con_tile_715>0)/4000 #96% prob

#Surface and Tile
sur_tile_715<-((exp(sur_715)*2.78)/4)-((exp(tile_715)*2.78)/4)
mean(sur_tile_715)
quantile(sur_tile_715, probs=c(0.025,0.975))
sum(sur_tile_715>0)/4000 #99% prob

  #Percent differences
  mean(exp(con_715))*2.78/4#6.5
  mean(exp(sur_715))*2.78/4#12.9
  mean(exp(tile_715))*2.78/4#2.6
  #Sur is 2 times higher than con
  #Con is 3 times higher than tile
  #Sur is 6 times higher than tile

#May 2016
con_516<-aduabundpost$b_Intercept+aduabundpost$b_mon_year2016May
sur_516<-aduabundpost$b_Intercept+aduabundpost$b_trt0fSurface+aduabundpost$b_mon_year2016May+aduabundpost$`b_trt0fSurface:mon_year2016May`
tile_516<-aduabundpost$b_Intercept+aduabundpost$b_trt0fTile+aduabundpost$b_mon_year2016May+aduabundpost$`b_trt0fTile:mon_year2016May`

#control and surface
con_sur_516<-((exp(sur_516)*2.78)/4)-((exp(con_516)*2.78)/4)
mean(con_sur_516)
quantile(con_sur_516, probs=c(0.025,0.975))
sum(con_sur_516>0)/4000 #57% prob

#control and tile
con_tile_516<-((exp(tile_516)*2.78)/4)-((exp(con_516)*2.78)/4)
mean(con_tile_516)
quantile(con_tile_516, probs=c(0.025,0.975))
sum(con_tile_516>0)/4000 #93% prob

#Surface and Tile
sur_tile_516<-((exp(tile_516)*2.78)/4)-((exp(sur_516)*2.78)/4)
mean(sur_tile_516)
quantile(sur_tile_516, probs=c(0.025,0.975))
sum(sur_tile_516>0)/4000 #91% prob

  #Percent differences
  mean(exp(con_516))*2.78/4#13.6
  mean(exp(sur_516))*2.78/4#14.9
  mean(exp(tile_516))*2.78/4#28.9
  #Tile abundance is twice that of the surface or control treatments

#June 2016
con_616<-aduabundpost$b_Intercept+aduabundpost$b_mon_year2016June
sur_616<-aduabundpost$b_Intercept+aduabundpost$b_mon_year2016June+aduabundpost$b_trt0fSurface+aduabundpost$`b_trt0fSurface:mon_year2016June`
tile_616<-aduabundpost$b_Intercept+aduabundpost$b_mon_year2016June+aduabundpost$b_trt0fTile+aduabundpost$`b_trt0fTile:mon_year2016June`

#control and surface
con_sur_616<-((exp(sur_616)*2.78)/4)-((exp(con_616)*2.78)/4)
mean(con_sur_616)
quantile(con_sur_616, probs=c(0.025,0.975))
sum(con_sur_616>0)/4000 #68.9 prob

#control and tile
con_tile_616<-((exp(tile_616)*2.78)/4)-((exp(con_616)*2.78)/4)
mean(con_tile_616)
quantile(con_tile_616, probs=c(0.025,0.975))
sum(con_tile_616>0)/4000 #53.8% prob

#Surface and Tile
sur_tile_616<-((exp(sur_616)*2.78)/4)-((exp(tile_616)*2.78)/4)
mean(sur_tile_616)
quantile(sur_tile_616, probs=c(0.025,0.975))
sum(sur_tile_616>0)/4000 #36.4% prob

mean(exp(con_616))*2.78/4#6.03
mean(exp(sur_616))*2.78/4#7.5
mean(exp(tile_616))*2.78/4#6.4

#July 2016
con_716<-aduabundpost$b_Intercept+aduabundpost$b_mon_year2016July
sur_716<-aduabundpost$b_Intercept+aduabundpost$b_mon_year2016July+aduabundpost$b_trt0fSurface+aduabundpost$`b_trt0fSurface:mon_year2016July`
tile_716<-aduabundpost$b_Intercept+aduabundpost$b_mon_year2016July+aduabundpost$b_trt0fTile+aduabundpost$`b_trt0fTile:mon_year2016July`

#control and surface
con_sur_716<-((exp(sur_716)*2.78)/4)-((exp(con_716)*2.78)/4)
mean(con_sur_716)
quantile(con_sur_716, probs=c(0.025,0.975))
sum(con_sur_716>0)/4000 #52.3 prob

#control and tile
con_tile_716<-((exp(con_716)*2.78)/4)-((exp(tile_716)*2.78)/4)
mean(con_tile_716)
quantile(con_tile_716, probs=c(0.025,0.975))
sum(con_tile_716>0)/4000 #38.9% prob

#Surface and Tile
sur_tile_716<-((exp(sur_716)*2.78)/4)-((exp(tile_716)*2.78)/4)
mean(sur_tile_716)
quantile(sur_tile_716, probs=c(0.025,0.975))
sum(sur_tile_716>0)/4000 #37.5% prob

mean(exp(con_716))*2.78/4#10.7
mean(exp(sur_716))*2.78/4#11.2
mean(exp(tile_716))*2.78/4#9.35

#####Adult Biomass####
tdbm<-brm(mg01~trt0f*mon_year+(1|site/trapid),data=eee2,family=Gamma(link="log"),
          prior=c(prior(normal(0,4),class="b"),
                  prior(normal(1,6),class="Intercept"),
                  prior(cauchy(0,1),class="sd")),
          chains=4,iter=2000)

pp_check(tdbm,type="boxplot") #Potentially need to restrict priors more

print(tdbm)

marginal_effects(tdbm, robust=FALSE)

tdbplot<-marginal_effects(tdbm, robust=FALSE)
tdbplot<-as.data.frame(tdbplot$`trt0f:mon_year`)
tdbplot$date<-factor(tdbplot$mon_year, levels = c("2015 May","2015 June", "2015 July", "2016 May", "2016 June", "2016 July"))
tdbplot$est2<-(tdbplot$estimate__)*2.78/4
tdbplot$upper2<-(tdbplot$upper__)*2.78/4
tdbplot$lower2<-(tdbplot$lower__)*2.78/4
tdbplot$year<-ifelse(tdbplot$mon_year=="2015 May", "2015",
                     ifelse(tdbplot$mon_year=="2015 June", "2015",
                            ifelse(tdbplot$mon_year=="2015 July", "2015", "2016")))
tdbplot$mth<-ifelse(tdbplot$mon_year=="2015 May", "May",
                    ifelse(tdbplot$mon_year=="2016 May", "May",
                           ifelse(tdbplot$mon_year=="2015 June", "June",
                                  ifelse(tdbplot$mon_year=="2016 June", "June", "July"))))
tdbplot$mth2<-factor(tdbplot$mth, levels=c("May", "June", "July"))

tdadultbioplot<-ggplot(tdbplot, aes(x=mth2, y=est2, ymin=lower2, ymax=upper2, color=trt0f, group=trt0f))+
  geom_point(position=pd, size=5)+
  geom_errorbar(width=0.1, position=pd)+
  geom_line(position=pd, size=5, alpha=0.2)+
  xlab("Date Sampled")+
  ylab(expression(paste("Adult Biomass (mgDM/m"^2,"/day)")))+
  scale_color_manual(name="Treatment", labels = c("Control", "Surface", "Tile"),  values = c("grey0","dodgerblue4","tomato4"))+
  theme_classic()+
  facet_wrap(~year)+
  theme(axis.text.x = element_text(color="black"))+
  theme(axis.text.y = element_text(color="black"))+
  theme(legend.title = element_text(face="bold", size="11"))

tdadultbioplot

ggsave("TilDrain_AdultBiomass.tiff", tdadultbioplot, dpi=400, width=7, height=4.5, units="in")

######Treatment and date comparisons#######
tdbpost<-posterior_samples(tdbm)
str(tdbpost)

#May 2015
conb_515<-tdbpost$b_Intercept+tdbpost$b_mon_year2015May
surb_515<-tdbpost$b_Intercept+tdbpost$b_mon_year2015May+tdbpost$`b_trt0fSurface:mon_year2015May`+tdbpost$b_trt0fSurface
tileb_515<-tdbpost$b_Intercept+tdbpost$b_mon_year2015May+tdbpost$`b_trt0fTile:mon_year2015May`+tdbpost$b_trt0fTile

#control and surface
bcon_sur_515<-((exp(surb_515)*2.78)/4)-((exp(conb_515)*2.78)/4)
mean(bcon_sur_515)
quantile(bcon_sur_515, probs=c(0.025,0.975))
sum(bcon_sur_515>0)/4000 #54.7% prob

#control and tile
bcon_tile_515<-((exp(tileb_515)*2.78)/4)-((exp(conb_515)*2.78)/4)
mean(bcon_tile_515)
quantile(bcon_tile_515, probs=c(0.025,0.975))
sum(bcon_tile_515>0)/4000 #78.9% prob

#Surface and Tile
bsur_tile_515<-((exp(tileb_515)*2.78)/4)-((exp(surb_515)*2.78)/4)
mean(bsur_tile_515)
quantile(bsur_tile_515, probs=c(0.025,0.975))
sum(bsur_tile_515>0)/4000 #75% prob

#ALL THE SAME

#June 2015
conb_615<-tdbpost$b_Intercept+tdbpost$b_mon_year2015June
surb_615<-tdbpost$b_Intercept+tdbpost$b_mon_year2015June+tdbpost$`b_trt0fSurface:mon_year2015June`+tdbpost$b_trt0fSurface
tileb_615<-tdbpost$b_Intercept+tdbpost$b_mon_year2015June+tdbpost$`b_trt0fTile:mon_year2015June`+tdbpost$b_trt0fTile

#control and surface
bcon_sur_615<-((exp(surb_615)*2.78)/4)-((exp(conb_615)*2.78)/4)
mean(bcon_sur_615)
quantile(bcon_sur_615, probs=c(0.025,0.975))
sum(bcon_sur_615>0)/4000 #75% prob

#control and tile
bcon_tile_615<-((exp(conb_615)*2.78)/4)-((exp(tileb_615)*2.78)/4)
mean(bcon_tile_615)
quantile(bcon_tile_615, probs=c(0.025,0.975))
sum(bcon_tile_615>0)/4000 #85.9% prob

#Surface and Tile
bsur_tile_615<-((exp(surb_615)*2.78)/4)-((exp(tileb_615)*2.78)/4)
mean(bsur_tile_615)
quantile(bsur_tile_615, probs=c(0.025,0.975))
sum(bsur_tile_615>0)/4000 #95.6% prob

exp(mean(conb_615))*2.78/4 #32
exp(mean(surb_615))*2.78/4 #57.9
exp(mean(tileb_615))*2.78/4 #12.5
#Tile is 78% lower than surface

#July 2015
conb_715<-tdbpost$b_Intercept
surb_715<-tdbpost$b_Intercept+tdbpost$b_trt0fSurface
tileb_715<-tdbpost$b_Intercept+tdbpost$b_trt0fTile

#control and surface
bcon_sur_715<-((exp(conb_715)*2.78)/4)-((exp(surb_715)*2.78)/4)
mean(bcon_sur_715)
quantile(bcon_sur_715, probs=c(0.025,0.975))
sum(bcon_sur_715>0)/4000 #78.3% prob

#control and tile
bcon_tile_715<-((exp(conb_715)*2.78)/4)-((exp(tileb_715)*2.78)/4)
mean(bcon_tile_715)
quantile(bcon_tile_715, probs=c(0.025,0.975))
sum(bcon_tile_715>0)/4000 #66.2% prob

#Surface and Tile
bsur_tile_715<-((exp(tileb_715)*2.78)/4)-((exp(surb_715)*2.78)/4)
mean(bsur_tile_715)
quantile(bsur_tile_715, probs=c(0.025,0.975))
sum(bsur_tile_715>0)/4000 #63.3% prob

#May 2016
conb_516<-tdbpost$b_Intercept+tdbpost$b_mon_year2016May
surb_516<-tdbpost$b_Intercept+tdbpost$b_mon_year2016May+tdbpost$b_trt0fSurface+tdbpost$`b_trt0fSurface:mon_year2016May`
tileb_516<-tdbpost$b_Intercept+tdbpost$b_mon_year2016May+tdbpost$b_trt0fTile+tdbpost$`b_trt0fTile:mon_year2016May`

#control and surface
bcon_sur_516<-((exp(conb_516)*2.78)/4)-((exp(surb_516)*2.78)/4)
mean(bcon_sur_516)
quantile(bcon_sur_516, probs=c(0.025,0.975))
sum(bcon_sur_516>0)/4000 #53.9% prob

#control and tile
bcon_tile_516<-((exp(conb_516)*2.78)/4)-((exp(tileb_516)*2.78)/4)
mean(bcon_tile_516)
quantile(bcon_tile_516, probs=c(0.025,0.975))
sum(bcon_tile_516>0)/4000 #86% prob

#Surface and Tile
bsur_tile_516<-((exp(surb_516)*2.78)/4)-((exp(tileb_516)*2.78)/4)
mean(bsur_tile_516)
quantile(bsur_tile_516, probs=c(0.025,0.975))
sum(bsur_tile_516>0)/4000 #82.9% prob

exp(mean(conb_516))*2.78/4 #37.2
exp(mean(surb_516))*2.78/4 #33.2
exp(mean(tileb_516))*2.78/4 #15.1

#ALL THE SAME

#June 2016
conb_616<-tdbpost$b_Intercept+tdbpost$b_mon_year2016June
surb_616<-tdbpost$b_Intercept+tdbpost$b_mon_year2016June+tdbpost$b_trt0fSurface+tdbpost$`b_trt0fSurface:mon_year2016June`
tileb_616<-tdbpost$b_Intercept+tdbpost$b_mon_year2016June+tdbpost$b_trt0fTile+tdbpost$`b_trt0fTile:mon_year2016June`

#control and surface
bcon_sur_616<-((exp(conb_616)*2.78)/4)-((exp(surb_616)*2.78)/4)
mean(bcon_sur_616)
quantile(bcon_sur_616, probs=c(0.025,0.975))
sum(bcon_sur_616>0)/4000 #86.4% prob

#control and tile
bcon_tile_616<-((exp(tileb_616)*2.78)/4)-((exp(conb_616)*2.78)/4)
mean(bcon_tile_616)
quantile(bcon_tile_616, probs=c(0.025,0.975))
sum(bcon_tile_616>0)/4000 #62.5% prob

#Surface and Tile
bsur_tile_616<-((exp(tileb_616)*2.78)/4)-((exp(surb_616)*2.78)/4)
mean(bsur_tile_616)
quantile(bsur_tile_616, probs=c(0.025,0.975))
sum(bsur_tile_616>0)/4000 #92% prob

exp(mean(conb_616))*2.78/4 #27.2
exp(mean(surb_616))*2.78/4 #11.1
exp(mean(tileb_616))*2.78/4 #36.8
#Surface is 70% lower than tile 

#July 2016
conb_716<-tdbpost$b_Intercept+tdbpost$b_mon_year2016July
surb_716<-tdbpost$b_Intercept+tdbpost$b_mon_year2016July+tdbpost$b_trt0fSurface+tdbpost$`b_trt0fSurface:mon_year2016July`
tileb_716<-tdbpost$b_Intercept+tdbpost$b_mon_year2016July+tdbpost$b_trt0fTile+tdbpost$`b_trt0fTile:mon_year2016July`

#control and surface
bcon_sur_716<-((exp(surb_716)*2.78)/4)-((exp(conb_716)*2.78)/4)
mean(bcon_sur_716)
quantile(bcon_sur_716, probs=c(0.025,0.975))
sum(bcon_sur_716>0)/4000 #54.6%

#control and tile
bcon_tile_716<-((exp(tileb_716)*2.78)/4)-((exp(conb_716)*2.78)/4)
mean(bcon_tile_716)
quantile(bcon_tile_716, probs=c(0.025,0.975))
sum(bcon_tile_716>0)/4000 #73% prob

#Surface and Tile
bsur_tile_716<-((exp(tileb_716)*2.78)/4)-((exp(surb_716)*2.78)/4)
mean(bsur_tile_716)
quantile(bsur_tile_716, probs=c(0.025,0.975))
sum(bsur_tile_716>0)/4000 #68% prob

exp(mean(conb_716))*2.78/4 #27.2
exp(mean(surb_716))*2.78/4 #11.1
exp(mean(tileb_716))*2.78/4 #36.8
