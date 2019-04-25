library(brms)
library(ggplot2)

#####Selenium in Water#####
se_water$mo_yr2<-as.factor(se_water$mo_yr)
se_water$trt0f<-as.factor(se_water$trt)

get_prior(se_01~trt0f*mo_yr2 + (1|site/rep), data=se_water, family=Gamma(link="log"))


sewaterm<-brm(se_01~trt0f*mo_yr2 + (1|site/rep), data=se_water, family=Gamma(link="log"),
              prior=c(prior(normal(0,4),class="b"),
                      prior(normal(0,2),class="Intercept"),
                      prior(cauchy(0,1), class="sd")),
              chains=4,iter=2000)

pp_check(sewaterm, type="boxplot")  #unsure what in the priors is allowing occasional insane samples; there are only two samples for june of 2015...
print(sewaterm)
marginal_effects(sewaterm)

seplot<-marginal_effects(sewaterm, robust=FALSE)
seplot<-as.data.frame(seplot$`trt0f:mo_yr2`)
seplot$mo_yr3<-factor(seplot$mo_yr2, levels=c("15-May", "15-Jun", "15-Jul", "16-May", "16-Jun", "16-Jul"))
seplot$yr<-ifelse(seplot$mo_yr2=="15-May", "2015",
                  ifelse(seplot$mo_yr2=="15-Jun", "2015",
                         ifelse(seplot$mo_yr2=="15-Jul", "2015", "2016")))
seplot$mth<-ifelse(seplot$mo_yr2=="15-May", "May",
                  ifelse(seplot$mo_yr2=="15-Jun", "June",
                         ifelse(seplot$mo_yr2=="16-Jun", "June",
                            ifelse(seplot$mo_yr2=="16-May", "May", "July"))))
seplot$mth2<-factor(seplot$mth, levels=c("May","June","July"))

senotile<-ggplot(data=subset(seplot, estimate__<20), aes(x=mth2, y=estimate__, ymin=lower__, ymax=upper__, color=trt0f))+
  geom_point(size=3, position=position_dodge(width=0.4))+
  geom_errorbar(position=position_dodge(width=0.4))+
  geom_hline(yintercept = 2, color="slategrey")+
  geom_hline(yintercept = 5, color="black", linetype="dashed")+
  scale_color_manual(name="Treatment", labels = c("Control", "Surface", "Tile"),  values = c("grey0","dodgerblue4","tomato4"))+
  facet_wrap(~yr)+
  theme_classic()+
  xlab("Date Sampled")+
  ylab(expression(paste("Dissolved Selenium (", mu,"g/l)")))+
  theme(axis.text.x = element_text(color="black"))+
  theme(axis.text.y = element_text(color="black"))+
  theme(legend.title = element_text(face="bold", size="11"))

senotile

ggsave("seleniu_nooutflow.tiff", senotile, dpi=400, width=7, height=4.5, units="in")

seall<-ggplot(seplot, aes(x=mth2, y=estimate__, ymin=lower__, ymax=upper__, color=trt0f))+
  geom_point(size=3, position=position_dodge(width=0.4))+
  geom_errorbar(position=position_dodge(width=0.4))+
  scale_color_manual(name="Treatment", labels = c("Control", "Surface", "Tile", "Outflow"),  values = c("grey0","dodgerblue4","tomato4", "goldenrod4"))+
  facet_wrap(~yr)+
  theme_classic()+
  xlab("Date Sampled")+
  ylab(expression(paste("Dissolved Selenium (", mu,"g/l)")))+
  theme(axis.text.x = element_text(color="black"))+
  theme(axis.text.y = element_text(color="black"))+
  scale_y_log10()+
  theme(legend.title = element_text(face="bold", size="11"))

seall

ggsave("selenium_alltrt.tiff", seall, dpi=800)


#########Differences by trt and month#######
sepost<-posterior_samples(sewaterm)
str(sepost)

#June 2015
scon_615<-sepost$b_Intercept+sepost$b_mo_yr215MJun
ssur_615<-sepost$b_Intercept+sepost$b_mo_yr215MJun+sepost$`b_trt0fSurface:mo_yr215MJun`+sepost$b_trt0fSurface
stile_615<-sepost$b_Intercept+sepost$b_mo_yr215MJun+sepost$`b_trt0fTile:mo_yr215MJun`+sepost$b_trt0fTile

#Control and Surface
secon_sur_615<-exp(ssur_615)-exp(scon_615)
mean(secon_sur_615)
quantile(secon_sur_615, probs=c(0.025,0.975))
sum(secon_sur_615>0)/4000 #90.6% prob

#Surface and Tile
setile_sur_615<-exp(ssur_615)-exp(stile_615)
mean(setile_sur_615)
quantile(setile_sur_615, probs=c(0.025,0.975))
sum(setile_sur_615>0)/4000 #60.6 prob

#Control and Tile
secon_tile_615<-exp(stile_615)-exp(scon_615)
mean(secon_tile_615)
quantile(secon_tile_615, probs=c(0.025,0.975))
sum(secon_tile_615>0)/4000 #85.98% prob

#Differences
exp(mean(scon_615))#0.56
exp(mean(ssur_615))#1.66
exp(mean(stile_615))#1.33
#Surface selenium is over 3 times higher than that of control sites

#May 2015
scon_515<-sepost$b_Intercept+sepost$b_mo_yr215MMay
ssur_515<-sepost$b_Intercept+sepost$b_mo_yr215MMay+sepost$`b_trt0fSurface:mo_yr215MMay`+sepost$b_trt0fSurface
stile_515<-sepost$b_Intercept+sepost$b_mo_yr215MMay+sepost$`b_trt0fTile:mo_yr215MMay`+sepost$b_trt0fTile

#Control and Surface
secon_sur_515<-exp(ssur_515)-exp(scon_515)
mean(secon_sur_515)
quantile(secon_sur_515, probs=c(0.025,0.975))
sum(secon_sur_515>0)/4000 #90.6% prob

#Surface and Tile
setile_sur_515<-exp(ssur_515)-exp(stile_515)
mean(setile_sur_515)
quantile(setile_sur_515, probs=c(0.025,0.975))
sum(setile_sur_515>0)/4000 #50.5 prob

#Control and Tile
secon_tile_515<-exp(stile_515)-exp(scon_515)
mean(secon_tile_515)
quantile(secon_tile_515, probs=c(0.025,0.975))
sum(secon_tile_515>0)/4000 #85.98% prob

#July 2015
scon_715<-sepost$b_Intercept
ssur_715<-sepost$b_Intercept+sepost$b_trt0fSurface
stile_715<-sepost$b_Intercept+sepost$b_trt0fTile
exp(mean(scon_715))
exp(mean(ssur_715))
exp(mean(stile_715))

#Control and Surface
secon_sur_715<-exp(ssur_715)-exp(scon_715)
mean(secon_sur_715)
quantile(secon_sur_715, probs=c(0.025,0.975))
sum(secon_sur_715>0)/4000 #90.6% prob

#Surface and Tile
setile_sur_715<-exp(ssur_715)-exp(stile_715)
mean(setile_sur_715)
quantile(setile_sur_715, probs=c(0.025,0.975))
sum(setile_sur_715>0)/4000 #50.5 prob

#Control and Tile
secon_tile_715<-exp(stile_715)-exp(scon_715)
mean(secon_tile_715)
quantile(secon_tile_715, probs=c(0.025,0.975))
sum(secon_tile_715>0)/4000 #85.98% prob





#June 2016
scon_616<-sepost$b_Intercept+sepost$b_mo_yr216MJun
ssur_616<-sepost$b_Intercept+sepost$b_mo_yr216MJun+sepost$`b_trt0fSurface:mo_yr216MJun`+sepost$b_trt0fSurface
stile_616<-sepost$b_Intercept+sepost$b_mo_yr216MJun+sepost$`b_trt0fTile:mo_yr216MJun`+sepost$b_trt0fTile

#Control and Surface
secon_sur_616<-exp(ssur_616)-exp(scon_616)
mean(secon_sur_616)
quantile(secon_sur_616, probs=c(0.025,0.975))
sum(secon_sur_616>0)/4000 #91% prob

#Surface and Tile
setile_sur_616<-exp(stile_616)-exp(ssur_616)
mean(setile_sur_616)
quantile(setile_sur_616, probs=c(0.025,0.975))
sum(setile_sur_616>0)/4000 #99.8% prob

#Control and Tile
secon_tile_616<-exp(stile_616)-exp(scon_616)
mean(secon_tile_616)
quantile(secon_tile_616, probs=c(0.025,0.975))
sum(secon_tile_616>0)/4000 #>99.9% prob

#Percent diff
exp(mean(stile_616))#5.55
exp(mean(scon_616))#0.17
exp(mean(ssur_616))#0.47
#Surface is more than 2.5 times higher than control
#Tile is over 32 times higher than control treatment
#Tile is 11.8 times higher than surface 

#May 2016
scon_516<-sepost$b_Intercept+sepost$b_mo_yr216MMay
ssur_516<-sepost$b_Intercept+sepost$b_mo_yr216MMay+sepost$`b_trt0fSurface:mo_yr216MMay`+sepost$b_trt0fSurface
stile_516<-sepost$b_Intercept+sepost$b_mo_yr216MMay+sepost$`b_trt0fTile:mo_yr216MMay`+sepost$b_trt0fTile

#Control and Surface
secon_sur_516<-exp(ssur_516)-exp(scon_516)
mean(secon_sur_516)
quantile(secon_sur_516, probs=c(0.025,0.975))
sum(secon_sur_516>0)/4000 #86.9% prob

#Surface and Tile
setile_sur_516<-exp(stile_516)-exp(ssur_516)
mean(setile_sur_516)
quantile(setile_sur_516, probs=c(0.025,0.975))
sum(setile_sur_516>0)/4000 #95.4% prob

#Control and Tile
bsecon_tile_516<-exp(stile_516)-exp(scon_516)
mean(secon_tile_516)
quantile(secon_tile_516, probs=c(0.025,0.975))
sum(secon_tile_516>0)/4000 #96.9% prob

#Percent diff
exp(mean(stile_516))#3.17
exp(mean(scon_516))#0.37
exp(mean(ssur_516))#0.89
#Tile is 8.5 times higher than control treatment
#Tile is 3.5 times higher than surface 

#July 2016
scon_716<-sepost$b_Intercept+sepost$b_mo_yr216MJul
ssur_716<-sepost$b_Intercept+sepost$b_mo_yr216MJul+sepost$`b_trt0fSurface:mo_yr216MJul`+sepost$b_trt0fSurface
stile_716<-sepost$b_Intercept+sepost$b_mo_yr216MJul+sepost$`b_trt0fTile:mo_yr216MJul`+sepost$b_trt0fTile

#Control and Surface
secon_sur_716<-exp(ssur_716)-exp(scon_716)
mean(secon_sur_716)
quantile(secon_sur_716, probs=c(0.025,0.975))
sum(secon_sur_716>0)/4000 #85.1% prob

#Surface and Tile
setile_sur_716<-exp(stile_716)-exp(ssur_716)
mean(setile_sur_716)
quantile(setile_sur_716, probs=c(0.025,0.975))
sum(setile_sur_716>0)/4000 #93.6% prob

#Control and Tile
secon_tile_716<-exp(stile_716)-exp(scon_716)
mean(secon_tile_716)
quantile(secon_tile_716, probs=c(0.025,0.975))
sum(secon_tile_716>0)/4000 #99.5% prob

#percent diff
exp(mean(stile_716))#1.88
exp(mean(scon_716))#0.265
exp(mean(ssur_716))#0.588
#Tile is 7 times higher than the control and 3 times higher than surface

#####clothianidin####

neonic$mo_yr<-as.factor(neonic$mth_yr)
neonic$trt0f<-as.factor(neonic$trt)

#clothianidin
clothm<-brm(cloth0001~trt0f*mo_yr + (1|SITE), data=neonic, family=Gamma(link="log"),
              prior=c(prior(normal(0,2),class="b"),
                      prior(normal(0,2),class="Intercept"),
                      prior(cauchy(0,1), class="sd")),
              chains=4,iter=2000)

pp_check(clothm, type="boxplot")
print(clothm)
marginal_effects(clothm, robust=FALSE)

clothplot<-marginal_effects(clothm)
clothplot<-as.data.frame(clothplot$`trt0f:mo_yr`)
clothplot$yr<-ifelse(clothplot$mo_yr=="May_2015", "2015",
                  ifelse(clothplot$mo_yr=="June_2015", "2015",
                         ifelse(clothplot$mo_yr=="July_2015", "2015", "2016")))
clothplot$mth<-ifelse(clothplot$mo_yr=="May_2015", "May",
                   ifelse(clothplot$mo_yr=="June_2015", "June",
                          ifelse(clothplot$mo_yr=="June_2016", "June",
                                 ifelse(clothplot$mo_yr=="May_2016", "May", "July"))))
clothplot$mth2<-factor(clothplot$mth, levels=c("May","June","July"))
clothplot$trt2<-factor(clothplot$trt0f, levels=c("Control", "Surface", "Tile", "Outflow"))
clothplot$moyr2<-ifelse(clothplot$mo_yr=="May_2015", "May 2015",
                        ifelse(clothplot$mo_yr=="June_2015", "June 2015",
                               ifelse(clothplot$mo_yr=="July_2015", "July 2015", "May 2016")))
clothplot$moyr3<-factor(clothplot$moyr2, levels=c("May 2015", "June 2015", "July 2015", "May 2016"))

clothp<-ggplot(clothplot, aes(x=moyr3, y=estimate__, ymin=lower__, ymax=upper__, color=trt2))+
  geom_point(size=3.5, position=position_dodge(width=0.4))+
  geom_errorbar(position=position_dodge(width=0.4))+
  scale_color_manual(name="Treatment", labels = c("Control", "Surface", "Tile", "Outflow"),  values = c("grey0","dodgerblue4","tomato4", "goldenrod4"))+
  theme_classic()+
  xlab("Date Sampled")+
  ylab(expression(paste("Clothianidin (", mu,"g/l)")))+
  theme(axis.text.x = element_text(color="black"))+
  theme(axis.text.y = element_text(color="black"))+
  theme(legend.title = element_text(face="bold", size="11"))

clothp

ggsave("Clothianidin.tiff", clothp, dpi=400, width=7, height=4.5, units="in")

######Differences Trt and date#######
clopost<-posterior_samples(clothm)
str(clopost)

#May 2015
cc_515<-clopost$b_Intercept+clopost$b_mo_yrMay_2015
cs_515<-clopost$b_Intercept+clopost$b_mo_yrMay_2015+clopost$`b_trt0fSurface:mo_yrMay_2015`+clopost$b_trt0fSurface
ct_515<-clopost$b_Intercept+clopost$b_mo_yrMay_2015+clopost$`b_trt0fTile:mo_yrMay_2015`+clopost$b_trt0fTile

#Control and Surface
cc_cs_515<-exp(cs_515)-exp(cc_515)
mean(cc_cs_515)
quantile(cc_cs_515, probs=c(0.025,0.975))
sum(cc_cs_515>0)/4000 #95% prob

#Surface and Tile
ct_cs_515<-exp(cs_515)-exp(ct_515)
mean(ct_cs_515)
quantile(ct_cs_515, probs=c(0.025,0.975))
sum(ct_cs_515>0)/4000 #95.5% prob

#Control and Tile
cc_ct_515<-exp(cc_515)-exp(ct_515)
mean(cc_ct_515)
quantile(cc_ct_515, probs=c(0.025,0.975))
sum(cc_ct_515>0)/4000 #56.85% prob

#Differences
exp(mean(cc_515))#0.00066
exp(mean(cs_515))#0.0046
exp(mean(ct_515))#0.00055
#Surface treatment is almost 8 times higher than the control treatment,
#and 8.4 times higher than the tile drainage

#June 2015
cc_615<-clopost$b_Intercept+clopost$b_mo_yrJune_2015
cs_615<-clopost$b_Intercept+clopost$b_mo_yrJune_2015+clopost$`b_trt0fSurface:mo_yrJune_2015`+clopost$b_trt0fSurface
ct_615<-clopost$b_Intercept+clopost$b_mo_yrJune_2015+clopost$`b_trt0fTile:mo_yrJune_2015`+clopost$b_trt0fTile
co_615<-clopost$b_Intercept+clopost$b_mo_yrJune_2015+clopost$`b_trt0fOutflow:mo_yrJune_2015`+clopost$b_trt0fOutflow

#Control and Surface
cc_cs_615<-exp(cc_615)-exp(cs_615)
mean(cc_cs_615)
quantile(cc_cs_615, probs=c(0.025,0.975))
sum(cc_cs_615>0)/4000 #50.8% prob

#Surface and Tile
ct_cs_615<-exp(ct_615)-exp(cs_615)
mean(ct_cs_615)
quantile(ct_cs_615, probs=c(0.025,0.975))
sum(ct_cs_615>0)/4000 #48.78% prob

#Control and Tile
cc_ct_615<-exp(ct_615)-exp(cc_615)
mean(cc_ct_615)
quantile(cc_ct_615, probs=c(0.025,0.975))
sum(cc_ct_615>0)/4000 #48.1% prob

#Control and OUtflow
cc_co_615<-exp(cc_615)-exp(co_615)
mean(cc_co_615)
quantile(cc_co_615, probs=c(0.025,0.975))
sum(cc_co_615>0)/4000 #42.6%

#Surface and Outflow
cs_co_615<-exp(co_615)-exp(cs_615)
mean(cs_co_615)
quantile(cs_co_615, probs=c(0.025,0.975))
sum(cs_co_615>0)/4000 #57.5%

#tile and outflow
ct_co_615<-exp(ct_615)-exp(co_615)
mean(ct_co_615)
quantile(ct_co_615, probs=c(0.025,0.975))
sum(ct_co_615>0)/4000 #41.4%

exp(mean(cc_615))#0.000329
exp(mean(cs_615))#0.000322
exp(mean(ct_615))#0.000306
exp(mean(co_615))#0.000502

#ALL THE SAME 

#July 2015
cc_715<-clopost$b_Intercept
cs_715<-clopost$b_Intercept+clopost$b_trt0fSurface
ct_715<-clopost$b_Intercept+clopost$b_trt0fTile
co_715<-clopost$b_Intercept+clopost$b_trt0fOutflow

#Control and Surface
cc_cs_715<-exp(cs_715)-exp(cc_715)
mean(cc_cs_715)
quantile(cc_cs_715, probs=c(0.025,0.975))
sum(cc_cs_715>0)/4000 #91.8% prob

#Surface and Tile
ct_cs_715<-exp(cs_715)-exp(ct_715)
mean(ct_cs_715)
quantile(ct_cs_715, probs=c(0.025,0.975))
sum(ct_cs_715>0)/4000 #91.4% prob

#Control and Tile
cc_ct_715<-exp(cc_715)-exp(ct_715)
mean(cc_ct_715)
quantile(cc_ct_715, probs=c(0.025,0.975))
sum(cc_ct_715>0)/4000 #57% prob

#Control and OUtflow
cc_co_715<-exp(co_715)-exp(cc_715)
mean(cc_co_715)
quantile(cc_co_715, probs=c(0.025,0.975))
sum(cc_co_715>0)/4000 #93.3%

#Surface and Outflow
cs_co_715<-exp(co_715)-exp(cs_715)
mean(cs_co_715)
quantile(cs_co_715, probs=c(0.025,0.975))
sum(cs_co_715>0)/4000 #65%

#tile and outflow
ct_co_715<-exp(co_715)-exp(ct_715)
mean(ct_co_715)
quantile(ct_co_715, probs=c(0.025,0.975))
sum(ct_co_715>0)/4000 #92.1%

(exp(mean(co_715)))#0.0047
(exp(mean(cs_715)))#0.0027
(exp(mean(ct_715)))#0.00055
(exp(mean(cc_715)))#0.00065
#Surface is 4 times higher than control and almost 5 times higher than tile drainage
#Outflow is comparable to surface but 7 times higher than control and 8.5 times higher than tile drainage

#May 2016

cc_516<-clopost$b_Intercept+clopost$b_mo_yrMay_2016
cs_516<-clopost$b_Intercept+clopost$b_mo_yrMay_2016+clopost$`b_trt0fSurface:mo_yrMay_2016`+clopost$b_trt0fSurface
ct_516<-clopost$b_Intercept+clopost$b_mo_yrMay_2016+clopost$`b_trt0fTile:mo_yrMay_2016`+clopost$b_trt0fTile
co_516<-clopost$b_Intercept+clopost$b_mo_yrMay_2016+clopost$`b_trt0fOutflow:mo_yrMay_2016`+clopost$b_trt0fOutflow

#Diff outflow and tile
co_ct516<-exp(co_516)-exp(ct_516)
mean(co_ct516)
quantile(co_ct516, probs=c(0.025,0.975))
sum(co_ct516>0)/4000 #80%

#########Imidacloprid#######
imidam<-brm(imid0001~trt0f*mo_yr + (1|SITE), data=neonic, family=Gamma(link="log"),
            prior=c(prior(normal(0,2),class="b"),
                    prior(normal(0,2),class="Intercept"),
                    prior(cauchy(0,1), class="sd")),
            chains=4,iter=2000)


pp_check(imidam, type="boxplot")
print(imidam)
marginal_effects(imidam, robust=FALSE)


imidaplot<-marginal_effects(imidam)
imidaplot<-as.data.frame(imidaplot$`trt0f:mo_yr`)
imidaplot$yr<-ifelse(imidaplot$mo_yr=="May_2015", "2015",
                    ifelse(imidaplot$mo_yr=="June_2015", "2015",
                           ifelse(imidaplot$mo_yr=="July_2015", "2015", "2016")))
imidaplot$mth<-ifelse(imidaplot$mo_yr=="May_2015", "May",
                      ifelse(imidaplot$mo_yr=="June_2015", "June",
                             ifelse(imidaplot$mo_yr=="June_2016", "June",
                                    ifelse(imidaplot$mo_yr=="May_2016", "May", "July"))))
imidaplot$mth2<-factor(imidaplot$mth, levels=c("May","June","July"))
imidaplot$trt2<-factor(imidaplot$trt0f, levels=c("Control", "Surface", "Tile", "Outflow"))
imidaplot$moyr2<-ifelse(imidaplot$mo_yr=="May_2015", "May 2015",
                        ifelse(imidaplot$mo_yr=="June_2015", "June 2015",
                               ifelse(imidaplot$mo_yr=="July_2015", "July 2015", "May 2016")))
imidaplot$moyr3<-factor(imidaplot$moyr2, levels=c("May 2015", "June 2015", "July 2015", "May 2016"))

imidap<-ggplot(imidaplot, aes(x=moyr3, y=estimate__, ymin=lower__, ymax=upper__, color=trt2))+
  geom_point(size=3, position=position_dodge(width=0.4))+
  geom_errorbar(position=position_dodge(width=0.4))+
  geom_hline(yintercept = 0.01, color="red3", linetype="dashed")+
  scale_color_manual(name="Treatment", labels = c("Control", "Surface", "Tile", "Outflow"),  values = c("grey0","dodgerblue4","tomato4", "goldenrod4"))+
  theme_classic()+
  xlab("Date Sampled")+
  ylab(expression(paste("Imidacloprid (", mu,"g/l)")))+
  theme(axis.text.x = element_text(color="black"))+
  theme(axis.text.y = element_text(color="black"))+
  theme(legend.title = element_text(face="bold", size="11"))

imidap

ggsave("Imidacloprid.tiff", imidap, dpi=400, width=7, height=4.5, units="in")



######Differences by trt and date######
impost<-posterior_samples(imidam)
str(impost)

#June 2015
ic_615<-impost$b_Intercept+impost$b_mo_yrJune_2015
is_615<-impost$b_Intercept+impost$b_mo_yrJune_2015+impost$b_trt0fSurface+impost$`b_trt0fSurface:mo_yrJune_2015`
it_615<-impost$b_Intercept+impost$b_mo_yrJune_2015+impost$b_trt0fTile+impost$`b_trt0fTile:mo_yrJune_2015`
io_615<-impost$b_Intercept+impost$b_mo_yrJune_2015+impost$b_trt0fOutflow+impost$`b_trt0fOutflow:mo_yrJune_2015`

#control and outflow
icio615<-exp(io_615)-exp(ic_615)
mean(icio615)
quantile(icio615, probs=c(0.025,0.975))
sum(icio615>0)/4000 #99.23%

#surface and outflow
isio615<-exp(io_615)-exp(is_615)
mean(isio615)
quantile(isio615, probs=c(0.025,0.975))
sum(isio615>0)/4000#99.275%

#tile and outflow
itio615<-exp(io_615)-exp(it_615)
mean(itio615)
quantile(itio615, probs=c(0.025,0.975))
sum(itio615>0)/4000#99.275%

exp(mean(ic_615))#.00019
exp(mean(is_615))#.000162
exp(mean(it_615))#.000189
exp(mean(io_615))#.0077
#Tile outflow is over 36 times higher than any of the other treatments

#July 2015
ic_715<-impost$b_Intercept
is_715<-impost$b_Intercept+impost$b_trt0fSurface
it_715<-impost$b_Intercept+impost$b_trt0fTile
io_715<-impost$b_Intercept+impost$b_trt0fOutflow

#control adnd tile
icit715<-exp(it_715)-exp(ic_715)
mean(icit715)
quantile(icit715, prob=c(0.025,0.975))
sum(icit715>0)/4000 #94.2%

#surface and tile
isit715<-exp(it_715)-exp(is_715)
mean(isit715)
quantile(isit715, prob=c(0.025,0.975))
sum(isit715>0)/4000 #87.8

#control and outflow
icio715<-exp(io_715)-exp(ic_715)
mean(icio715)
quantile(icio715, probs=c(0.025,0.975))
sum(icio715>0)/4000 #96.9%

#surface and outflow
isio715<-exp(io_715)-exp(is_715)
mean(isio715)
quantile(isio715, probs=c(0.025,0.975))
sum(isio715>0)/4000#95.4%

#tile and outflow
itio715<-exp(io_715)-exp(it_715)
mean(itio715)
quantile(itio715, probs=c(0.025,0.975))
sum(itio715>0)/4000#74.4%

exp(mean(ic_715))#0.00082
exp(mean(is_715))#0.00115
exp(mean(it_715))#0.0035
exp(mean(io_715))#0.0079
#Tile drainage is 4 times higher than control
#Outflow is over 6 times higher than surface and 9.5 times higher than control treatment

#May 2016
ic_516<-impost$b_Intercept+impost$b_mo_yrMay_2016
is_516<-impost$b_Intercept+impost$b_trt0fSurface+impost$b_mo_yrMay_2016+impost$`b_trt0fSurface:mo_yrMay_2016`
it_516<-impost$b_Intercept+impost$b_trt0fTile+impost$b_mo_yrMay_2016+impost$`b_trt0fTile:mo_yrMay_2016`
io_516<-impost$b_Intercept+impost$b_trt0fOutflow+impost$b_mo_yrMay_2016+impost$`b_trt0fOutflow:mo_yrMay_2016`

#control and surface
icis516<-exp(ic_516)-exp(is_516)
mean(icis516)
quantile(icis516, probs=c(0.025,0.975))
sum(icis516>0)/4000#64.8%

#control adnd tile
icit516<-exp(it_516)-exp(ic_516)
mean(icit516)
quantile(icit516, prob=c(0.025,0.975))
sum(icit516>0)/4000 #87.4%

#surface and tile
isit516<-exp(it_516)-exp(is_516)
mean(isit516)
quantile(isit516, prob=c(0.025,0.975))
sum(isit516>0)/4000 #90.35

#control and outflow
icio516<-exp(io_516)-exp(ic_516)
mean(icio516)
quantile(icio516, probs=c(0.025,0.975))
sum(icio516>0)/4000 #87.2

#surface and outflow
isio516<-exp(io_516)-exp(is_516)
mean(isio516)
quantile(isio516, probs=c(0.025,0.975))
sum(isio516>0)/4000#90.5%

#tile and outflow
itio516<-exp(io_516)-exp(it_516)
mean(itio516)
quantile(itio516, probs=c(0.025,0.975))
sum(itio516>0)/4000#60%

exp(mean(ic_516))#0.0146
exp(mean(is_516))#0.0100
exp(mean(it_516))#0.0412
exp(mean(io_516))#0.0569
#Tile drainage is 4 times higher than the surface treatment
#Tile outflow is 5.6 times higher than surface treatment