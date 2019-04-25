library(ggplot2)
library(RColorBrewer)


####brm_model###
library(brms)
##2016 larvae_adults (marg2+sebugs_raw)####
m1<-brm(se_ugg~taxon*stage*trt+(1|site/rep),data=sebbugs_raw,family=Gamma(link="log"),
        prior=c(prior(normal(0,2),class="b"),
                prior(cauchy(0,2),class="sd")))
m1
conditions=data.frame(stage=c("Adult","Larvae"))
plot(marginal_effects(m1,effects='taxon:trt',robust=FALSE,conditions=conditions),points=TRUE,position=position_dodge(width=3))
marg<-marginal_effects(m1,effects='taxon:trt',robust=FALSE,conditions=conditions)
marg2<-as.data.frame(marg$`taxon:trt`)

m1post<-posterior_samples(m1)

#####Chiros
con_cha<-exp(m1post$b_Intercept)
tile_cha<-exp(m1post$b_Intercept+m1post$b_trttile)
surf_cha<-exp(m1post$b_Intercept+m1post$b_trtsurface)

quantile(con_cha,probs=c(0.025,0.975))
quantile(surf_cha,probs=c(0.025,0.975))
quantile(tile_cha,probs=c(0.025,0.975))

sum(tile_cha>3)/4000

foldtc<-tile_cha/con_cha
mean(foldtc)
quantile(foldtc,probs=c(0.025,0.975))

foldts<-tile_cha/surf_cha
mean(foldts)
quantile(foldts,probs=c(0.025,0.975))
sum(foldts>1)/4000


sebugs_raw$estimate__<-sebugs_raw$se_ugg


concern <- data.frame( x = c(-Inf, Inf), y = 3, cutoff = factor(3) )
p + geom_text(data = ann_text,label = "Text" )

p5<-ggplot()+
  geom_pointrange(data=marg2,aes(x=trt,y=estimate__,ymax=upper__,ymin=lower__,color=stage,group=stage),
                            position=position_dodge(width=-.3),size=.6,alpha=.8,shape=16)+
  facet_wrap(~taxon)+
  guides(col = guide_legend(reverse = TRUE))+
  theme_classic()+
  geom_point(data=sebbugs_raw,aes(x=trt,y=se_ugg,color=stage),
                                 position=position_dodge(width=-.3),size=.5,alpha=.7,shape=1)+
  ylab("selenium in insect tissue (ug/g)")+
  theme(text=element_text(size=12))+
  geom_hline(yintercept=c(3,7),linetype='dotted',color="gray",size=.5,show.legend = TRUE)+
  scale_color_manual(values=c("#d8b365","#5ab4ac"))+
  theme(axis.text=element_text(color="black"))+
  scale_y_continuous(minor_breaks=seq(0,23,by=1),breaks=seq(0,20,by=2))
p5
ggsave("p5.pdf",p5,dpi=600,width=5,height=3,units="in")

###2016 larvae_adults plots####

m2p<-ggplot(marg2,aes(x=reorder(month,order),y=estimate__,ymin=lower__,ymax=upper__,color=trt,group=trt))+
  geom_point(position=position_dodge(width=0.35),size=6,alpha=1)+
  geom_linerange(position=position_dodge(width=0.35),size=1.25,alpha=0.7)+
  
  facet_wrap(~year)+
  scale_color_manual(values=c("#696969","#2e75b6","#B62E75"))+
  theme_classic()+
  geom_line(position=position_dodge(width=0.35),size=1.3,alpha=0.8)+
  geom_hline(yintercept=1.2)

m2p 


# Create a new powerpoint document
doc <- pptx()
# Add a new slide into the ppt document 
doc <- addSlide(doc, "Two Content" )
# add a slide title
doc<- addTitle(doc, "" )


# Add an editable box plot
doc <- addPlot(doc, function() print(m12p), vector.graphic = TRUE )
# Add a raster box plot
doc <- addPlot(doc, function() print(m12p), vector.graphic = FALSE )
# write the document to a file
writeDoc(doc, file = "m12p.pptx")

##### adults 1516 (marg4+se_bugsad)####
se_bugsad<-subset(se_bugs1516_2,stage=="adults")
se_bugsad<-subset(se_bugsad,taxon=="Chironomidae"|taxon=="Coenagrionidae")

m4<-brm(se_ugg~taxon*mon_yr*trt+(1|mon_yr/site/rep),data=se_bugsad,family=Gamma(link="log"),
        prior=c(prior(normal(0,2),class="Intercept"),
                prior(normal(0,1),class="b"),
                prior(cauchy(0,1),class="sd")))
m4
pp_check(m4)
conditions=data.frame(mon_yr=c("June_2015","June_2016"))
plot(marginal_effects(m4,effects='taxon:trt',robust=FALSE,conditions=conditions),points=TRUE,position=position_dodge(width=3))
marg4<-marginal_effects(m4,effects='taxon:trt',robust=FALSE,conditions=conditions)
marg4<-as.data.frame(marg4$`taxon:trt`)






#### 1516 overall mean larvae adults
m4b<-brm(se_ugg~taxon*stage*trt+(1|year/site/rep),data=se_bugs1516,family=Gamma(link="log"),
        prior=c(prior(normal(0,2),class="Intercept"),
                prior(normal(0,1),class="b"),
                prior(cauchy(0,1),class="sd")))
m4b
conditions<-data.frame(stage=c("adults","larvae"))
marginal_effects(m4b,effects="taxon:trt",conditions=conditions,robust=FALSE)
m4marg<-marginal_effects(m4b,effects="taxon:trt",conditions=conditions,robust=FALSE)
m4marg<-as.data.frame(m4marg$`taxon:trt`)
colnames(m4marg)

m4p<-ggplot(subset(m4marg,taxon=="Chironomidae"),aes(x=trt,y=estimate__,ymin=lower__,ymax=upper__,color=stage,group=stage))+
  geom_point(position=position_dodge(width=-0.35),size=6,alpha=1)+
  geom_linerange(position=position_dodge(width=-0.35),size=1.25,alpha=0.7)+
  
  #facet_wrap(~year)+
  scale_color_manual(values=c("#696969","#2e75b6"))+
  theme_classic()+
  geom_hline(yintercept=c(3,7),color='red2',size=1)+
  theme(legend.title=element_blank(),
        text=element_text(size=20),
        axis.title.x=element_blank())+
  ylab(expression(paste("Selenium in insects"," (",mu*"g/g dry)")))

m4p 
ggsave("m4p.jpg",m4p,dpi=800)


m4p2<-ggplot(subset(m4marg,taxon=="Chironomidae"),
             aes(x=trt,y=estimate__,ymin=lower__,ymax=upper__,color=stage,group=stage))+
  geom_point(position=position_dodge(width=-0.35),size=6,alpha=1)+
  geom_linerange(position=position_dodge(width=-0.35),size=1.25,alpha=0.7)+
  
  #facet_wrap(~year)+
  scale_color_manual(values=c("#696969","white"))+
  theme_classic()+
  geom_hline(yintercept=c(3,7),color='red2',size=1)+
  theme(legend.title=element_blank(),
        text=element_text(size=20),
        axis.title.x=element_blank())+
  ylab(expression(paste("Selenium in insects"," (",mu*"g/g dry)")))

m4p2 
ggsave("m4p2.jpg",m4p2,dpi=800)


