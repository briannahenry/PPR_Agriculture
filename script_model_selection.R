######selenium in bugs selection#####
install.packages("compare")
library(compare)
library(tidyr)
library(dplyr)

se_bugs1516_2<-merge(se_bugs1516b,sites_only)

write.csv(se_bugs1516,"se_bugs1516.csv")
write.csv(se_bugs1516_2,"se_bugs1516_2.csv")
se_bugs1516_2$min_area


m4b<-brm(se_ugg~taxon*stage*trt+(1|year/site/rep),data=se_bugs1516_2,family=Gamma(link="log"),
         prior=c(prior(normal(0,2),class="Intercept"),
                 prior(normal(0,1),class="b"),
                 prior(cauchy(0,1),class="sd")))


######tot emerge mg####

m11<-brm(totmg~date*trt+(1|site/trap2/taxon),data=emerge2,family=Gamma(link="log"),
         prior=c(prior(normal(0,4),class="b"),
                 prior(normal(1,6),class="Intercept"),
                 prior(cauchy(0,1),class="sd")),
         chains=4,iter=2000)

m11b<-brm(totmg~date*trt*min_area+(1|site/trap2/taxon),data=emerge2,family=Gamma(link="log"),
         prior=c(prior(normal(0,4),class="b"),
                 prior(normal(1,6),class="Intercept"),
                 prior(cauchy(0,1),class="sd")),
         chains=4,iter=2000)

m11c<-brm(totmg~date*trt*depth_min+(1|site/trap2/taxon),data=emerge2,family=Gamma(link="log"),
          prior=c(prior(normal(0,4),class="b"),
                  prior(normal(1,6),class="Intercept"),
                  prior(cauchy(0,1),class="sd")),
          chains=4,iter=2000)

m11d<-brm(totmg~date*min_area+(1|site/trap2/taxon),data=emerge2,family=Gamma(link="log"),
         prior=c(prior(normal(0,4),class="b"),
                 prior(normal(1,6),class="Intercept"),
                 prior(cauchy(0,1),class="sd")),
         chains=4,iter=2000)

m11e<-brm(totmg~date*depth_min+(1|site/trap2/taxon),data=emerge2,family=Gamma(link="log"),
          prior=c(prior(normal(0,4),class="b"),
                  prior(normal(1,6),class="Intercept"),
                  prior(cauchy(0,1),class="sd")),
          chains=4,iter=2000)

emerge2$depth_min_st<-emerge2$depth_min-mean(emerge2$depth_min)/sd(emerge2$depth_min)
m11<-brm(totmg~date*trt+(1|site/trap2/taxon),data=emerge2,family=Gamma(link="log"),
         prior=c(prior(normal(0,4),class="b"),
                 prior(normal(1,6),class="Intercept"),
                 prior(cauchy(0,1),class="sd")),
         chains=4,iter=2000)

WAIC(m11,m11b,m11c,m11d,m11e)

marginal_effects(m11d,effects='min_area')
marginal_effects(m11e,effects='depth_min')

post_m11e<-posterior_samples(m11e)
m11e_sl<-post_m11e$b_depth_min
quantile(m11e_sl,probs=c(0.025,0.975))