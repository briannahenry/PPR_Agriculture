library(brms)
library(ggplot2)

write.csv(ysi, file="ysi2.csv")

ggplot(ysi2, aes(x=trt, y=do, color=site))+
  geom_point()+
  theme_classic()

ggplot(ysi2, aes(x=trt, y=temp, color=site))+
  geom_point()+
  theme_classic()

ggplot(ysi2, aes(x=trt, y=cond, color=site))+
  geom_point()+
  theme_classic()

ggplot(ysi2, aes(x=trt, y=ph, color=site))+
  geom_point()+
  theme_classic()

ggplot(ysi2, aes(x=trt, y=turb, color=site))+
  geom_point()+
  theme_classic()
