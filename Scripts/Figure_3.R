## =======================================================================================
##
##   The Impacts of Different Trophic Status on Community Structure 
##	 and Metabolic Potential of Planktonic Microbiota in Lakes on Yun-Gui plateau of China
##
##   *Figure 3. Environmental drivers of microbial community composition.
##
##  |  2019-05-07
## 
##   Shen Mengyuan, mengyuanshen@126.com
##
## =======================================================================================

# Loading R packages --------------------------------------------------------------------
library(tidyverse)
library(vegan)
library(ggrepel)
library(ggord)
library(cowplot)

# Loading data ----------------------------------------------------------------
## Taxonomic composition
kaiju <- read.csv(file = "../Datasets/kaiju_all_lake_phylum_count.txt",sep = "\t",stringsAsFactors = F,header = F)
colnames(kaiju) <- c("lake","count","phylum")
kaiju <- spread(kaiju,lake,count) %>% column_to_rownames(var="phylum") 
kaiju[is.na(kaiju)] <- 0
colnames(kaiju) <- c("DCL-1","DCL-2","XYL-1","XYL-2","EL-1","EL-2","FXL-1",
                     "FXL-2","LGL-1","LGL-2")

kaiju <- t(kaiju) %>% decostand(method = "hellinger")

## Environmental factors
env <- read.csv(file = "../Datasets/env.txt",stringsAsFactors = F,sep="\t") %>%
  column_to_rownames(var="id") 
env.row <- select(env, c(PH,WT,area,volume,water.level))
env.log2 <- log2(select(env,c(TP,TN,average.water.depth)))
env.norm <- cbind(env.row,env.log2)
env.z <- decostand(env.norm, "standardize")
# RDA ---------------------------------------------------------------------
## Only the significant environmental variables 
## were determined (p < 0.05) for RDA
WT.rda<-rda(kaiju ~ env.z$WT);anova(WT.rda, step =1000)
# Pr(>F):0.045*
PH.rda<-rda(kaiju ~ env.z$PH);anova(PH.rda, step =1000)
# Pr(>F):0.316
TP.rda<-rda(kaiju ~ env.z$TP);anova(TP.rda, step =1000)
# Pr(>F):0.001***
TN.rda<-rda(kaiju ~ env.z$TN);anova(TN.rda, step =1000)
# Pr(>F):0.003**
area.rda<-rda(kaiju ~ env.z$area);anova(area.rda, step =1000)
# Pr(>F):0.821
average.water.depth.rda<-rda(kaiju ~ env.norm$average.water.depth);anova(average.water.depth.rda, step =1000)
# Pr(>F):0.007**
volume.rda<-rda(kaiju ~ env.norm$volume);anova(volume.rda, step =1000)
# Pr(>F):0.137
water.level.rda<-rda(kaiju ~ env.norm$water.level);anova(water.level.rda, step =1000)
# Pr(>F):0.232

otu.tab.1<- rda(kaiju ~ WT+TN+TP+average.water.depth, env.z)
## to avoid the col-linearity among environmental variables,
## high variance inflation factors (VIF >20) were eliminated
vif.cca(otu.tab.1)
summary(otu.tab.1,scaling = 2)

sp <-  as.data.frame(scores(otu.tab.1, choices = 1:2, scaling = 2, display = 'sp'))
sp$LEN <- (sp$RDA1)^2+(sp$RDA2)^2
new_sp <- sp[order(sp[,3],decreasing=T),][1:8,]
st <-  as.data.frame(scores(otu.tab.1, choices = 1:2, scaling = 2, display = 'wa'))
env <- as.data.frame(scores(otu.tab.1, choices = 1:2, scaling = 2, display = 'bp'))

group_1<-rep(c("A","B"),c(4,6))
group_2 <- rep(c("Eutrophic", "Mesoeutrophic",
                 "Oligomesotrophic", "Oligotrophic"), c(4,2,2,2))
st$group <- group_1
cols <- c("Eutrophic"="#E41A1C", "Mesoeutrophic"="#377EB8",
          "Oligomesotrophic"="#4DAF4A","Oligotrophic"="#984EA3",
          "A"="#A65628","B"="#F781BF")

p_env <- ggplot(st, aes(x=RDA1, y=RDA2, color=group)) + 
  geom_point(alpha=.7, size=2,show.legend=F) +
  theme_classic()+
  stat_ellipse(level=0.75,show.legend=F,linetype=3,size=0.75) +
  theme_bw()+theme(panel.grid=element_blank()) +
  geom_hline(yintercept=0,linetype=3,size=1, colour = "darkgrey") + 
  geom_vline(xintercept=0,linetype=3,size=1, colour = "darkgrey") +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  labs(x="RDA 1 (79.65%)", sep="",y="RDA 1 (5.66%)", sep="") +
  geom_text_repel(data = st, aes(RDA1,RDA2,label=row.names(st),color=group_2),show.legend=F) +
  geom_text_repel(data = new_sp,aes(x=RDA1,y=RDA2,label=row.names(new_sp)),colour = "black") +
  geom_text_repel(data = env,aes(x=RDA1,y=RDA2,label=row.names(env)),colour = "black") +
  geom_point(data = new_sp, aes(x=RDA1, y=RDA2), shape= 2, size = 2,colour = "black") +
  geom_segment(data = env,aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
               arrow = arrow(angle=30,length = unit(0.2,"cm"),
                             type = "open"),linetype=1, size=0.5,colour = "black")

p_env

ggsave(file = "Phylum_env_RDA_ellipse_0.75.pdf",p_env,width=5, height=5)


# variation partitioning --------------------------------------------------
anova(otu.tab.1, step = 1000)
spe.env0.rda <- rda(kaiju ~ 1, env.z )
step.forward <-ordistep(spe.env0.rda,scope=formula(otu.tab.1),
                        direction='forward', pstep=1000)
step.forward 
re_vp <- varpart(kaiju, ~TN, ~TP, data=env.z)
plot(re_vp,digits =2,Xnames = c('TN','TP'), bg = c('blue', 'red'))

