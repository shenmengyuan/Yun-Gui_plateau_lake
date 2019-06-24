## =======================================================================================
##
##   The Impacts of Different Trophic Status on Community Structure 
##	 and Metabolic Potential of Microbiota in the Lakes of Yun-Gui plateau of China
##
##   *Figure 3. Environmental drivers of microbial community composition
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
kaiju <- read.csv(file = "../Datasets/kaiju_all_lake_phylum_count.txt",sep = "\t",stringsAsFactors = F,header = F)
colnames(kaiju) <- c("lake","count","phylum")
kaiju <- spread(kaiju,lake,count) %>% column_to_rownames(var="phylum") 
kaiju[is.na(kaiju)] <- 0
colnames(kaiju) <- c("DCL-1","DCL-2","XYL-1","XYL-2","EL-1","EL-2","FXL-1",
                     "FXL-2","LGL-1","LGL-2")
					 
kaiju <- t(kaiju) %>% decostand(method = "hellinger")

env <- read.csv(file = "../Datasets/env.txt",stringsAsFactors = F,sep="\t") %>%
  column_to_rownames(var="id") 
env.z <- decostand(env, "standardize")

# RDA
env <- env.z
otu.tab.1<- rda(kaiju ~ T+TN+TP+average.water.depth, env)

# otu.tab.1<- rda(kaiju ~ T+TP+PH+average.water.depth+water.level, env)
# vif.cca(otu.tab.1)

sp=data.frame(otu.tab.1$CCA$v[, 1:2])
sp$LEN <- (sp$RDA1)^2+(sp$RDA2)^2
new_sp <- sp[order(sp[,3],decreasing=T),][1:5,]
st=data.frame(otu.tab.1$CCA$wa[, 1:2])
yz=data.frame(otu.tab.1$CCA$biplot[, 1:2])

rownames(st) <- c("DCL-1","DCL-2","XYL-1","XYL-2","EL-1","EL-2","FXL-1",
                  "FXL-2","LGL-1","LGL-2")
group<-rep(c("A","B"),c(4,6))
group_2 <- rep(c("Eutrophic", "Mesoeutrophic",
                 "Oligomesotrophic", "Oligotrophic"), c(4,2,2,2))
cols <- c("Eutrophic"="#E41A1C", "Mesoeutrophic"="#377EB8",
          "Oligomesotrophic"="#4DAF4A","Oligotrophic"="#984EA3",
          "A"="#F781BF","B"="#A65628")

p_env <- ggord(otu.tab.1,repel=T,grp_in = group,parse = TRUE,
      ellipse=T,size=2,addcol='white',ellipse_pro = 0.80) +
  theme_bw()+theme(panel.grid=element_blank()) +
  geom_hline(yintercept=0,linetype=3,size=1) + 
  geom_vline(xintercept=0,linetype=3,size=1) +
  geom_text_repel(data = st,aes(x=RDA1,y=RDA2,label=row.names(st),color=group_2)) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  geom_segment(data = new_sp,aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
               arrow = arrow(angle=30,length = unit(0.4,"cm"),
                             type = "open"),linetype=1, size=0.5,colour = "black")+
  geom_text_repel(data = new_sp,aes(x=RDA1,y=RDA2,label=row.names(new_sp)),colour = "black")

p_env