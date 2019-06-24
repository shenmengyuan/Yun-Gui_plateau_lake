## =======================================================================================
##
##   The Impacts of Different Trophic Status on Community Structure 
##	 and Metabolic Potential of Microbiota in the Lakes of Yun-Gui plateau of China
##
##   *Figure 1. Structure and composition of the microbial communities along with the different trophic types
##
##  |  2019-05-09
## 
##   Shen Mengyuan, mengyuanshen@126.com
##
## =======================================================================================


# Loading R packages --------------------------------------------------------------------
library(dendextend)
library(vegan)
library(tidyverse)
library(reshape2)    
library(cowplot)
library(RColorBrewer)
library(ggsci)
library(dendsort)
library(xlsx)

# Loading data ----------------------------------------------------------------
## kaiju phylum
kaiju <- read.csv(file = "../Datasets/kaiju_all_lake_phylum_count.txt",sep = "\t",stringsAsFactors = F,header = F)
colnames(kaiju) <- c("lake","count","phylum")
kaiju <- spread(kaiju,lake,count) %>% column_to_rownames(var="phylum") 
kaiju[is.na(kaiju)] <- 0
colnames(kaiju) <- c("DCL-1","DCL-2","XYL-1","XYL-2","EL-1","EL-2","FXL-1",
                     "FXL-2","LGL-1","LGL-2")

# Phyla with relative abundance not in the top ten are shown as “Other”.
taxonomy_counts <- kaiju
for (i in 1:ncol(taxonomy_counts)) {
  taxonomy_counts[,i] <- taxonomy_counts[,i]/sum(taxonomy_counts[,i])
}

norm_tax_df = as.data.frame(taxonomy_counts)
norm_tax_df$abundance = rowSums(norm_tax_df)/ncol(norm_tax_df)
norm_tax_df <- norm_tax_df[order(norm_tax_df$abundance, decreasing = TRUE), ]
phylum_top10 <- norm_tax_df[1:10, -ncol(norm_tax_df)]
phylum_top10['Others', ] <- 1 - colSums(phylum_top10)


phylum_top10 <- select(phylum_top10,
                                 c("EL-1","EL-2","FXL-2",
                                   "FXL-1","LGL-2","LGL-1",
                                   "DCL-1","DCL-2","XYL-1",
                                   "XYL-2"))

# Barplot -----------------------------------------------------------------
phylum_top10_new <- phylum_top10
phylum_top10_new$Taxonomy <- factor(rownames(phylum_top10), levels = rev(rownames(phylum_top10)))
phylum_top10_new <- melt(phylum_top10_new, id = 'Taxonomy')

p <- ggplot(phylum_top10_new, aes(variable, 100 * value, fill = Taxonomy)) +
  geom_col(position = 'stack', width = 0.6,color="black",size=0.1) +
  labs(x = '', y = 'Relative Abundance(%)') +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 13)) +
  theme(legend.text = element_text(size = 11)) +
  coord_flip() +
  scale_fill_manual(values =  rev(c("#D0D1E6","#A6BDDB",
    "#67A9CF","#3690C0","#02818A","#016C59",
    "#014636","#74C476","#A1D99B","#C7E9C0","#E5F5E0","#F7FCF5")))+
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  theme(legend.title = element_blank()) + 
  theme_bw()+theme(panel.grid=element_blank())
p
ggsave(filename = "../图表/raw_20190509/001_累积柱状图_two_color_add.pdf",p,
       width = 6,height = 4)

# Hierarchical clustering -----------------------------------------------------------------
data <- kaiju %>% t() %>% decostand(method = "hellinger")
col_dend = as.dendrogram(dendsort(hclust(vegdist(data,method = "bray"), method ="average" )))
col_dend = color_branches(col_dend, k = 2,col=c("#F781BF","#A65628"))
p <- dendlist(col_dend,col_dend)
tanglegram(p,sort = TRUE,lwd = 1.5)

# PcoA -----------------------------------------------------------------
kaiju <- t(kaiju) %>% decostand(method = "hellinger")
m = "bray_curtis_Phylum"
distance <- vegdist(kaiju, method = 'bray')
pcoa <- cmdscale(distance, k = 4, eig = T)
points = as.data.frame(pcoa$points) 
eig = pcoa$eig
# rename group name
colnames(points) = c("PC1", "PC2", "PC3", "PC4")
points$group<-rep(c("A","B"),c(4,6))
points$group_2 <- rep(c("Eutrophic", "Mesoeutrophic",
                        "Oligomesotrophic", "Oligotrophic"), c(4,2,2,2))
cols <- c("Eutrophic"="#E41A1C", "Mesoeutrophic"="#377EB8",
          "Oligomesotrophic"="#4DAF4A","Oligotrophic"="#984EA3",
          "A"="#A65628","B"="#F781BF")
p_phylum <- ggplot(points, aes(x=PC1, y=PC2, color=group)) + 
  geom_point(alpha=.7, size=2,show.legend=F) +
  labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")
  #     title=paste(m," PCoA",sep="")
  ) +theme_classic()+
  stat_ellipse(level=0.75,show.legend=F,linetype=3,size=0.8) + 
  theme_bw()+theme(panel.grid=element_blank()) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  geom_text_repel(aes(PC1,PC2,label=row.names(points),color=group_2),show.legend=F)
p_phylum

