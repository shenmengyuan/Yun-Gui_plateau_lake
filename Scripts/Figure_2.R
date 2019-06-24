## =======================================================================================
##
##   The Impacts of Different Trophic Status on Community Structure 
##	 and Metabolic Potential of Microbiota in the Lakes of Yun-Gui plateau of China
##
##   *Figure 2. Metabolism function and composition of the microbial communities along with the different habitat types
##
##  |  2019-05-07
## 
##   Shen Mengyuan, mengyuanshen@126.com
##
## =======================================================================================

# Loading R packages --------------------------------------------------------------------
library(tidyverse)
library(pheatmap)
library(vegan)

# Loading data ----------------------------------------------------------------
kegg_TPM <- read.csv(file = "../Datasets/gene_TPM_name_filter_190506.csv",stringsAsFactors = F,header = T)
kegg_TPM[is.na(kegg_TPM)] <- 0
colnames(kegg_TPM) <- c("Gene_level","DCL-1","DCL-2","XYL-1","XYL-2","EL-1","EL-2",
                        "FXL-1","FXL-2","LGL-1","LGL-2")
kegg_gene_ann <- read.csv(file = "../Datasets/gene.anno.txt",sep = "\t",header = T,stringsAsFactors = T)
kegg_TPM_annotion <- inner_join(kegg_TPM,kegg_gene_ann,by="Gene_level")

Metabolism <- filter(kegg_TPM_annotion,Pathway1=="Metabolism") %>%
              select(c("Gene_level","DCL-1","DCL-2","XYL-1","XYL-2","EL-1","EL-2",
                       "FXL-1","FXL-2","LGL-1","LGL-2")) %>% unique()
row.names(Metabolism) <- Metabolism[,1];Metabolism <- Metabolism[,-1]


level2 <- filter(kegg_TPM_annotion,Pathway1=="Metabolism") %>%
  group_by(Pathway2) %>% 
  summarise(`DCL-1`=sum(`DCL-1`),`DCL-2`=sum(`DCL-2`),
            `XYL-1`=sum(`XYL-1`),`XYL-2`=sum(`XYL-2`),
            `EL-1`=sum(`EL-1`),`EL-2`=sum(`EL-2`),
            `FXL-1`=sum(`FXL-1`),`FXL-2`=sum(`FXL-2`),
            `LGL-1`=sum(`LGL-1`),`LGL-2`=sum(`LGL-2`) ) %>%
  column_to_rownames(var = "Pathway2")

# PcoA ------------------------------------------------------------------
Metabolism <- t(Metabolism) %>% decostand(method = "hellinger")
m = "bray_curtis_Metabolism"
distance <- vegdist(Metabolism, method = 'bray')
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
p_kegg <- ggplot(points, aes(x=PC1, y=PC2, color=group)) + 
  geom_point(alpha=.7, size=2,show.legend=F) +
  labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""),
       title=paste(m," PCoA",sep="")
  )+
  theme_classic()+
  stat_ellipse(level=0.90,show.legend=F,linetype=3,size=0.8) + 
  theme_bw()+theme(panel.grid=element_blank()) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  geom_text_repel(aes(PC1,PC2,label=row.names(points),color=group_2),show.legend=F)
p_kegg

ggsave(file = "KEGG_Metabolism_PcoA_ellipse_0.90_20190506.pdf",p_kegg,width=4, height=4)

# PERMANOVA -------------------------------------------------------------
group <- read.delim('../Datasets/group.txt', sep = '\t', stringsAsFactors = FALSE)
# kegg
adonis_result_dis_kegg <- adonis(distance~group, group, permutations = 9999) 
adonis_result_dis_kegg

  
# heatmap ----------------------------------------------------------------------

percentage <- level2
for (i in 1:ncol(percentage)) {
  percentage[,i] <- percentage[,i]/sum(percentage[,i])
}

col_anno = data.frame(Habitat_types=c(rep(c("Eutrophic", "Mesoeutrophic",
                                           "Oligomesotrophic", "Oligotrophic"),c(4,2,2,2))),
                      row.names=colnames(level2),
                      Group =rep(c("Group I","Group II"),c(4,6)))
ann_colors = list(Habitat_type=c(Eutrophic="#E41A1C",
                                 Mesoeutrophic="#377EB8",
                                 Oligomesotrophic="#4DAF4A",
                                 Oligotrophic="#984EA3"),
                  Group=c(`Group I`="#A65628",`Group II`="#F781BF"))
col_anno

pheatmap(percentage*100,cluster_rows = T,
               cluster_cols = T,
               annotation_col=col_anno, 
               annotation_colors = ann_colors,
               annotation_names_col = T,
               clustering_method = "average",
               display_numbers = T,
               main = "Metabolism KEGG level2 heatmap",
               filename = "kegg_level2_heatmap_20190506.pdf",
                width=9,height=6)


# simper ------------------------------------------------------------------
Metabolism_level2 <- t(level2) %>% decostand(method = "hellinger") %>% as.matrix()
group <- rep(c("Group I","Group II"),c(4,6))
Simper_result <- simper(Metabolism_level2,group,permutations=9999)
summary_result <-  summary(Simper_result)$`Group I_Group II` %>% as.data.frame()
summary_result$adjusted.res2.p = p.adjust(summary_result$p, method = "fdr")





