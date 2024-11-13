setwd("~/Desktop/NYU_Projects/Amina_GeneExp/data/")

library(rstatix)
library(ggpubr)
library(ggh4x)
library(data.table)
library(stringr)

#Fig1
t = read.table("TADs/tad.nontad.br_5kb_intersect80.genes.allwithin", header = T, sep="\t")
tmp = t
tmp$Identity = str_replace(tmp$Identity, "tadbr\\b", "TADbr")
tmp$Identity = str_replace(tmp$Identity, "nonTAD\\b", "nonTADbody")
tmp$Identity = str_replace(tmp$Identity, "TAD\\b", "TADbody")
tmp$GenDen = (tmp$n/tmp$Width)*1000
stat.test <- tmp %>% t_test(GenDen ~ Identity, p.adjust.method = "bonferroni")
#stat.test <- stat.test %>% add_xy_position(x = "Identity", dodge = .5)
p1 = ggboxplot(tmp, x = "Identity", y = "GenDen",
               fill = "Identity",palette = c("#00AFBB", "#E7B800", "#FC4E07"),
               outlier.shape = NA,
               bxp.errorbar = T, bxp.errorbar.width = 0.3,
               xlab="", ylab = "GeneDensity") + 
  scale_y_continuous(breaks=seq(0,0.4,0.2), limits=c(0,0.41)) + 
  theme(legend.position="none", 
        axis.text = element_text(size=14),
        axis.title = element_text(size=18, color="black"),
        axis.line = element_line(color = 'black', linewidth = 0.65),
        axis.ticks.length=unit(0.3, "cm")) + 
  stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.01, y.position = c(0.32, 0.36, 0.40)) +  
  guides(x = "axis_truncated", y = "axis_truncated") 

t = read.table("TADs/tad.nontad.br_5kb_intersect80.genicfeatures.txt", header = T, sep="\t")
tmp = t; tmp$Identity = str_replace(tmp$Identity, "tadbr\\b", "TADbr")
tmp$Identity = str_replace(tmp$Identity, "TAD\\b", "TADbody")
tmp$Identity = str_replace(tmp$Identity, "nonTAD\\b", "nonTADbody")
tmp[,-1] %>% group_by(Identity) %>% summarize(meanTrLen = mean(TrLen), meanGeneGC = mean(geneGC))

stat.test <- tmp %>% t_test(TrLen ~ Identity, p.adjust.method = "bonferroni")
tmp$TrLen = tmp$TrLen/1000
p2 = ggboxplot(tmp, x = "Identity", y = "TrLen",
               fill = "Identity",palette = c("#00AFBB", "#E7B800", "#FC4E07"),
               outlier.shape = NA,
               bxp.errorbar = T, bxp.errorbar.width = 0.3,
               xlab="", ylab = "Transcript Length (Kb)") + 
  scale_y_continuous(breaks=seq(0,4.5,1.5), limits=c(0,4.6)) + 
  theme(legend.position="none", 
        axis.text = element_text(size=14),
        axis.title = element_text(size=18, color="black"),
        axis.line = element_line(color = 'black', linewidth = 0.65),
        axis.ticks.length=unit(0.3, "cm")) + 
  stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.01, y.position = c(3.9, 4.2, 4.6)) +  
  guides(x = "axis_truncated", y = "axis_truncated") 

stat.test <- tmp %>% t_test(geneGC ~ Identity, p.adjust.method = "bonferroni")
tmp$geneGC = tmp$geneGC*100
p3 = ggboxplot(tmp, x = "Identity", y = "geneGC",
               fill = "Identity",palette = c("#00AFBB", "#E7B800", "#FC4E07"),
               outlier.shape = NA,
               bxp.errorbar = T, bxp.errorbar.width = 0.3,
               xlab="", ylab = "Gene GC Content (%)") + 
  scale_y_continuous(breaks=seq(10,90,40), limits=c(10,100)) + 
  theme(legend.position="none", 
        axis.text = element_text(size=14),
        axis.title = element_text(size=18, color="black"),
        axis.line = element_line(color = 'black', linewidth = 0.65),
        axis.ticks.length=unit(0.3, "cm")) + 
  stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.01, y.position = c(88, 93, 98)) +  
  guides(x = "axis_truncated", y = "axis_truncated") 


t = fread("TADs/met_phast_rho.tad.nontad.br.txt", sep = "\t", header = T)
tmp = t; tmp$Identity = str_replace(tmp$Identity, "tadbr\\b", "TADbr")
tmp$Identity = str_replace(tmp$Identity, "TAD\\b", "TADbody")
tmp$Identity = str_replace(tmp$Identity, "nonTAD\\b", "nonTADbody")
tmp$Identity = as.factor(tmp$Identity)
#tmp %>%  group_by(Identity) %>%  get_summary_stats(meanMetNormLen, type = "mean")
#ggplot(tmp, aes(x = meanMetNormLen, colour = Identity)) + geom_density()

stat.test <- tmp %>% t_test(meanMetNormLen ~ Identity, p.adjust.method = "bonferroni")
p4 = ggboxplot(tmp, x = "Identity", y = "meanMetNormLen",
               fill = "Identity",palette = c("#00AFBB", "#E7B800", "#FC4E07"),
               outlier.shape = NA,
               bxp.errorbar = T, bxp.errorbar.width = 0.3,
               xlab="", ylab = "Methylation Score") + 
  scale_y_continuous(breaks=seq(0,30,10), limits=c(0,32)) + 
  theme(legend.position="none", 
        axis.text = element_text(size=14),
        axis.title = element_text(size=18, color="black"),
        axis.line = element_line(color = 'black', linewidth = 0.65),
        axis.ticks.length=unit(0.3, "cm")) + 
  stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.01, y.position = c(24, 26.5, 29)) +  
  guides(x = "axis_truncated", y = "axis_truncated") 

ggarrange(p1, p2, p3, p4, ncol=2, nrow=2, labels = c('A', 'B', 'C', 'D'), common.legend = FALSE)



#FigS2
t = read.table("TADs/tad.nontad.br_nuc_composition.txt", header = T, sep="\t")
tmp = t; tmp = tmp[,c(6,8)]
names(tmp) = c("Identity","GCcontent")
tmp$Identity = str_replace(tmp$Identity, "tadbr\\b", "TADbr")
tmp$Identity = str_replace(tmp$Identity, "TAD\\b", "TADbody")
tmp$Identity = str_replace(tmp$Identity, "nonTAD\\b", "nonTADbody")
tmp$Identity = as.factor(tmp$Identity)
tmp %>% group_by(Identity) %>% summarize(meanGC = mean(GCcontent))
stat.test <- tmp %>% t_test(GCcontent ~ Identity, p.adjust.method = "bonferroni")
tmp$GCcontent = tmp$GCcontent*100
s2 = ggboxplot(tmp, x = "Identity", y = "GCcontent",
               fill = "Identity",palette = c("#00AFBB", "#E7B800", "#FC4E07"),
               outlier.shape = NA,
               bxp.errorbar = T, bxp.errorbar.width = 0.3,
               xlab="", ylab = "GC Content (%)") + 
  scale_y_continuous(breaks=seq(30,60,10), limits=c(30,63)) + 
  theme(legend.position="none", 
        axis.text = element_text(size=14),
        axis.title = element_text(size=18, color="black"),
        axis.line = element_line(color = 'black', linewidth = 0.65),
        axis.ticks.length=unit(0.3, "cm")) + 
  stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.01, y.position = c(56, 59, 62)) +  
  guides(x = "axis_truncated", y = "axis_truncated") 

s2


#Fig 2
t1 = read.table("TADs/CV_greenhouse.txt", header = T, sep="\t")
t2 = read.table("TADs/CV_field.txt", header = T, sep="\t")
tmp1 = t1[,c(2,8,13)]
tmp1 %>% group_by(Identity) %>% summarize(meanCV_W1 = mean(CV_bwGene_W1), meanCV_S1 = mean(CV_bwGene_S1))


CV_betweenGenes = c(tmp1$CV_bwGene_W1, tmp1$CV_bwGene_S1)
DevTime = c(rep("Normal", times=nrow(tmp1)), rep("Saline", times=nrow(tmp1)))
Identity = c(rep(tmp1$Identity, times=2))
plotDat = data.frame(CV_betweenGenes, DevTime, Identity)
plotDat$Identity = as.factor(plotDat$Identity)
plotDat$DevTime = as.factor(plotDat$DevTime)
stat.test <- plotDat %>% group_by(DevTime) %>% t_test(CV_betweenGenes ~ Identity) %>% add_significance("p") %>% add_xy_position(x = "DevTime", dodge = 1)
p2 = ggboxplot(plotDat, x = "DevTime", y = "CV_betweenGenes",
               fill = "Identity",palette = c("#00AFBB", "#E7B800"),
               outlier.shape = NA,
               bxp.errorbar = T, bxp.errorbar.width = 0.3,
               xlab="", ylab = "Coefficient of variation") + 
  scale_y_continuous(breaks=seq(0,5,2.5), limits=c(0,5)) + 
  theme(axis.text = element_text(size=16),
        axis.title = element_text(size=20, color="black"),
        axis.line = element_line(color = 'black', linewidth = 0.65),
        axis.ticks.length=unit(0.3, "cm")) + 
  stat_pvalue_manual(stat.test, label = "p.signif", size = 10,  tip.length = 0.01, y.position = 4.75) +  
  guides(x = "axis_truncated", y = "axis_truncated") 



tmp2 = t1[,c(2,9:12,14:17)]
CV_betweenGenes = c(tmp2$CV_bwGene_W4, tmp2$CV_bwGene_W6, 
                    tmp2$CV_bwGene_W7, tmp2$CV_bwGene_W8)
DevTime = c(rep("60 min", times=nrow(tmp2)), rep("180 min", times=nrow(tmp2)), 
            rep("240 min", times=nrow(tmp2)), rep("5 days", times=nrow(tmp2)))
Identity = c(rep(tmp2$Identity, times=4))
plotDat = data.frame(CV_betweenGenes, DevTime, Identity)
plotDat$Identity = as.factor(plotDat$Identity)
#plotDat %>% group_by(Identity) %>% summarize(meanCV = mean(CV_bwGene_W1))
stat.test <- plotDat %>% group_by(DevTime) %>% t_test(CV_betweenGenes ~ Identity) %>% add_significance("p") %>% add_xy_position(x = "DevTime", dodge = 1)
s3.1 = ggboxplot(plotDat, x = "DevTime", y = "CV_betweenGenes",
               fill = "Identity",palette = c("#00AFBB", "#E7B800"),
               outlier.shape = NA,
               bxp.errorbar = T, bxp.errorbar.width = 0.3,
               xlab="", ylab = "Coefficient of variation") + 
  scale_y_continuous(breaks=seq(0,5,2.5), limits=c(0,5)) + ggtitle("Normal Greenhouse") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.text = element_text(size=16),
        axis.title.x = element_text(size=20, color="black"),
        axis.title.y = element_text(size=16, color="black"),
        axis.line = element_line(color = 'black', linewidth = 0.65),
        axis.ticks.length=unit(0.3, "cm")) + 
  stat_pvalue_manual(stat.test, label = "p.signif", size = 10,  tip.length = 0.01, y.position = 4.75) +  
  guides(x = "axis_truncated", y = "axis_truncated") 

CV_betweenGenes = c(tmp2$CV_bwGene_S4, tmp2$CV_bwGene_S6, 
                    tmp2$CV_bwGene_S7, tmp2$CV_bwGene_S8)
DevTime = c(rep("60 min", times=nrow(tmp2)), rep("180 min", times=nrow(tmp2)), 
            rep("240 min", times=nrow(tmp2)), rep("5 days", times=nrow(tmp2)))
Identity = c(rep(tmp2$Identity, times=4))
plotDat = data.frame(CV_betweenGenes, DevTime, Identity)
plotDat$Identity = as.factor(plotDat$Identity)
#plotDat %>% group_by(Identity) %>% summarize(meanCV = mean(CV_bwGene_W1))
stat.test <- plotDat %>% group_by(DevTime) %>% t_test(CV_betweenGenes ~ Identity) %>% add_significance("p") %>% add_xy_position(x = "DevTime", dodge = 1)
s3.2 = ggboxplot(plotDat, x = "DevTime", y = "CV_betweenGenes",
                 fill = "Identity",palette = c("#00AFBB", "#E7B800"),
                 outlier.shape = NA,
                 bxp.errorbar = T, bxp.errorbar.width = 0.3,
                 xlab="", ylab = "Coefficient of variation") + 
  scale_y_continuous(breaks=seq(0,5,2.5), limits=c(0,5)) + ggtitle("Saline Greenhouse") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.text = element_text(size=16),
        axis.title.x = element_text(size=20, color="black"),
        axis.title.y = element_text(size=16, color="black"),
        axis.line = element_line(color = 'black', linewidth = 0.65),
        axis.ticks.length=unit(0.3, "cm")) + 
  stat_pvalue_manual(stat.test, label = "p.signif", size = 10,  tip.length = 0.01, y.position = 4.75) +  
  guides(x = "axis_truncated", y = "axis_truncated") 


str(t2)
library(tidyr)
tmp1 = t2 %>% drop_na(CV_bwGene_Wd)
CV_betweenGenes = c(tmp1$CV_bwGene_Wd, tmp1$CV_bwGene_D)
DevTime = c(rep("Normal", times=nrow(tmp1)), rep("Drought", times=nrow(tmp1)))
Identity = c(rep(tmp1$Identity, times=2))
plotDat = data.frame(CV_betweenGenes, DevTime, Identity)
plotDat$Identity = as.factor(plotDat$Identity)
plotDat$DevTime <- factor(plotDat$DevTime , levels=c("Normal", "Drought"))
stat.test <- plotDat %>% group_by(DevTime) %>% t_test(CV_betweenGenes ~ Identity) %>% add_significance("p") %>% add_xy_position(x = "DevTime", dodge = 1)
s3.3 = ggboxplot(plotDat, x = "DevTime", y = "CV_betweenGenes",
               fill = "Identity",palette = c("#00AFBB", "#E7B800"),
               outlier.shape = NA,
               bxp.errorbar = T, bxp.errorbar.width = 0.3,
               xlab="", ylab = "Coefficient of variation") + 
  scale_y_continuous(breaks=seq(0,5,2.5), limits=c(0,5)) + ggtitle("Drought Field") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.text = element_text(size=16),
        axis.title.x = element_text(size=20, color="black"),
        axis.title.y = element_text(size=16, color="black"),
        axis.line = element_line(color = 'black', linewidth = 0.65),
        axis.ticks.length=unit(0.3, "cm")) +
  stat_pvalue_manual(stat.test, label = "p.signif", size = 10,  tip.length = 0.01, y.position = 4.75) +  
  guides(x = "axis_truncated", y = "axis_truncated") 

tmp1 = t2 %>% drop_na(CV_bwGene_Ws)
CV_betweenGenes = c(tmp1$CV_bwGene_Wd, tmp1$CV_bwGene_D)
DevTime = c(rep("Normal", times=nrow(tmp1)), rep("Saline", times=nrow(tmp1)))
Identity = c(rep(tmp1$Identity, times=2))
plotDat = data.frame(CV_betweenGenes, DevTime, Identity)
plotDat$Identity = as.factor(plotDat$Identity)
plotDat$DevTime <- factor(plotDat$DevTime , levels=c("Normal", "Saline"))
stat.test <- plotDat %>% group_by(DevTime) %>% t_test(CV_betweenGenes ~ Identity) %>% add_significance("p") %>% add_xy_position(x = "DevTime", dodge = 1)
s3.4 = ggboxplot(plotDat, x = "DevTime", y = "CV_betweenGenes",
                 fill = "Identity",palette = c("#00AFBB", "#E7B800"),
                 outlier.shape = NA,
                 bxp.errorbar = T, bxp.errorbar.width = 0.3,
                 xlab="", ylab = "Coefficient of variation") + 
  scale_y_continuous(breaks=seq(0,5,2.5), limits=c(0,5)) + ggtitle("Saline Field") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.text = element_text(size=16),
        axis.title.x = element_text(size=20, color="black"),
        axis.title.y = element_text(size=16, color="black"),
        axis.line = element_line(color = 'black', linewidth = 0.65),
        axis.ticks.length=unit(0.3, "cm")) +
  stat_pvalue_manual(stat.test, label = "p.signif", size = 10,  tip.length = 0.01, y.position = 4.75) +  
  guides(x = "axis_truncated", y = "axis_truncated") 


s3 = ggarrange(s3.1, s3.3, s3.2, s3.4, ncol=2, nrow=2, labels = c('A','C','B','D'), widths = c(1.7, 1), common.legend = TRUE)
