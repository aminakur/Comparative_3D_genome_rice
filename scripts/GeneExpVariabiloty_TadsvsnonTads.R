setwd("~/Desktop/NYU_Projects/Amina_GeneExp/data/")

library(dplyr)
library(tidyr)
library(stringr)
genexp = read.table("GeneExp/Normal_Salt_mean.expd8sample.txt", header = T, sep="\t")
genexp = genexp[grepl("Os", genexp$TrID),]
genexp = genexp %>% separate(TrID, into = c("geneID", NA), sep="-", remove = F)
genexp$geneID = str_replace_all(genexp$geneID, 't', 'g')
#34609 genes


tadgenes = read.table("TADs/TADs_5kb_intersect80.genes.allwithin", header = T, sep="\t")
tadgenes$tadID = paste("T", c(1:nrow(tadgenes)), sep="")
ntgenes = read.table("TADs/nonTADs_5kb_intersect80.genes.allwithin", header = T, sep="\t")

tadgenes$genes = sapply(tadgenes$genes, function(v){strsplit(v, ',')[[1]]}, USE.NAMES = FALSE)
ntgenes$genes = sapply(ntgenes$genes, function(v){strsplit(v, ',')[[1]]}, USE.NAMES = FALSE)

tadgenes$Identity = rep("TAD", times = nrow(tadgenes))
ntgenes$Identity = rep("nonTAD", times = nrow(ntgenes))

tmp = ntgenes
colnames(tmp) = colnames(tadgenes)
allgenes = rbind(tadgenes, tmp)
allgenes = allgenes[,c(1,7,4,2,3,5,6)]

genexpd = list()
N_genexpd = c()
for (i in 1:nrow(allgenes)){
  genexpd[[i]] = intersect(allgenes$genes[[i]], genexp$geneID)
  N_genexpd = append(N_genexpd, length(intersect(allgenes$genes[[i]], genexp$geneID)))
}
allgenes$genexpd = genexpd 
allgenes$N_genexpd = N_genexpd 

dim(allgenes[which(allgenes$N_genexpd > 4),])   #1247 Tads and nonTads that have at least 5 genes expressed
dim(allgenes[which(allgenes$N_genexpd > 4 & allgenes$Identity == "TAD"),])  #678 Tads
dim(allgenes[which(allgenes$N_genexpd > 4 & allgenes$Identity == "nonTAD"),])     #569 nonTads

gene_new = (allgenes[which(allgenes$N_genexpd > 4),])
summary(gene_new[which(gene_new$Identity == "TAD"),]$N_genexpd)
summary(gene_new[which(gene_new$Identity == "nonTAD"),]$N_genexpd)

tmp = gene_new[which(gene_new$N_genexpd < 33),]
summary(tmp[which(tmp$Identity == "TAD"),]$N_genexpd)
summary(tmp[which(tmp$Identity == "nonTAD"),]$N_genexpd)

#######################################################################################################
#1. Are genes within Tads more correlated (within a timepoint) than those within nonTads

#Trying CV within one timepoint of genes inside Tad and nonTad boundaries
df = data.frame(matrix(ncol = 10, nrow = 0))
colnames(df) = c("CV_bwGene_S1", "CV_bwGene_S4", "CV_bwGene_S6", "CV_bwGene_S7", "CV_bwGene_S8",
                 "CV_bwGene_W1", "CV_bwGene_W4", "CV_bwGene_W6", "CV_bwGene_W7", "CV_bwGene_W8")
for (i in 1:nrow(gene_new)){
  t = as.character(gene_new$genexpd[[i]])
  t2 = genexp[which(genexp$geneID %in% t),]
  for (j in 3:ncol(t2)){
    cv = sd(t2[,j], na.rm = T)/mean(t2[,j], na.rm = T)
    df[i,j-2] = cv
  }
}
tmp = cbind(gene_new, df)
gene_new = tmp
t.test(gene_new$CV_bwGene_W1 ~ gene_new$Identity)
t.test(gene_new$CV_bwGene_W4 ~ gene_new$Identity)
t.test(gene_new$CV_bwGene_W6 ~ gene_new$Identity)
t.test(gene_new$CV_bwGene_W7 ~ gene_new$Identity)
t.test(gene_new$CV_bwGene_W8 ~ gene_new$Identity)
t.test(gene_new$CV_bwGene_S1 ~ gene_new$Identity)
t.test(gene_new$CV_bwGene_S4 ~ gene_new$Identity)
t.test(gene_new$CV_bwGene_S6 ~ gene_new$Identity)
t.test(gene_new$CV_bwGene_S7 ~ gene_new$Identity)
t.test(gene_new$CV_bwGene_S8 ~ gene_new$Identity)
###Significant difference for just normalized mean data
#so, genes within TADs tend to be more correlated in expression amongst them compared to genes lying outside TADs
#####No significant diff for log transformed data

#Plotting
library(ggpubr)
CV_betweenGenes = c(gene_new$CV_bwGene_W1, gene_new$CV_bwGene_W4, gene_new$CV_bwGene_W6, 
                    gene_new$CV_bwGene_W7, gene_new$CV_bwGene_W8)
DevTime = c(rep("0 min", times=nrow(gene_new)), rep("60 min", times=nrow(gene_new)), rep("180 min", times=nrow(gene_new)), 
            rep("240 min", times=nrow(gene_new)), rep("5 days", times=nrow(gene_new)))
Group = c(rep(gene_new$Identity, times=5))
plotDat = data.frame(CV_betweenGenes, DevTime, Group)
CVbetweenGenesW_pertimepoint_plot <- ggboxplot(plotDat, x = "DevTime", y = "CV_betweenGenes",
               color = "Group",  xlab = "Time (Normal)")

CV_betweenGenes = c(gene_new$CV_bwGene_S1, gene_new$CV_bwGene_S4, gene_new$CV_bwGene_S6, 
                    gene_new$CV_bwGene_S7, gene_new$CV_bwGene_S8)
DevTime = c(rep("0 min", times=nrow(gene_new)), rep("60 min", times=nrow(gene_new)), rep("180 min", times=nrow(gene_new)), 
            rep("240 min", times=nrow(gene_new)), rep("5 days", times=nrow(gene_new)))
Group = c(rep(gene_new$Identity, times=5))
plotDat = data.frame(CV_betweenGenes, DevTime, Group)
CVbetweenGenesS_pertimepoint_plot <- ggboxplot(plotDat, x = "DevTime", y = "CV_betweenGenes",
                                              color = "Group",  xlab = "Time (Salt)")

CVbetweenGenes_pertimepoint_plot = ggarrange(CVbetweenGenesW_pertimepoint_plot, CVbetweenGenesS_pertimepoint_plot, 
                                             ncol=1, nrow=2, labels = c('A', 'B'), common.legend = TRUE)
CVbetweenGenes_pertimepoint_plot

#Plotting by taking the mean of betweenCVs across conditions
tmp = gene_new[,c(2,10:19)]
tmp$meanCVbetweenGenes_Normal = rowMeans(tmp[,c(7:11)])
tmp$meanCVbetweenGenes_Salt = rowMeans(tmp[,c(2:6)])
meanCVbetweenGenes = c(tmp$meanCVbetweenGenes_Normal, tmp$meanCVbetweenGenes_Salt)
trt = c(rep("Normal", times=nrow(tmp)), rep("Salt", times=nrow(tmp)))
Group = c(rep(gene_new$Identity, times=2))
plotDat = data.frame(meanCVbetweenGenes, trt, Group)
meanCVbetweenGenes_plot <- ggboxplot(plotDat, x = "trt", y = "meanCVbetweenGenes",
                                               color = "Group", xlab = "")
meanCVbetweenGenes_plot

#Testing whether the higher CV for nonTADs is due to the fact that nonTADs contain high number of genes
#And Plotting the results
tmp = gene_new[which(gene_new$N_genexpd < 33),]     #max number of genes exp in TAD is 32
t.test(tmp$CV_bwGene_S1 ~ tmp$Identity); #p-value = 1.64e-12
t.test(tmp$CV_bwGene_S4 ~ tmp$Identity); #p-value = 6.489e-12
t.test(tmp$CV_bwGene_S6 ~ tmp$Identity); #p-value = 1.53e-10
t.test(tmp$CV_bwGene_S7 ~ tmp$Identity); #p-value = 2.548e-09
t.test(tmp$CV_bwGene_S8 ~ tmp$Identity); #p-value = 7.41e-13
t.test(tmp$CV_bwGene_W1 ~ tmp$Identity); #p-value = 2.881e-12
t.test(tmp$CV_bwGene_W4 ~ tmp$Identity); #p-value = 1.068e-11
t.test(tmp$CV_bwGene_W6 ~ tmp$Identity); #p-value = 6.199e-12 
t.test(tmp$CV_bwGene_W7 ~ tmp$Identity); #p-value = 7.955e-08
t.test(tmp$CV_bwGene_W8 ~ tmp$Identity)  #p-value = 1.683e-08

CV_betweenGenes = c(tmp$CV_bwGene_W1, tmp$CV_bwGene_W4, 
                    tmp$CV_bwGene_W6, tmp$CV_bwGene_W7, tmp$CV_bwGene_W8)
DevTime = c(rep("0 min", times=nrow(tmp)), rep("60 min", times=nrow(tmp)), 
            rep("180 min", times=nrow(tmp)),rep("240 min", times=nrow(tmp)), 
            rep("5 days", times=nrow(tmp)))
Group = c(rep(tmp$Identity, times=5))
plotDat = data.frame(CV_betweenGenes, DevTime, Group)
CVbetweenGenesW_pertimepoint_plot2 <- ggboxplot(plotDat, x = "DevTime", y = "CV_betweenGenes",
                                                color = "Group", xlab = "Time (Normal)")

CV_betweenGenes = c(tmp$CV_bwGene_S1, tmp$CV_bwGene_S4, 
                    tmp$CV_bwGene_S6, tmp$CV_bwGene_S7, tmp$CV_bwGene_S8)
DevTime = c(rep("0 min", times=nrow(tmp)), rep("60 min", times=nrow(tmp)), 
            rep("180 min", times=nrow(tmp)),rep("240 min", times=nrow(tmp)), 
            rep("5 days", times=nrow(tmp)))
Group = c(rep(tmp$Identity, times=5))
plotDat = data.frame(CV_betweenGenes, DevTime, Group)
CVbetweenGenesS_pertimepoint_plot2 <- ggboxplot(plotDat, x = "DevTime", y = "CV_betweenGenes",
                                                color = "Group", xlab = "Time (Salt)")


CVbetweenGenes_pertimepoint_plot2 = ggarrange(CVbetweenGenesW_pertimepoint_plot, CVbetweenGenesS_pertimepoint_plot, 
                                              ncol=1, nrow=2, labels = c('A', 'B'), common.legend = TRUE)
CVbetweenGenes_pertimepoint_plot2

#######################################################################################################
#2. Are genes within Tads more expressed than within  nonTads

#Trying mean within one timepoint of genes inside Tad and nonTad boundaries
df = data.frame(matrix(ncol = 10, nrow = 0))
colnames(df) = c("m_bwGene_S1", "m_bwGene_S4", "m_bwGene_S6", "m_bwGene_S7", "m_bwGene_S8",
                 "m_bwGene_W1", "m_bwGene_W4", "m_bwGene_W6", "m_bwGene_W7", "m_bwGene_W8")
for (i in 1:nrow(gene_new)){
  t = as.character(gene_new$genexpd[[i]])
  t2 = genexp[which(genexp$geneID %in% t),]
  for (j in 3:ncol(t2)){
    m = mean(t2[,j], na.rm = T)
    df[i,j-2] = m
  }
}
tmp = cbind(gene_new, df)
gene_new = tmp
t.test(gene_new$m_bwGene_W1 ~ gene_new$Identity)
t.test(gene_new$m_bwGene_W4 ~ gene_new$Identity)
t.test(gene_new$m_bwGene_W6 ~ gene_new$Identity)
t.test(gene_new$m_bwGene_W7 ~ gene_new$Identity)
t.test(gene_new$m_bwGene_W8 ~ gene_new$Identity)
t.test(gene_new$m_bwGene_S1 ~ gene_new$Identity)
t.test(gene_new$m_bwGene_S4 ~ gene_new$Identity)
t.test(gene_new$m_bwGene_S6 ~ gene_new$Identity)
t.test(gene_new$m_bwGene_S7 ~ gene_new$Identity)
t.test(gene_new$m_bwGene_S8 ~ gene_new$Identity)
###No Significant difference for just normalized mean data
#so, genes within TADs dont tend to be more highly expressed

#Just testing whether this is due to a few extreme mean expression values
tmp = tmp[which(tmp$m_bwGene_S1 < 10000),]
t.test(tmp$m_bwGene_S1 ~ tmp$Identity)
#Clearly doesnt make a diff, both Tads and nonTads have similar expression levels.


#Trying the same thing as above but now we have three nonoverlapping groups -- tad, nontad, and tad boundary
#Trying mean within one timepoint of genes inside Tad and nonTad boundaries
gene3group = read.table("TADs/tad.nontad.br_5kb_intersect80.genes.allwithin", header = T, sep="\t")
gene3group$genes = sapply(gene3group$genes, function(v){strsplit(v, ',')[[1]]}, USE.NAMES = FALSE)

genexpd = list()
N_genexpd = c()
for (i in 1:nrow(gene3group)){
  genexpd[[i]] = intersect(gene3group$genes[[i]], genexp$geneID)
  N_genexpd = append(N_genexpd, length(intersect(gene3group$genes[[i]], genexp$geneID)))
}
gene3group$genexpd = genexpd 
gene3group$N_genexpd = N_genexpd 
gene3group_new = (gene3group[which(gene3group$N_genexpd > 0),])

df = data.frame(matrix(ncol = 10, nrow = 0))
colnames(df) = c("m_bwGene_S1", "m_bwGene_S4", "m_bwGene_S6", "m_bwGene_S7", "m_bwGene_S8",
                 "m_bwGene_W1", "m_bwGene_W4", "m_bwGene_W6", "m_bwGene_W7", "m_bwGene_W8")
for (i in 1:nrow(gene3group_new)){
  t = as.character(gene3group_new$genexpd[[i]])
  t2 = genexp[which(genexp$geneID %in% t),]
  for (j in 3:ncol(t2)){
    m = mean(t2[,j], na.rm = T)
    df[i,j-2] = m
  }
}
tmp = cbind(gene3group_new, df)
tmp %>%  group_by(Identity) %>%
  get_summary_stats(m_bwGene_W1, type = "mean")

tmp %>% pairwise_t_test(m_bwGene_W1 ~ Identity, p.adjust.method = "bonferroni")
1#Still no significant difference

#######################################################################################################
#3. Are genes within Tads more correlated across developmental timepoints than within nonTads

#Trying CV between timepoints for each gene with TADs (or nonTads) and then averaging over the set with boundaries
meanCVwiGenes1 = meanCVwiGenes2 = meanCVwiGenes3 = meanCVwiGenes4 = c()
for (i in 1:nrow(gene_new)){
  t = as.character(gene_new$genexpd[[i]])
  t2 = genexp[which(genexp$geneID %in% t),]
  cv1 = cv2 = cv3 = cv4 = c()
  for (j in 1:nrow(t2)){
    tmp1 = as.numeric(t2[j, c(3:7)])
    tmp2 = sd(tmp1, na.rm = T)/mean(tmp1, na.rm = T)
    cv1 = append(cv1, tmp2)
    tmp1 = as.numeric(t2[j, c(8:12)])
    tmp2 = sd(tmp1, na.rm = T)/mean(tmp1, na.rm = T)
    cv2 = append(cv2, tmp2)
    tmp1 = as.numeric(t2[j, c(3,7)])
    tmp2 = sd(tmp1, na.rm = T)/mean(tmp1, na.rm = T)
    cv3 = append(cv3, tmp2)
    tmp1 = as.numeric(t2[j, c(8,12)])
    tmp2 = sd(tmp1, na.rm = T)/mean(tmp1, na.rm = T)
    cv4 = append(cv4, tmp2)
    
  }
  meanCVwiGenes1 = append(meanCVwiGenes1, mean(cv1, na.rm = T))
  meanCVwiGenes2 = append(meanCVwiGenes2, mean(cv2, na.rm = T))
  meanCVwiGenes3 = append(meanCVwiGenes3, mean(cv3, na.rm = T))
  meanCVwiGenes4 = append(meanCVwiGenes4, mean(cv4, na.rm = T))
}
gene_new$WmeanCV_wiGene = meanCVwiGenes2
gene_new$SmeanCV_wiGene = meanCVwiGenes1
gene_new$WmeanCV1.8_wiGene = meanCVwiGenes4
gene_new$SmeanCV1.8_wiGene = meanCVwiGenes3
t.test(gene_new$WmeanCV_wiGene ~ gene_new$Identity)
t.test(gene_new$SmeanCV_wiGene ~ gene_new$Identity)
t.test(gene_new$WmeanCV1.8_wiGene ~ gene_new$Identity)
t.test(gene_new$SmeanCV1.8_wiGene ~ gene_new$Identity)
#No diff when I tested for all times separately in normal and salt, or even when I only take the dev time as 0min and 5 days

#boxplot(gene_new[which(gene_new$Identity == "TAD"),]$meanCV_wiGene, gene_new[which(gene_new$Identity == "nonTAD"),]$meanCV_wiGene)
###across developmental stage though the TAD and nonTad genes on average display the same CV. Check to see whether this is true for boundary genes also

CV_wiGenes = c(gene_new$WmeanCV_wiGene, gene_new$SmeanCV_wiGene, gene_new$WmeanCV1.8_wiGene, gene_new$SmeanCV1.8_wiGene)
Env = c(rep("W_all", times=nrow(gene_new)), rep("S_all", times=nrow(gene_new)),
        rep("W_0min.5days", times=nrow(gene_new)), rep("S_0min.5days", times=nrow(gene_new)))
Group = c(rep(gene_new$Identity, times=4))
plotDat = data.frame(CV_wiGenes, Env, Group)
CVwiGenes_throughDevTime_plot <- ggboxplot(plotDat, x = "Env", y = "CV_wiGenes",
                                               color = "Group")
CVwiGenes_throughDevTime_plot

############################################################################################################################################
#4. What is the relationship (if any) between genes that are within Tads (and nonTads) of FC and CVofFC between conditions
  
#Trying CV between W and S for the last timepoint (T8 = 5 days) for each gene with TADs (or nonTads)
CV_of_FC = c(); meanFC = c()
absCV_of_FC = c(); absmeanFC = c()
FC_T8_genes = list()
for (i in 1:nrow(gene_new)){
  t = as.character(gene_new$genexpd[[i]])
  t2 = genexp[which(genexp$geneID %in% t),]
  fc = c()
  for (j in 1:nrow(t2)){
    tmp1 = log(abs(t2[j,12]-t2[j,7]), 2)
    tmp2 = ifelse((t2[j,12]-t2[j,7])<0, 0-tmp1, tmp1)
    fc = append(fc, tmp2)
  }
  FC_T8_genes[[i]] = fc
  fc = fc[!is.infinite(fc)]
  tmp1 = sd(fc, na.rm = T)/mean(fc, na.rm = T)
  CV_of_FC = append(CV_of_FC, tmp1)
  tmp1 = sd(abs(fc), na.rm = T)/mean(abs(fc), na.rm = T)
  absCV_of_FC = append(absCV_of_FC, tmp1)
  meanFC = append(meanFC, mean(fc, na.rm = T))
  absmeanFC = append(absmeanFC, mean(abs(fc), na.rm = T))
}

bwCond_comp = data.frame(gene_new$Identity, meanFC, absmeanFC, CV_of_FC, absCV_of_FC)
tmp = cbind(gene_new, bwCond_comp[,-1])
gene_new = tmp
bwCond_comp$FC_T8_genes = FC_T8_genes
bwCond_comp$Nexp = gene_new$N_genexpd
names(bwCond_comp)[1] = "Group"
t.test(bwCond_comp$meanFC ~ bwCond_comp$Group)
#summary(bwCond_comp[which(bwCond_comp$Group == "TAD"),]$absmeanFC);
#summary(bwCond_comp[which(bwCond_comp$Group == "nonTAD"),]$absmeanFC);  

t.test(bwCond_comp$CV_of_FC ~ bwCond_comp$Group)
#summary(bwCond_comp[which(bwCond_comp$Group == "TAD"),]$CV_of_FC);
#summary(bwCond_comp[which(bwCond_comp$Group == "nonTAD"),]$CV_of_FC);  
p1 = ggboxplot(bwCond_comp, x = "Group", y = "CV_of_FC",
          color = "Group")
ggpar(p1, ylim=c(-600, 800))

t.test(bwCond_comp$absmeanFC ~ bwCond_comp$Group)
#t = -2.4629, df = 1205.7, p-value = 0.01392
t.test(bwCond_comp$absmeanFC ~ bwCond_comp$Group, alternative = "less")
#t = -2.4629, df = 1205.7, p-value = 0.00696
summary(bwCond_comp[which(bwCond_comp$Group == "TAD"),]$absmeanFC);
summary(bwCond_comp[which(bwCond_comp$Group == "nonTAD"),]$absmeanFC);  

t.test(bwCond_comp$absCV_of_FC ~ bwCond_comp$Group)
#t = 3.0384, df = 1236, p-value = 0.002428
t.test(bwCond_comp$absCV_of_FC ~ bwCond_comp$Group, alternative = "greater")
#t = 3.0384, df = 1236, p-value = 0.001214
summary(bwCond_comp[which(bwCond_comp$Group == "TAD"),]$absCV_of_FC);
summary(bwCond_comp[which(bwCond_comp$Group == "nonTAD"),]$absCV_of_FC);

names(bwCond_comp)[1] = "Group"
##Plotting 
library(rstatix)
stat.test <- bwCond_comp %>% t_test(absmeanFC ~ Group)
stat.test <- stat.test %>% add_xy_position(x = "Group")
p1 = ggboxplot(bwCond_comp, x = "Group", y = "absmeanFC",
                         color = "Group", xlab="", ylab = "abs(meanFC)") +
      stat_pvalue_manual(stat.test, label = "p", tip.length = 0.01)

stat.test <- bwCond_comp %>% t_test(absCV_of_FC ~ Group)
stat.test <- stat.test %>% add_xy_position(x = "Group")
p2 = ggboxplot(bwCond_comp, x = "Group", y = "absCV_of_FC",
               color = "Group", xlab="", ylab = "abs(coefVarFC)") +
  stat_pvalue_manual(stat.test, label = "p", tip.length = 0.01)

bwCond_FC = ggarrange(p1, p2, ncol=2, nrow=1, labels = c('A', 'B'), common.legend = TRUE)

###Doing the same thing for abs variables above taking into account the number of genes 
#Will take the cutoff of max genes expressed at 32
tmp = bwCond_comp[which(bwCond_comp$Nexp <= 32),]
t.test(tmp$absmeanFC ~ tmp$Group, alterative = "greater")
#t = -2.8274, df = 805.76, p-value = 0.004808

t.test(tmp$absCV_of_FC ~ tmp$Group)
#t = 2.9313, df = 717.2, p-value = 0.003482

##Plotting 
stat.test <- tmp %>% t_test(absmeanFC ~ Group)
stat.test <- stat.test %>% add_xy_position(x = "Group")
p1 = ggboxplot(tmp, x = "Group", y = "absmeanFC",
               color = "Group", xlab="", ylab = "abs(meanFC)") +
  stat_pvalue_manual(stat.test, label = "p", tip.length = 0.01)

stat.test <- tmp %>% t_test(absCV_of_FC ~ Group)
stat.test <- stat.test %>% add_xy_position(x = "Group")
p2 = ggboxplot(tmp, x = "Group", y = "absCV_of_FC",
               color = "Group", xlab="", ylab = "abs(coefVarFC)") +
  stat_pvalue_manual(stat.test, label = "p", tip.length = 0.01)

bwCond_FC2 = ggarrange(p1, p2, ncol=2, nrow=1, labels = c('A', 'B'), common.legend = TRUE)

##################################################################################################################################
#5. Testing for difference in specificity, GeneGC, TrLen between Tad and nonTad genes

tao = read.table("~/Desktop/NYU_Projects/RiceSalinity/Data/QuanGen/Ind_Wet_covariates.all.txt", header = T, sep = "\t")
tao = tao[,c(1,6:8)]
names(tao)[1] = "TrID"
tao$TrID = str_replace_all(tao$TrID, 'OS', 'Os')
tao$TrID = str_replace_all(tao$TrID, 'T', 't')
tao$geneID = str_replace_all(tao$TrID, 't', 'g')
tao = tao %>% separate(geneID, into = c("geneID", "Var"), sep="-", remove = F)

tmp = read.table("~/Desktop/NYU_Projects/Amina_GeneExp/data/AzucenaRef/OsNIP_emart.txt", header = T, sep = "\t")
tmp2 = merge(tao[,c(1,5,4)], tmp[,c(6,1,2)], by="TrID")
tmp2 = tao
#tmp2 = tao[!duplicated(tao$geneID), ]
tmp2 = tao[,c(5,2,3,4)] %>% group_by(geneID) %>% summarise(across(everything(), ~ mean(.x, na.rm = TRUE)))
tao = as.data.frame(tmp2)

df = data.frame(matrix(ncol = 3, nrow = 0))
colnames(df) = c("TrLen_bwGene", "TrGC_bwGene", "tao_bwGene")
for (i in 1:nrow(gene_new)){
  t = as.character(gene_new$genexpd[[i]])
  t2 = tao[which(tao$geneID %in% t),]
  for (j in 2:4){
    m = mean(t2[,j], na.rm = T)
    df[i,j-1] = m
  }
}
tmp = cbind(gene_new, df)
gene_new = tmp
t.test(tmp$tao_bwGene ~ tmp$Identity)
t.test(tmp$TrGC_bwGene ~ tmp$Identity)
#mean of TADs =47.12, nonTADs = 46.62
t.test(tmp$TrLen_bwGene ~ tmp$Identity)
#gene GC content is sig higher (p-val = 0.02823) and Tr len is sig lower (p-val = 0.0004632) lower in TADs
#No diff in specificity b/w Tads and non-Tads

names(tmp)[2] = "Group"
##Plotting 
stat.test <- tmp %>% t_test(TrGC_bwGene ~ Group)
stat.test <- stat.test %>% add_xy_position(x = "Group")
p1 = ggboxplot(tmp, x = "Group", y = "TrGC_bwGene",
               color = "Group", xlab="", ylab = "TrGC") +
  stat_pvalue_manual(stat.test, label = "p", tip.length = 0.01)

stat.test <- tmp %>% t_test(TrLen_bwGene ~ Group)
stat.test <- stat.test %>% add_xy_position(x = "Group")
p2 = ggboxplot(tmp, x = "Group", y = "TrLen_bwGene",
               color = "Group", xlab="", ylab = "TrLen") +
  stat_pvalue_manual(stat.test, label = "p", tip.length = 0.01)

TR_GC_len = ggarrange(p1, p2, ncol=2, nrow=1, labels = c('A', 'B'), common.legend = TRUE)

##################################################################################################################################
#DO NOT RUN THIS
###7. Testing for differences in in DNA methylation between Tad boundaries, domains (minus boundary) and nonTad domains
#NOTE: This block does NOT need to be run (see the end of block), run the next block instead. 

methyl = fread("TADs/methylation_from_zoe.bedgraph", header = F, sep="\t")

#I will take the mean Met score and also teh sum of the scores and divide that by the number of sites for which we have the Met score (in case the number of sites vary)
#Else the score will be inflated due to len (we know nonTADs seems to be much longer than TADs)
meanMetNormLen = c()
meanMet = c()
for (i in 1:nrow(gene_new)){
  t = methyl[which(methyl$V1 == gene_new$Chr[i] & methyl$V2 >= gene_new$tadStart[i] & methyl$V3 <= gene_new$tadStop[i])]
  t$V5 = t$V3-t$V2
  meanMet[i] = mean(t$V4)
  meanMetNormLen[i] = sum(t$V4)/sum(t$V5)
}

gene_new$domainMeanMet_NormLen = meanMetNormLen
gene_new$domainMeanMet = meanMet
t.test(meanMet ~ gene_new$Identity)
mean(methyl$V4)
#Background whole genome mean is 18.4314, nonTAD mean is 17.03169, TAD mean is 14.41779
t.test(meanMetNormLen ~ gene_new$Identity)
mean(sum(methyl$V4)/sum(methyl$V3-methyl$V2))
#Background whole genome mean by normalized length is 12.64794, nonTAD mean is 11.642934, TAD mean is 9.646157
#From both these measures it seems like TADs are undermethylated.

#Estimating what the meanMet is for the boundaries -- would expect them to be lower
#I am defining boundary vs nonboundary genes as boundary +/- 2.5kb and then gene is gene + 1000bp promoter (taking care of directionality). Any overlap would give you genes within boundaries
#I. Using Cris' data
tadbr2.5kb_genes = read.table("TADs/TADboundarys_5kb_intersect80.genes.any", header = T, sep="\t")
tmp = unique(tadbr2.5kb_genes[,c(1:3)])

meanMetNormLen = c()
meanMet = c()
for (i in 1:nrow(tmp)){
  t = methyl[which(methyl$V1 == tmp$Chr[i] & methyl$V2 >= tmp$tadbrStart[i] & methyl$V3 <= tmp$tadbrStop[i])]
  t$V5 = t$V3-t$V2
  meanMet[i] = mean(t$V4)
  meanMetNormLen[i] = sum(t$V4)/sum(t$V5)
}
summary(meanMet)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.3078  5.5885  9.2099 11.0153 15.1280 45.0970 
summary(meanMetNormLen)
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.1416  3.4209  5.6629  7.1330  9.6533 35.0406 


#So it seems like that TAD genes have more CG content, but the TAD domains are undermethylated
#Can this be preferential undermethylation?

#Also maybe check GC content over the TAD and nonTAD domains just to be sure.
#On server used the commands below
#>bedtools nuc -fi /scratch/sg7367/AminaGeneExp/data/Azuref/Azucena.chr.fasta -bed TADs_5kb_intersect80.bed > TAD_nuc_composition.txt
#>bedtools nuc -fi /scratch/sg7367/AminaGeneExp/data/Azuref/Azucena.chr.fasta -bed nonTADs_5kb_intersect80.bed > nonTAD_nuc_composition.txt
#remove #from the beginning of first line

t1 = read.table("TADs/TAD_nuc_composition.txt", header = T, sep="\t")
t2 = read.table("TADs/nonTAD_nuc_composition.txt", header = T, sep="\t")
t.test(t1$X6_pct_gc, t2$X7_pct_gc)
#mean of t1 =0.4312878, t2 = 0.4351971 

#I will take the mean Met score and also teh sum of the scores and divide that by the number of sites for which we have the Met score (in case the number of sites vary)
#Else the score will be inflated due to len (we know nonTADs seems to be much longer than TADs)
gene = read.table("AzucenaRef/OSJAZ_OSNIP.gene.bed", header = F, sep = "\t")
df = data.frame(matrix(ncol = 2, nrow = 0))
colnames(df) = c("meanMet", "meanMetNormLen")
for (i in 1:nrow(gene_new)){
  t = as.character(gene_new$genexpd[[i]])
  t2 = gene[which(gene$V5 %in% t),]
  meanMet = c()
  meanMetNormLen = c()
  for(j in 1:nrow(t2)){
    t3 = methyl[which(methyl$V1 == t2$V1[j] & methyl$V2 >= t2$V2[j] & methyl$V3 <= t2$V3[j])]
    t3$V5 = t3$V3-t3$V2
    meanMet[j] = mean(t3$V4)
    meanMetNormLen[j] = sum(t3$V4)/sum(t3$V5)
  }
  df[i,1] = mean(meanMet)
  df[i,2] = mean(meanMetNormLen)
}

gene_new$geneMeanMet = df[,1]
gene_new$geneMeanMet_NormLen = df[,2]
t.test(meanMet ~ gene_new$Identity)
mean(methyl$V4)

#After talking to Michael about this, I think it might be a good idea to compare nonTADs vs TAD boundaries vs TAD domain without TAD boundaries. 
#Also since we have the genome and methylation data available, maybe we should just do this by domain locations instead of genes within domains. 
#Compare GC pct and methylation directly for the three categories.

##################################################################################################################################
###7. Testing for differences in in DNA methylation between Tad boundaries, domains (minus boundary) and nonTad domains

#On server used the commands below
#>bedtools nuc -fi /scratch/sg7367/AminaGeneExp/data/Azuref/Azucena.chr.fasta -bed tad.nontad.br.bed > tad.nontad.br_nuc_composition.txt
#remove #from the beginning of first line

#First I will test whether there is any diff in methylation between these
library(data.table)
library(rstatix)
gc = fread("TADs/tad.nontad.br_nuc_composition.txt", header = T, sep="\t")
gc = gc[,c(1:3,6,8)]
names(gc) = c("Chr", "Start", "Stop", "Identity", "GC")
gc$Identity = as.factor(gc$Identity)

gc %>%  group_by(Identity) %>%
  get_summary_stats(GC, type = "mean")
#mean of tadbr = 0.445, tad = 0.43, and nonTAD = 0.435
gc %>% pairwise_t_test(GC ~ Identity, p.adjust.method = "bonferroni")
#p-val: nonTAD-TAD = 0.02, nonTAD-tadbr = 1.3e-08, TAD-tadbr = 7.69e-27

#go to "/Users/sonalgupta/Desktop/NYU_Projects/Amina_GeneExp/scripts" and run Rscript tad.nontad.br.methyl.phast.rho.R
dat = fread("TADs/met_phast_rho.tad.nontad.br.txt", sep = "\t", header = T)
dat$Identity = as.factor(dat$Identity)
dat %>%  group_by(Identity) %>%
  get_summary_stats(meanMet, type = "mean")

dat %>%  group_by(Identity) %>%
  get_summary_stats(meanMetNormLen, type = "mean")
dat %>% pairwise_t_test(meanMetNormLen ~ Identity, p.adjust.method = "bonferroni")

dat %>%  group_by(Identity) %>%
  get_summary_stats(meanRho, type = "mean")
dat %>% pairwise_t_test(meanRho ~ Identity, p.adjust.method = "bonferroni")

dat %>%  group_by(Identity) %>%
  get_summary_stats(meanPhast, type = "mean")
dat %>% pairwise_t_test(meanPhast ~ Identity, p.adjust.method = "bonferroni")

##################################################################################################################################
###6. Testing for differences in gene pairs separated by TAD boundary vs not

######i
#I am defining boundary vs nonboundary genes asthose gene pairs identified to be either T2T, T2H, H2T and then if there is a TAD boundary (absolute value -- example: chr01 8.5kb) between them, then they are boundary genes
tadbrgenes = read.table("~/Desktop/NYU_Projects/Amina_GeneExp/data/TADs/TADs_5kb_intersect80.boundarysep.genes", header = T, sep = "\t")
nontadbrgenes = read.table("~/Desktop/NYU_Projects/Amina_GeneExp/data/TADs/5kb_80_nonboundarysep.genes", header = T, sep="\t")

tadbrgenes_exp= tadbrgenes[which(tadbrgenes$Gene1 %in% genexp$geneID & tadbrgenes$Gene2 %in% genexp$geneID),]
nontadbrgenes_exp = nontadbrgenes[which(nontadbrgenes$Gene1 %in% genexp$geneID & nontadbrgenes$Gene2 %in% genexp$geneID),]
#1334 pairs tads, and 15129 pair nontad sep are expressed.
table(tadbrgenes_exp$GenePair)
#H2H T2H T2T = 379 629 333 
table(nontadbrgenes_exp$GenePair)
##H2H T2H T2T = 3680 7488 4074 

tadbrgenes_exp$Group = rep("TADsep", times = nrow(tadbrgenes_exp))
nontadbrgenes_exp$Group = rep("nonTADsep", times = nrow(nontadbrgenes_exp))
genepairs = rbind(tadbrgenes_exp[,-c(1,2)], nontadbrgenes_exp)
tmp = c(1:nrow(genepairs))
genepairs$PairID = paste("Pair", tmp, sep="")

##I will now test whether the CV between boundary genes (TAD vs nonTAD) are any different
#Cv between gene pairs 
df = data.frame(matrix(ncol = 20, nrow = 0))
colnames(df) = c("mean_bwGene_S1", "mean_bwGene_S4", "mean_bwGene_S6", "mean_bwGene_S7", "mean_bwGene_S8", 
                 "mean_bwGene_W1", "mean_bwGene_W4", "mean_bwGene_W6", "mean_bwGene_W7", "mean_bwGene_W8",
                 "CV_bwGene_S1", "CV_bwGene_S4", "CV_bwGene_S6", "CV_bwGene_S7", "CV_bwGene_S8", 
                 "CV_bwGene_W1", "CV_bwGene_W4", "CV_bwGene_W6", "CV_bwGene_W7", "CV_bwGene_W8")
for (i in 1:nrow(genepairs)){
  t = c(genepairs[i,1], genepairs[i,2])
  t2 = genexp[which(genexp$geneID %in% t),]
  for (j in 3:ncol(t2)){
    m = mean(t2[,j], na.rm = T)
    df[i,j-2] = m
  }
  for (j in 3:ncol(t2)){
    cv = sd(t2[,j], na.rm = T)/mean(t2[,j], na.rm = T)
    df[i,j+10-2] = cv
  }
}

tmp = cbind(genepairs[,c(13,12,11)], df)
t.test(tmp$mean_bwGene_S1~tmp$Group); t.test(tmp$mean_bwGene_S4~tmp$Group); t.test(tmp$mean_bwGene_S6~tmp$Group); t.test(tmp$mean_bwGene_S7~tmp$Group); t.test(tmp$mean_bwGene_S8~tmp$Group)
t.test(tmp$CV_bwGene_S1~tmp$Group); t.test(tmp$CV_bwGene_S4~tmp$Group); t.test(tmp$CV_bwGene_S6~tmp$Group); t.test(tmp$CV_bwGene_S7~tmp$Group); t.test(tmp$CV_bwGene_S8~tmp$Group)

t.test(tmp$mean_bwGene_W1~tmp$Group); t.test(tmp$mean_bwGene_W4~tmp$Group); t.test(tmp$mean_bwGene_W6~tmp$Group); t.test(tmp$mean_bwGene_W7~tmp$Group); t.test(tmp$mean_bwGene_W8~tmp$Group)
t.test(tmp$CV_bwGene_W1~tmp$Group); t.test(tmp$CV_bwGene_W4~tmp$Group); t.test(tmp$CV_bwGene_W6~tmp$Group); t.test(tmp$CV_bwGene_W7~tmp$Group); t.test(tmp$CV_bwGene_W8~tmp$Group)
#None are significant without taking the type of pairs

boundpair_meanCV = tmp
boundaryGenesmeanbetweenGenesW1_plot <- ggboxplot(boundpair_meanCV, x = "GenePair", y = "mean_bwGene_W1",
                                                color = "Group")
boundaryGenesmeanbetweenGenesW4_plot <- ggboxplot(boundpair_meanCV, x = "GenePair", y = "mean_bwGene_W4",
                                                color = "Group")
boundaryGenesmeanbetweenGenesW6_plot <- ggboxplot(boundpair_meanCV, x = "GenePair", y = "mean_bwGene_W6",
                                                color = "Group")
boundaryGenesmeanbetweenGenesW7_plot <- ggboxplot(boundpair_meanCV, x = "GenePair", y = "mean_bwGene_W7",
                                                color = "Group")
boundaryGenesmeanbetweenGenesW8_plot <- ggboxplot(boundpair_meanCV, x = "GenePair", y = "mean_bwGene_W8",
                                                color = "Group")

boundaryGenesCVbetweenGenesW1_plot <- ggboxplot(boundpair_CVbetween, x = "GenePair", y = "CV_bwGene_W1",
                                               color = "Group")
boundaryGenesCVbetweenGenesW4_plot <- ggboxplot(boundpair_CVbetween, x = "GenePair", y = "CV_bwGene_W4",
                                                color = "Group")
boundaryGenesCVbetweenGenesW6_plot <- ggboxplot(boundpair_CVbetween, x = "GenePair", y = "CV_bwGene_W6",
                                                color = "Group")
boundaryGenesCVbetweenGenesW7_plot <- ggboxplot(boundpair_CVbetween, x = "GenePair", y = "CV_bwGene_W7",
                                                color = "Group")
boundaryGenesCVbetweenGenesW8_plot <- ggboxplot(boundpair_CVbetween, x = "GenePair", y = "CV_bwGene_W8",
                                                color = "Group")
boundaryGenesCVbetweenGenesW_plot = ggarrange(boundaryGenesmeanbetweenGenesW1_plot, boundaryGenesmeanbetweenGenesW4_plot, boundaryGenesmeanbetweenGenesW6_plot,
                                              boundaryGenesmeanbetweenGenesW7_plot, boundaryGenesmeanbetweenGenesW8_plot,
                                              boundaryGenesCVbetweenGenesW1_plot, boundaryGenesCVbetweenGenesW4_plot, boundaryGenesCVbetweenGenesW6_plot,
                                              boundaryGenesCVbetweenGenesW7_plot, boundaryGenesCVbetweenGenesW8_plot,
                                              ncol=5, nrow= 2, common.legend = TRUE, labels = c('A', 'B', 'C', 'D', 'E', 'F','G', 'H', 'I', 'J'))


##Looking at these plots it doesnt look like that within divergent, convergent and tandem also there is any diff.
#testing now
divergent = boundpair_CVbetween[which(boundpair_CVbetween$GenePair == "H2H"),]
t.test(divergent$CV_bwGene_W4~divergent$Group); 
convergent = boundpair_CVbetween[which(boundpair_CVbetween$GenePair == "T2T"),]
t.test(convergent$CV_bwGene_W4~convergent$Group); 
tandem = boundpair_CVbetween[which(boundpair_CVbetween$GenePair == "T2H"),]
t.test(tandem$CV_bwGene_W4~tandem$Group); 
#No diff


########################
######ii
#I am defining boundary vs nonboundary genes as boundary +/- 2.5kb and then gene is gene + 1000bp promoter (taking care of directionality). Any overlap would give you genes within boundaries

#I. Using Cris' data
tadbr2.5kb_genes = read.table("TADs/TADboundarys_5kb_intersect80.genes.any", header = T, sep="\t")
nottadbr2.5kb_genes = read.table("TADs/notTADboundarys_5kb_intersect80.genes.any", header = T, sep="\t")

exp_tadbr = genexp[which(genexp$geneID %in% tadbr2.5kb_genes$geneID),]
exp_nottadbr = genexp[which(genexp$geneID %in% nottadbr2.5kb_genes$x),]
t1 = log10(exp_tadbr$W1)
t2 = log10(exp_nottadbr$W1)

#salt
t.test(exp_tadbr[,3], exp_nottadbr[,3])
t.test(exp_tadbr[,4], exp_nottadbr[,4])
t.test(exp_tadbr[,5], exp_nottadbr[,5])
t.test(exp_tadbr[,6], exp_nottadbr[,6])
t.test(exp_tadbr[,7], exp_nottadbr[,7])
#normal
t.test(exp_tadbr[,8], exp_nottadbr[,8])
t.test(exp_tadbr[,9], exp_nottadbr[,9])
t.test(exp_tadbr[,10], exp_nottadbr[,10])
t.test(exp_tadbr[,11], exp_nottadbr[,11])
t.test(exp_tadbr[,11], exp_nottadbr[,11])

#No sig difference -- could this be due to the sample size difference (~3000 vs ~31500)? 
#Do bootstrapping -- randomly select ~3000 observations from nottadbr estimate mean and compare with tadbr for 1000 rounds??

#Permutation tests
m = 10000
n = 3000
set.seed(1234)
nulldist = c()
for (i in 1:m){
  nulldist[i] = mean(sample(exp_nottadbr[,9], size = n, replace=TRUE))
}
OGteststat =  mean(exp_tadbr[,9])
hist(nulldist)
abline(v=OGteststat, col="red")
lowtail = sum(nulldist <= OGteststat) + 1
uptail = sum(nulldist >= OGteststat) + 1
numerator = min(lowtail, uptail)
denom = m+1
pval = 2*numerator/denom


t_mean = colMeans(exp_tadbr[sapply(exp_tadbr, is.numeric)])
t_cv = sapply(exp_tadbr[,-(1:2)], sd)/t_mean

nt_mean = colMeans(exp_nottadbr[sapply(exp_nottadbr, is.numeric)])
nt_cv = sapply(exp_nottadbr[,-(1:2)], sd)/nt_mean

Mean = c(t_mean, nt_mean)
CV = c(t_cv, nt_cv)
Type = c(rep("Boundary", times=length(t_mean)),rep("nonBoundary", times=length(nt_mean)))
Cond = rep(colnames(genexp[,-(1:2)]), times = 2)
plotmeanexp = data.frame(Type, Cond, Mean, CV)



########################
######iii
#I am defining boundary vs nonboundary genes as boundary +/- 2.5kb and then gene is just the 1000bp promoter (without the gene body; taking care of directionality). Any overlap would give you genes within boundaries

#I. Using Cris' data
tad.br.nontad.promoterwithin = read.table("TADs/tad.nontad.br_5kb_intersect80.promoters.allwithin.txt", header = T, sep="\t")
tad.br.nontad.promoterany = read.table("TADs/tad.nontad.br_5kb_intersect80.promoters.anyoverlap.txt", header = T, sep="\t")

exp_tadbr = genexp[which(genexp$geneID %in% tad.br.nontad.promoterwithin[which(tad.br.nontad.promoterwithin$Identity == "tadbr"),]$geneID),]
exp_nottadbr = genexp[which(genexp$geneID %in% tad.br.nontad.promoterwithin[which(tad.br.nontad.promoterwithin$Identity != "tadbr"),]$geneID),]
exp_tadbody = genexp[which(genexp$geneID %in% tad.br.nontad.promoterwithin[which(tad.br.nontad.promoterwithin$Identity == "TAD"),]$geneID),]
exp_nottadbody = genexp[which(genexp$geneID %in% tad.br.nontad.promoterwithin[which(tad.br.nontad.promoterwithin$Identity == "nonTAD"),]$geneID),]

#salt
t.test(exp_tadbr[,3], exp_nottadbr[,3])
t.test(exp_tadbr[,4], exp_nottadbr[,4])
t.test(exp_tadbr[,5], exp_nottadbr[,5])
t.test(exp_tadbr[,6], exp_nottadbr[,6])
t.test(exp_tadbr[,7], exp_nottadbr[,7])
#normal
t.test(exp_tadbr[,8], exp_nottadbr[,8])
t.test(exp_tadbr[,9], exp_nottadbr[,9])
t.test(exp_tadbr[,10], exp_nottadbr[,10])
t.test(exp_tadbr[,11], exp_nottadbr[,11])
t.test(exp_tadbr[,11], exp_nottadbr[,11])

#No sig difference -- could this be due to the sample size difference (~3000 vs ~31500)? 
########################################################################################################

##8. Does gene exp vary within TAD boundaries with higher/lower insulation score
#I am defining boundary vs nonboundary genes as boundary +/- 2.5kb and then gene is gene + 1000bp promoter (taking care of directionality). Any overlap would give you genes within boundaries
#Of the 2414 TAD boundaries, 1818 boundaries are unique, and further 1607 boundaries have genes within them

tadbrscore_tmp = read.table("TADs/az_5kb_boundaries_insulationScore.bed", header = F, sep="\t")
names(tadbrscore_tmp) = c("Chr", "Start", "Stop", "ID", "InsuSc", "Rand")
summary(tadbrscore_tmp$InsuSc)
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-1.60145 -0.40831 -0.25114 -0.21327 -0.05085  2.06774 

tmp = tadbrscore_tmp[which(tadbrscore_tmp$InsuSc <= -0.40831 | tadbrscore_tmp$InsuSc >= -0.05085),]
tmp$InScoreType = ifelse(tmp$InsuSc <= -0.40831, "Low", "High")

tadbr2.5kb_genes = read.table("TADs/TADboundarys_5kb_intersect80.genes.any", header = T, sep="\t")
names(tadbr2.5kb_genes)[2] = "Start"; names(tadbr2.5kb_genes)[3] = "Stop"
tmp2 = merge(tadbr2.5kb_genes, tmp[,c(1:3,5,7)], by=c("Chr", "Start", "Stop"))
tadbr2.5kb_genes_insulationType = tmp2
#taking genebody and promoter with any overlap, we have 1358 genes

tadbrgenes = tad.br.nontad.promoterwithin[which(tad.br.nontad.promoterwithin$Identity == "tadbr"),]
tmp2 = merge(tadbrgenes, tmp[,c(1:3,5,7)], by=c("Chr", "Start", "Stop"))
tadbrgeneswi_insulationType = tmp2
#taking only promoters within the boundary we have a total of 699 genes

tadbrgenes = tad.br.nontad.promoterany[which(tad.br.nontad.promoterany$Identity == "tadbr"),]
tmp2 = merge(tadbrgenes, tmp[,c(1:3,5,7)], by=c("Chr", "Start", "Stop"))
tadbrgenesany_insulationType = tmp2
#taking only promoters with any overlap the boundary we have a total of 987 genes

exp = genexp[which(genexp$geneID %in% tadbr2.5kb_genes_insulationType$geneID),]
test = merge(tadbr2.5kb_genes_insulationType[,c(1:3,8:10)], exp[,-1], by="geneID")
test = test[,c(2:6,1,7:16)]
tadbr2.5kb_genesExp = test
t.test(tadbr2.5kb_genesExp$W4 ~ tadbr2.5kb_genesExp$InScoreType)

exp = genexp[which(genexp$geneID %in% tadbrgeneswi_insulationType$geneID),]
test = merge(tadbrgeneswi_insulationType[,c(1:3,11:13)], exp[,-1], by="geneID")
test = test[,c(2:6,1,7:16)]
tadbrgeneswi_Exp = test
t.test(tadbrgeneswi_Exp$W4 ~ tadbrgeneswi_Exp$InScoreType)

exp = genexp[which(genexp$geneID %in% tadbrgenesany_insulationType$geneID),]
test = merge(tadbrgenesany_insulationType[,c(1:3,11:13)], exp[,-1], by="geneID")
test = test[,c(2:6,1,7:16)]
tadbrgenesany_Exp = test
t.test(tadbrgenesany_Exp$W4 ~ tadbrgenesany_Exp$InScoreType)
#Taking the genebody and complete promoters within the high vs low insulation types boundaries have no diff in their mean expression
#But taking any promoter overlap, there seems to be almost singnificant difference bwteeen high and low, with high group have low overall expression.

#I am now taking meanExp of the genes within the same boundary
tmp = test %>%
  group_by(Chr, Start, Stop, InsuSc, InScoreType) %>% 
  summarise_if(is.numeric, mean) 
#summarize(genes = paste(geneID, collapse = ',')) %>% 
tadbrgenesany_meanExp = tmp
t.test(tadbrgenesany_meanExp$W1~tadbrgenesany_meanExp$InScoreType)
#The difference reduces when comparing the means of meanExp


