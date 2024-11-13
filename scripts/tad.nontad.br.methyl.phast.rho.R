#Will compare GC and methylation 

setwd("~/Desktop/NYU_Projects/Amina_GeneExp/data/")
library(GenomicRanges)

#1. intersect the genome co-ordinates with TADs (+/- 2.5kb) to get the non-TAD domain (minus the TAD boundary) bed
gen = read.table("AzucenaRef/AzuChr.bed", header = F, sep="\t")
tad = read.table("TADs/TADs_5kb_intersect80.bed", header = F, sep = "\t")
gengr = GRanges(seqnames=gen$V1,
                IRanges(start=gen$V2,end=gen$V3))
tadgr = GRanges(seqnames=tad$V1,
                IRanges(start=tad$V2-2500,end=tad$V3+2500), tadid = tad$V4)
nontadgr = GenomicRanges::setdiff(gengr, tadgr)
nontad = as(nontadgr, "data.frame")
t = c(1:nrow(nontad))
tmp = paste("NT", t, sep = "")
nontad$ID = tmp
remove(t); remove(tmp)
nontad = nontad[,-5]
names(nontad) = c("Chr","Start","Stop","Width","ID")
nontad$Identity = "nonTAD"

#2. #Get the Tad boundary bed
tadbr = tad[,c(1,2,4)]
tmp = tad[,c(1,3,4)]
colnames(tadbr) = c("Chr", "Pos", "ID")
colnames(tmp) = c("Chr", "Pos", "ID")
t = rbind(tadbr, tmp)
tadbr = t
#So there are a total of 2414 TAD boundaries (for 1207 tads), but some of these boundaries are replicated, so I remove these
#Removing Duplicates (1818 boundaries remain)
t = unique(tadbr[,c(1,2)])
#boundary is defined as the tad domain boundary +/- 2500, 
tadbr2.5kb = t
tadbr2.5kb$Start = t$Pos - 2500; tadbr2.5kb$Stop = t$Pos + 2500;
tadbr2.5kb$Width = tadbr2.5kb$Stop - tadbr2.5kb$Start
t = c(1:nrow(tadbr2.5kb))
tmp = paste("TadBr", t, sep = "")
tadbr2.5kb$ID = tmp
tadbr2.5kb = tadbr2.5kb[,-2]
tadbr2.5kb$Identity = "tadbr"

#3. Removing the TAD boundaries from TAD domains
tmp  = tad
tmp$V2 = tmp$V2+2500
tmp$V3 = tmp$V3-2500
tmp$Width = tmp$V3-tmp$V2
tmp = tmp[which(tmp$Width > 0),] #making sure that all TAds are still greater than 0 in length
tmp = tmp[,c(1:3,5,4)]
names(tmp) = c("Chr","Start","Stop","Width","ID")
tmp$Identity = "TAD"

tad.br.nontad = rbind(tmp, tadbr2.5kb, nontad)
write.table(tad.br.nontad, "tad.nontad.br.bed", quote = F, row.names = F, col.names = T, sep="\t")


library(data.table)
methyl = fread("TADs/methylation_from_zoe.bedgraph", header = F, sep="\t")
tmp = tad.br.nontad
meanMetNormLen = c()
meanMet = c()
for (i in 1:nrow(tmp)){
  t = methyl[which(methyl$V1 == tmp$Chr[i] & methyl$V2 >= tmp$Start[i] & methyl$V3 <= tmp$Stop[i])]
  t$V5 = t$V3-t$V2
  meanMet[i] = mean(t$V4)
  meanMetNormLen[i] = sum(t$V4)/sum(t$V5)
}
tmp$meanMet = meanMet
tmp$meanMetNormLen = meanMetNormLen
write.table(tmp, "met.tad.nontad.br.txt", quote = F, row.names = F, col.names = T, sep="\t")


phast = fread("TADs/PhastCons8wayChrAll.bedgraph", header = F, sep="\t", colClasses=c("V2" = "numeric"))
phast$V2 = as.integer(phast$V2)
meanPhastNormLen = c()
meanPhast = c()
for (i in 1:nrow(tmp)){
  t = phast[which(phast$V1 == tmp$Chr[i] & phast$V2 >= tmp$Start[i] & phast$V3 <= tmp$Stop[i])]
  t$V5 = t$V3-t$V2
  meanPhast[i] = mean(t$V4)
  meanPhastNormLen[i] = sum(t$V4)/sum(t$V5)
}
tmp$meanPhast = meanPhast
tmp$meanPhastNormLen = meanPhastNormLen
write.table(tmp, "met_phast.tad.nontad.br.txt", quote = F, row.names = F, col.names = T, sep="\t")


fit = fread("TADs/rho_from_zoe.bedgraph", header = F, sep="\t")
meanRhoNormLen = c()
meanRho = c()
for (i in 1:nrow(tmp)){
  t = fit[which(fit$V1 == tmp$Chr[i] & fit$V2 >= tmp$Start[i] & fit$V3 <= tmp$Stop[i])]
  t$V5 = t$V3-t$V2
  meanRho[i] = mean(t$V4)
  meanRhoNormLen[i] = sum(t$V4)/sum(t$V5)
}
tmp$meanRho = meanRho
tmp$meanRhoNormLen = meanRhoNormLen

write.table(tmp, "met_phast_rho.tad.nontad.br.txt", quote = F, row.names = F, col.names = T, sep="\t")





