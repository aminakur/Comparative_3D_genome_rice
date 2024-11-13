### This script will do multiple things
#1. intersect the genome co-ordinates with TADs to get the non-TAD bed
#2. Get gene density for plotting - Amina (but she had already done this, so not needed)
#3. identify genes that are within the TAD (within TAD domains)
#4. identify genes that are within the nonTAD domains 
#5. identify genes within at TAD boundaries (TAD start/stop +/- 2.5kb) and not within boundaries
#6. identify all pairs of neighboring genes, and then classify them as separated by TAD boundary vs not
#7. Classify the genome as tad boundary, tad domain (minus the boundary) and nontad domains, and identify genes within them

setwd("~/Desktop/NYU_Projects/Amina_GeneExp/data/")
library(GenomicRanges)

#1. intersect the genome co-ordinates with TADs to get the non-TAD bed
gen = read.table("AzucenaRef/AzuChr.bed", header = F, sep="\t")
tad = read.table("TADs/TADs_5kb_intersect80.bed", header = F, sep = "\t")
gengr = GRanges(seqnames=gen$V1,
             IRanges(start=gen$V2,end=gen$V3))
tadgr = GRanges(seqnames=tad$V1,
                IRanges(start=tad$V2,end=tad$V3), tadid = tad$V4)
#hits <- findOverlaps(tadgr, gengr, type="within")
nontadgr = GenomicRanges::setdiff(gengr, tadgr)
nontad = as(nontadgr, "data.frame")
t = c(1:nrow(nontad))
tmp = paste("NT", t, sep = "")
nontad$ID = tmp
remove(t); remove(tmp)

str(nontad)
sum(nontad$width); summary(nontad$width)
sum(tad$V3-tad$V2); summary(tad$V3-tad$V2)

nontad = nontad[,-5]
write.table(nontad, "TADs/nonTADs_5kb_intersect80.bed", quote = F, row.names = F, col.names = T, sep="\t")

##################################################################################################################################
#2. Getting gene density -- proportion of genes every 1kb thorughout the genome
gen = read.table("AzucenaRef/AzuChr.bed", header = F, sep="\t")
seqlengths = gen$V3
names(seqlengths) = gen$V1
genTile = tileGenome(seqlengths, tilewidth = 1000, cut.last.tile.in.chrom = T)

gene = read.table("AzucenaRef/OSJAZ_OSNIP.gene.bed", header = F, sep = "\t")
genegr = GRanges(seqnames=gene$V1,
                 IRanges(start=gene$V2,end=gene$V3), geneid = gene$V5)

hits = GenomicRanges::findOverlaps(genegr, genTile, type="any")
overlaps <- pintersect(genegr[queryHits(hits)], genTile[subjectHits(hits)])
percentOverlap <- width(overlaps) / width(genTile[subjectHits(hits)])

tmp = as.data.frame(genTile)
tmp2 =  as.data.frame(genTile[subjectHits(hits)])
tmp2$GenDenperKb = percentOverlap 
genDen = merge(tmp[,c(1:3)], tmp2[,c(1:3,6)], by=c("seqnames", "start", "end"), all.x = T)
str(genDen)
write.table(genDen, "AzucenaRef/GeneDesityPerKb.txt", row.names = F, col.names = F, quote = F, sep="\t")


##################################################################################################################################
#3. identify genes that are within the TAD (within TAD domains)
###Now I will identify genes that are within the TAD  -- completely within (not on the boundary)
library(stringr)
library(dplyr)
gene = read.table("AzucenaRef/OSJAZ_OSNIP.gene.bed", header = F, sep = "\t")
genegr = GRanges(seqnames=gene$V1,
                 IRanges(start=gene$V2,end=gene$V3), geneid = gene$V5)

hits = findOverlaps(genegr, tadgr, type="within")
df_ann <- cbind(tad[subjectHits(hits),],gene[queryHits(hits),])
names(df_ann) = c("Chr", "tadStart", "tadStop", "tadID", "geneChr", "geneStart", "geneStop", "geneStrand", "geneID")
tmp = df_ann %>%
  group_by(Chr, tadStart, tadStop, tadID) %>%
  summarize(genes = paste(geneID, collapse = ','))      #To collapse genes in one row
tmp2 = df_ann %>%
  group_by(Chr, tadStart, tadStop, tadID) %>% tally()   #To count the number of genes
            
t = merge(tmp, tmp2, by=c("Chr", "tadStart", "tadStop", "tadID"))
write.table(t, "TADs/TADs_5kb_intersect80.genes.allwithin", quote = F, col.names = T, row.names = F, sep = "\t")


##################################################################################################################################
#4. identify genes that are within the nonTAD domains 
###Now I will identify genes that are within the nonTAD domains -- completely within (not on the boundary)
hits = findOverlaps(genegr, nontadgr, type="within")
df_ann <- cbind(nontad[subjectHits(hits),],gene[queryHits(hits),])
df_ann = df_ann[,-4]
names(df_ann) = c("Chr", "nontadStart", "nontadStop", "nontadID", "geneChr", "geneStart", "geneStop","geneStrand", "geneID")
tmp = df_ann %>%
  group_by(Chr, nontadStart, nontadStop, nontadID) %>%
  summarize(genes = paste(geneID, collapse = ','))      #To collapse genes in one row
tmp2 = df_ann %>%
  group_by(Chr, nontadStart, nontadStop, nontadID) %>% tally()   #To count the number of genes

nt = merge(tmp, tmp2, by=c("Chr", "nontadStart", "nontadStop", "nontadID"))
write.table(nt, "TADs/nonTADs_5kb_intersect80.genes.allwithin", quote = F, col.names = T, row.names = F, sep = "\t")

names(t) = c("Chr", "Start", "Stop", "ID", "geneName", "freqGenes")
names(nt) = c("Chr", "Start", "Stop", "ID", "geneName", "freqGenes")
tmp = rbind(t, nt)
write.table(tmp, "TADs/TAD-nonTAD.genes.allwithin", quote = F, row.names = F, col.names = T, sep = "\t")


##################################################################################################################################
#5. identify genes within at TAD boundaries (TAD start/stop +/- 2.5kb) and not within boundaries
#making the tad boundary dataframe
tadbr = tad[,c(1,2,4)]
tmp = tad[,c(1,3,4)]
colnames(tadbr) = c("Chr", "Pos", "ID")
colnames(tmp) = c("Chr", "Pos", "ID")
t = rbind(tadbr, tmp)
tadbr = t
#So there are a total of 2414 TAD boundaries (for 1207 tads), but some of these boundaries are replicated, so I remove these
#Removing Duplicates (1818 boundaries remain)
t = unique(tadbr[,c(1,2)])

#Now I will identify the gene given within boundary 
#boundary is defined as the tad domain boundary +/- 2500, 
#Genes will have 1000bp of promoter region included
tadbr2.5kb = t
tadbr2.5kb$Start = t$Pos - 2500; tadbr2.5kb$Stop = t$Pos + 2500;
tadbr2.5kb = tadbr2.5kb[,-2]
tadbr2.5kbgr = GRanges(seqnames=tadbr2.5kb$Chr,
                  IRanges(start=tadbr2.5kb$Start,end=tadbr2.5kb$Stop))
tmp = gene
tmp$V2 = ifelse(tmp$V4 == '+', tmp$V2-1000, tmp$V2)
tmp$V3 = ifelse(tmp$V4 == '-', tmp$V3+1000, tmp$V3)
geneprgr = GRanges(seqnames=tmp$V1,
                 IRanges(start=tmp$V2,end=tmp$V3), geneid = tmp$V5, strand = tmp$V4)
hits = findOverlaps(geneprgr, tadbr2.5kbgr, type="any")
df_ann <- cbind(tadbr2.5kb[subjectHits(hits),],tmp[queryHits(hits),])
names(df_ann) = c("Chr", "tadbrStart", "tadbrStop", "geneChr", "geneStart", "geneStop", "geneStrand", "geneID")

tmp = setdiff(gene$V5, df_ann$geneID)
write.table(df_ann, "TADs/TADboundarys_5kb_intersect80.genes.any", quote = F, col.names = T, row.names = F, sep = "\t")
write.table(tmp, "TADs/notTADboundarys_5kb_intersect80.genes.any", quote = F, col.names = T, row.names = F, sep = "\t")


##################################################################################################################################
#6. Next thing I want to identify all pairs of neighboring genes, and then classify them as separated by TAD boundary vs not
#install.packages("devtools")
#devtools::install_github("pgpmartin/GeneNeighborhood", force = TRUE)
library(GeneNeighborhood)
gene$Gnumber = c(1:nrow(gene))
genegr = GRanges(seqnames=gene$V1,
                 IRanges(start=gene$V2,end=gene$V3, names = gene$Gnumber), strand = gene$V4)

GeneNeighbors <- getGeneNeighborhood(genegr)
#The above function idenitfies genes lying before/after with/without strand info, and identifies overlapping gene
table(GeneNeighbors$GenePair)
#So we have a total of 9621 H2H (divergent), 9187 T2T (convergent) and 17760 T2H (tandem)
##NOTE:: in the manual it says this is the gene in the follow column, but I checked manually and it is with gene in the "precede" col

pairs = GeneNeighbors %>% filter(GenePair == "H2H" | GenePair == "T2T" | GenePair == "T2H")
pairs = as.data.frame(pairs[,c(1,4,20)]); names(pairs)[2] = "PairName"
pairs$GeneName = as.integer(pairs$GeneName); pairs$PairName = as.integer(pairs$PairName)
tmp = gene; 
names(tmp)[6] = "GeneName"
tmp1 = merge(pairs, tmp, by = "GeneName")
names(tmp)[6] = "PairName"
tmp2 = merge(pairs, tmp, by = "PairName")
tmp = cbind(tmp1, tmp2)
tmp = tmp[,c(8,16,4,5,6,7,12,13,14,15,3)]
names(tmp) = c("Gene1","Gene2",
               "G1_Chr","G1_Start","G1_Stop","G1_Strand",
               "G2_Chr","G2_Start","G2_Stop","G2_Strand",
               "GenePair")
pairs = tmp
#This has all the pairs in the genome classified (total 36k pairs)
write.table(pairs, "~/Desktop/NYU_Projects/Amina_GeneExp/data/AzucenaRef/GenePairs_all.txt", quote = F, row.names = F, col.names = T, sep = "\t")

##Next I will determine which of these pairs have Tad boubdaries between them 
#For this I will use G1_Start and G2_Stop
library(stringr)
#making the tad boundary dataframe
tadbr = tad[,c(1,2,4)]
tmp = tad[,c(1,3,4)]
colnames(tadbr) = c("Chr", "Pos", "ID")
colnames(tmp) = c("Chr", "Pos", "ID")
t = rbind(tadbr, tmp)
tadbr = t
#So there are a total of 2414 TAD boundaries (for 1207 tads), but some of these boundaries are replicated, so I remove these
#Removing Duplicates (1818 boundaries remain)
t = unique(tadbr[,c(1,2)])

#Now I will identify the before and after gene given the boundary
tadbrgr = GRanges(seqnames=t$Chr,
                IRanges(start=t$Pos,end=t$Pos))
pairgr = GRanges(seqnames = pairs$G1_Chr,
                 IRanges(start=pairs$G1_Start,end=pairs$G2_Stop))

hits = findOverlaps(pairgr, tadbrgr, type=c("any"))
df_ann <- cbind(t[subjectHits(hits),],pairs[queryHits(hits),])
#there are 2347 unique gene pairs that are separated by TAD boundaries (10 pairs are replicated across boundaries)

tadbrgenes = df_ann
nontadbrgenes = setdiff(pairs, tadbrgenes[,-c(1,2)])
#33861 pairs

write.table(tadbrgenes, "TADs/TADs_5kb_intersect80.boundarysep.genes", quote = F, row.names = F, col.names = T, sep = "\t")
write.table(nontadbrgenes, "TADs/5kb_80_nonboundarysep.genes", quote = F, row.names = F, col.names = T, sep = "\t")

dim(tadbrgenes[which(tadbrgenes$Gene1 %in% Wgenexp$geneID & tadbrgenes$Gene2 %in% Wgenexp$geneID),])
dim(nontadbrgenes[which(nontadbrgenes$Gene1 %in% Wgenexp$geneID & nontadbrgenes$Gene2 %in% Wgenexp$geneID),])

#divergent = GeneNeighbors %>% filter(GenePair == "H2H")
#convergent = GeneNeighbors %>% filter(GenePair == "T2T")
#tandem = GeneNeighbors %>% filter(GenePair == "T2H")


##################################################################################################################################
#7. Classify the genome as tad boundary, tad domain (minus the boundary) and nontad domains, and identify genes within them

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


tad.br.nontad.gr = GRanges(seqnames = tad.br.nontad$Chr,
                           IRanges(start = tad.br.nontad$Start, end = tad.br.nontad$Stop), ID = tad.br.nontad$ID)
gene = read.table("AzucenaRef/OSJAZ_OSNIP.gene.bed", header = F, sep = "\t")

genegr = GRanges(seqnames=gene$V1,
                 IRanges(start=gene$V2,end=gene$V3), geneid = gene$V5)

hits = findOverlaps(genegr, tad.br.nontad.gr, type="within")
df_ann <- cbind(tad.br.nontad[subjectHits(hits),],gene[queryHits(hits),])
names(df_ann) = c("Chr", "Start", "Stop", "Width", "ID", "Identity", "geneChr", "geneStart", "geneStop", "geneStrand", "geneID")
tmp = df_ann %>%
  group_by(Chr, Start, Stop, Width, ID, Identity) %>%
  summarize(genes = paste(geneID, collapse = ','))      #To collapse genes in one row
tmp2 = df_ann %>%
  group_by(Chr, Start, Stop, Width, ID, Identity) %>% tally()   #To count the number of genes

t = merge(tmp, tmp2, by=c("Chr", "Start", "Stop", "Width", "ID", "Identity"))
write.table(t, "TADs/tad.nontad.br_5kb_intersect80.genes.allwithin", quote = F, col.names = T, row.names = F, sep = "\t")

table(tad.br.nontad$Identity)
table(t$Identity)

tmp = t[which(t$Identity == "TAD"),]
summary(tmp$n)
sum(tmp$n)


###identifying genes with only promoters (completely within/any overlap) within these three groups
tad.br.nontad = read.table("TADs/tad.nontad.br.bed", header = T, sep = "\t")
tad.br.nontad.gr = GRanges(seqnames = tad.br.nontad$Chr,
                           IRanges(start = tad.br.nontad$Start, end = tad.br.nontad$Stop), ID = tad.br.nontad$ID)
gene = read.table("AzucenaRef/OSJAZ_OSNIP.gene.bed", header = F, sep = "\t")
tmp = gene
tmp$V2 = ifelse(gene$V4 == '+', gene$V2-1000, gene$V3)
tmp$V3 = ifelse(gene$V4 == '+', gene$V2, gene$V3+1000)
geneprgr = GRanges(seqnames=tmp$V1,
                   IRanges(start=tmp$V2,end=tmp$V3), geneid = tmp$V5, strand = tmp$V4)

hits = findOverlaps(geneprgr, tad.br.nontad.gr, type="within")
df_ann <- cbind(tad.br.nontad[subjectHits(hits),],tmp[queryHits(hits),])
names(df_ann) = c(names(tad.br.nontad), "geneChr", "geneStart", "geneStop", "geneStrand", "geneID")
write.table(df_ann, "TADs/tad.nontad.br_5kb_intersect80.promoters.allwithin.txt", quote = F, col.names = T, row.names = F, sep = "\t")

hits = findOverlaps(geneprgr, tad.br.nontad.gr, type="any")
df_ann <- cbind(tad.br.nontad[subjectHits(hits),],tmp[queryHits(hits),])
names(df_ann) = c(names(tad.br.nontad), "geneChr", "geneStart", "geneStop", "geneStrand", "geneID")
write.table(df_ann, "TADs/tad.nontad.br_5kb_intersect80.promoters.anyoverlap.txt", quote = F, col.names = T, row.names = F, sep = "\t")

##################################################################################################################################

###Ignore below
#Now I will identify the before and after gene given the boundary
tadbrgr = GRanges(seqnames=tadbr$Chr,
                  IRanges(start=tadbr$Start,end=tadbr$Stop), tadid = tadbr$ID)
tmp1 = follow(tadbrgr, genegr, select=c("last"), ignore.strand=TRUE)
tmp2 = precede(tadbrgr, genegr, select=c("first"), ignore.strand=TRUE)
t1 = gene[tmp1,]
t2 = gene[tmp2,]
tadbr_genes = cbind(tadbr,t1,t2)
#But there also could be gene that intersects with boundary (lets call this the "at gene") -- these would not have been counted in before and after.
#In this case, before and at gene, and after and at gene are both cases that are separated by a tad boundary 
hits = findOverlaps(genegr, tadbrgr, type=c("any"))
df_ann <- cbind(tadbr[subjectHits(hits),],gene[queryHits(hits),])
names(df_ann) = c("Chr", "tadStart", "tadStop", "tadID", "geneChr", "geneStart", "geneStop", "geneID")

tmp = df_ann %>%
  group_by(Chr, tadStart, tadStop, tadID) %>%
  summarize(genes = paste(geneID, collapse = ','))      #To collapse genes in one row


w = (tad$V3-tad$V2)/1000
t = rep("TAD", times = length(w))
w = append(w, (nontad$width/1000))
t = append(t, rep("nonTAD", times = length(nontad$width)))
df = data.frame(w,t)
ggplot(df, aes(x=w, fill=t)) +
  geom_histogram(alpha=0.4)



