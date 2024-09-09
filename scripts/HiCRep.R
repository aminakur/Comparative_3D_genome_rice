setwd("~/Downloads/Comparative Hi-C project/hic maps and genome tracks")

#install devtools

# install hicrep
devtools::install_github("qunhualilab/hicrep")
library(hicrep)

# To work with .hic file, the package 'strawr' is needed. Install it:
remotes::install_github("aidenlab/straw/R")
library(strawr)

#hic2mat() function directly read in .hic format file and convert it to squared matrix
mat1 <- hic2mat("azucena-1.hic", chromosome1 = "CM020634.1", chromosome2 = "CM020634.1", resol =   50000, method = "NONE")
mat2 <- hic2mat("azucena-2.hic", chromosome1 = "CM020634.1", chromosome2 = "CM020634.1", resol =   50000, method = "NONE")

#train the smoothing parameter h
h_x <- htrain(mat1, mat2, 50000, 1000000, 0:150)

# Compute SCC for one chromosome
scc.out = get.scc(mat1, mat2, resol = 5000, h = 30, lbr = 0, ubr = 1000000)
scc.out

# Compute all chromosomes (CM020876.1 format)
all.scc <- list()
for (i in paste0("HG4171", c(as.character(65.1:76.1)))) {
  mat1.chr <- hic2mat("ruf-1.hic", chromosome1 = i, chromosome2 = i, resol =   50000, method = "NONE")
  mat2.chr <- hic2mat("ruf-2.hic", chromosome1 = i, chromosome2 = i, resol =   50000, method = "NONE")
  all.scc[[i]] = get.scc(mat1.chr, mat2.chr, resol = 5000, h = 30, lbr = 0, ubr = 1000000)
}

#mean SCC for all chromosomes
meanSCC = (all.scc[["HG417165.1"]][["scc"]]+all.scc[["HG417166.1"]][["scc"]]+all.scc[["HG417167.1"]][["scc"]]+all.scc[["HG417168.1"]][["scc"]]+all.scc[["HG417169.1"]][["scc"]]+all.scc[["HG417170.2"]][["scc"]]+all.scc[["HG417171.2"]][["scc"]]+all.scc[["HG417172.1"]][["scc"]]+all.scc[["HG417173.1"]][["scc"]]+all.scc[["HG417174.1"]][["scc"]]+all.scc[["HG417175.1"]][["scc"]]+all.scc[["HG417176.1"]][["scc"]])/12
meanSCC

# Compute all chromosomes (chr01-chr12 format)
all.scc <- list()
for (i in paste0("chr", sprintf("%02d", (1:12)))) {
  mat1.chr <- hic2mat("NPB-43.hic", chromosome1 = i, chromosome2 = i, resol =   50000, method = "NONE")
  mat2.chr <- hic2mat("NPB-2.hic", chromosome1 = i, chromosome2 = i, resol =   50000, method = "NONE")
  all.scc[[i]] = get.scc(mat1.chr, mat2.chr, resol = 50000, h = 30, lbr = 0, ubr = 1000000)
}

#Compute SCC for all O.rufipogon chromosomes
ruf_names <- c("HG417165.1","HG417166.1","HG417167.1","HG417168.1","HG417169.1","HG417170.2","HG417171.2","HG417172.1","HG417173.1","HG417174.1","HG417175.1","HG417176.1")

all.scc <- list()
for (i in ruf_names) {
  mat1.chr <- hic2mat("ruf-1.hic", chromosome1 = i, chromosome2 = i, resol =   50000, method = "NONE")
  mat2.chr <- hic2mat("ruf-2.hic", chromosome1 = i, chromosome2 = i, resol =   50000, method = "NONE")
  all.scc[[i]] = get.scc(mat1.chr, mat2.chr, resol = 5000, h = 100, lbr = 0, ubr = 1000000)
}

#mean SCC for all chromosomes
meanSCC = (all.scc[["HG417165.1"]][["scc"]]+all.scc[["HG417166.1"]][["scc"]]+all.scc[["HG417167.1"]][["scc"]]+all.scc[["HG417168.1"]][["scc"]]+all.scc[["HG417169.1"]][["scc"]]+all.scc[["HG417170.2"]][["scc"]]+all.scc[["HG417171.2"]][["scc"]]+all.scc[["HG417172.1"]][["scc"]]+all.scc[["HG417173.1"]][["scc"]]+all.scc[["HG417174.1"]][["scc"]]+all.scc[["HG417175.1"]][["scc"]]+all.scc[["HG417176.1"]][["scc"]])/12
meanSCC

#for two different varieties (NPB and Az in this example)
matNPB_names <- c("chr01","chr02","chr03","chr04","chr05","chr06","chr07","chr08","chr09","chr10","chr11","chr12")
matAz_names <- c("CM020633.1", "CM020634.1","CM020635.1","CM020636.1","CM020637.1","CM020638.1","CM020639.1","CM020640.1","CM020641.1","CM020642.1","CM020643.1","CM020644.1")

all.scc <- list()
for(i in 1:12){
  chrname = matNPB_names[i]
  mat1.chr <- hic2mat("NPB-2.hic", chromosome1 = chrname, chromosome2 = chrname, resol =   50000, method = "NONE")
  chrname1 = matAz_names[i]
  mat2.chr <- hic2mat("azucena-1.hic", chromosome1 = chrname1, chromosome2 = chrname1, resol =   50000, method = "NONE")
  all.scc[[i]] = get.scc(mat1.chr, mat2.chr, resol = 5000, h = 30, lbr = 0, ubr = 1000000)
}

#Equalize the total number of reads
data(mat1)
#check total number of reads before adjustment
sum(mat1[,-c(1:3)])
DS_mat1 <- depth.adj(mat1, 5782958) 
#check total number of reads after adjustment
sum(DS_mat1[,-c(1:3)])















