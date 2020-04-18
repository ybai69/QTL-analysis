setwd("~/Dropbox (UFL)/pch6088_proj_code")
rm(list=ls())



#Loading the packages
library("qtl") 
library("ggplot2")
library("bestNormalize")
library(dplyr)

#Read the data from .csv files
geno <- read.csv(file = 'Clean/genotypes_clean.csv')
sex <- read.csv(file = 'Clean/gender.csv')
islet <- read.csv(file = 'Clean/islet_mlratio_clean.csv')
dim(islet)
dim(geno)
dim(sex)
annot <- read.csv(file = 'Clean/microarray_annot.csv')
sum(is.na(annot$pos.cM))
#[1] 30
a_gene_id.new <-annot$a_gene_id[!(is.na(annot$pos.cM))]
a_gene_id.new <- paste0("X",a_gene_id.new)


batch.islet <- read.csv(file = 'Clean/islet_batch_clean.csv')
batch.islet.sorted <- batch.islet[order(batch.islet$mouseNum),]
MouseNum.batch.islet <- batch.islet.sorted$mouseNum
batch.islet.idx <-MouseNum.batch.islet[! MouseNum.batch.islet %in% MouseNum.intersect ]
batch.islet.sorted <- batch.islet.sorted [ !(batch.islet.sorted$mouseNum %in% batch.islet.idx ), ]
write.csv(batch.islet.sorted,'islet/batch_islet_sorted.csv')


#Clean the data
MouseNum.islet <- islet$MouseNum
MouseNum.geno <- sex$MouseNum

MouseNum.intersect <- intersect(MouseNum.islet,MouseNum.geno )
length(MouseNum.intersect)
geno.idx <- MouseNum.geno[! MouseNum.geno %in% MouseNum.intersect ]
islet.idx <-MouseNum.islet[! MouseNum.islet %in% MouseNum.intersect ]

geno.new <- geno[ !(geno$MouseNum %in% geno.idx ), ]
dim(geno.new)
geno.new.sorted <- geno.new[order(geno.new$MouseNum),]
sex.new <- sex[!(sex$MouseNum %in% geno.idx ),]
islet.new <- islet[ !(islet$MouseNum %in% islet.idx ), ]
dim(islet.new)

write.csv(islet.new,'islet_new.csv')
write.csv(geno.new,'geno_new.csv')
write.csv(sex.new,'sex_new.csv')
write.csv(geno.new.sorted,'islet/geno_new_sorted.csv', row.names=FALSE)
islet.sex.new <- merge(islet.new, sex.new,  by.x="MouseNum", by.y="MouseNum", 
                       sort = FALSE)
islet.sex.new.sorted <- merge(islet.new, sex.new, 
                              by.x="MouseNum", by.y="MouseNum")
write.csv(islet.sex.new.sorted,'islet/pheno_new_sorted.csv', row.names=FALSE)

islet.colnames <- colnames(islet.sex.new.sorted)[-c(1,40574,40575)]
length(islet.colnames)
length(islet.colnames[islet.colnames %in% a_gene_id.new ])
#[1] 37796
islet.colnames.new <- islet.colnames[islet.colnames %in% a_gene_id.new ]

islet.sex.new.sorted.sub <-subset(islet.sex.new.sorted, 
                                  select = c("MouseNum",islet.colnames.new, 
                                             "Sex" ,"pgm"))
write.csv(islet.sex.new.sorted.sub,'islet/pheno_new_sorted_sub.csv', 
          row.names=FALSE)



#read data for a QTL experiment
dat <- read.cross("csvs", dir="islet", genotypes=c("RR","BR","BB"), 
                  alleles=c("R","B"), genfile="geno_new_sorted.csv", 
                  phefile="pheno_new_sorted_sub.csv")


#test <- read.cross("csvs", genfile="listeria_gen.csv", phefile="listeria_phe.csv")

#########summary of cross data
summary(dat)
nind(dat)
nphe(dat)
nchr(dat)
totmar(dat)
nmar(dat)
plotMissing(dat, chr= c(1:19,"X"))
plotMap(dat, chr= c(1:19,"X"))


pheno_col <- 2
plotPheno(dat, pheno.col=pheno_col, xlab =paste("phe", pheno_col-1))

pheno_col <- 201
plotPheno(dat, pheno.col=pheno_col, xlab =paste("phe", pheno_col-1))

pheno_col <- 2001
plotPheno(dat, pheno.col=pheno_col, xlab =paste("phe", pheno_col-1))


pheno_col <- 20001
plotPheno(dat, pheno.col=pheno_col, xlab =paste("phe", pheno_col-1))

pheno_col <- 30001
plotPheno(dat, pheno.col=pheno_col, xlab =paste("phe", pheno_col-1))


#Plot the genotype data
geno.image(dat,chr= c(1:19,"X"))



# copy old version before transformation
dat.old <- dat
#######normalization transformation 
orderNorm_obj <- apply(dat$pheno[,-c(1,37798,37799)],2, orderNorm)
for( i in 2:37797){
    dat$pheno[,i] <-  orderNorm_obj[[i-1]]$x.t
}
    
#check the distribution of pheno after transformation
pheno_col <- 2
plotPheno(dat, pheno.col=pheno_col, xlab =paste("phe", pheno_col-1))


pheno_col <- 21
plotPheno(dat, pheno.col=pheno_col, xlab =paste("phe", pheno_col-1))


pheno_col <- 201
plotPheno(dat, pheno.col=pheno_col, xlab =paste("phe", pheno_col-1))

pheno_col <- 2001
plotPheno(dat, pheno.col=pheno_col, xlab =paste("phe", pheno_col-1))


pheno_col <- 20001
plotPheno(dat, pheno.col=pheno_col, xlab =paste("phe", pheno_col-1))

pheno_col <- 30001
plotPheno(dat, pheno.col=pheno_col, xlab =paste("phe", pheno_col-1))



#estimated recombination fractions between markers
dat <- est.rf(dat)
plotRF(dat,chr= c(1:19,"X"), col.scheme=c("redblue"))





dat<- calc.genoprob(dat, step=1,map.function="c-f",error.prob=0.002)
# Genotype probabilities for EM and H-K
pheno_col <- 201
out.em <- scanone(dat, pheno.col=pheno_col,method="em")
out.hk <- scanone(dat, pheno.col=pheno_col, method="hk")
plot(out.em,out.hk, chr = c(1:19,"X") ,col = c("red", " blue"),main = paste("phe",pheno_col-1 ) )
# Add a legend
#legend(list(x = 1500,y = 3.2), legend=c("EM", "HK"), col=c("red", "blue"),lwd=c(2.5,2.5),)


# Genotype probabilities for EM and H-K
pheno_col <- c(1990)
out.em <- scanone(dat, pheno.col=pheno_col,method="em")
out.hk <- scanone(dat, pheno.col=pheno_col, method="hk")
plot(out.em,out.hk,chr = c(1:19,"X") ,col = c("red", " blue"),main = paste("phe",pheno_col-1 ) ,ylab="LOD score")


############################################
N <- 37796
lod.res <- vector("list", N);
pos.res <- vector("list",N);
flag <- rep(0,N);
pheno.idx <- which(flag < 1)
pheno.idx <- pheno.idx +1
for(pheno_col in pheno.idx ){
    print(pheno_col-1)
    out.hk <- scanone(dat, chr = c(1:19,"X") ,pheno.col=pheno_col, method="hk")
    mylist <- split(out.hk, out.hk$chr)
    lod.group.max <- c()
    pos.group.max <- c()
    for(i in c(1: length(unique(out.hk$chr))) ){
        lod.max.temp <- max(mylist[[i]]$lod)
        lod.group.max <- c(lod.group.max, lod.max.temp) 
        pos.group.max <- c(pos.group.max, mylist[[i]][which.max(mylist[[i]]$lod),2] )
    }
    lod.res[[pheno_col-1]] <- lod.group.max
    pos.res[[pheno_col-1]] <- pos.group.max
    flag[pheno_col-1] <- 1
    #plot(out.hk,col = c("red", " blue"), main = paste("phe",pheno_col-1 ) ,ylab="LOD score")
}



############################################
N <- 37796
lod.cov.res <- vector("list", N);
pos.cov.res <- vector("list",N);
flag <- rep(0,N);
pheno.idx <- which(flag < 1)
pheno.idx <- pheno.idx +1
sex <- as.numeric(pull.pheno(dat,"Sex") == "Male")
for(pheno_col in pheno.idx ){
    print(pheno_col-1)
    out.covar <- scanone(dat, chr = c(1:19,"X") ,pheno.col=pheno_col, 
                         method="hk", addcovar=sex, intcovar=sex)
    mylist <- split(out.covar, out.covar$chr)
    lod.group.max <- c()
    pos.group.max <- c()
    for(i in c(1: length(unique(out.covar$chr))) ){
        lod.max.temp <- max(mylist[[i]]$lod)
        lod.group.max <- c(lod.group.max, lod.max.temp) 
        pos.group.max <- c(pos.group.max, mylist[[i]][which.max(mylist[[i]]$lod),2] )
    }
    lod.cov.res[[pheno_col-1]] <- lod.group.max
    pos.cov.res[[pheno_col-1]] <- pos.group.max
    flag[pheno_col-1] <- 1
}

############################################
annot.sub <- subset(annot, select = c("a_gene_id","pos.cM","chr"))
annot.sub$a_gene_id <- paste0("X",annot.sub$a_gene_id)
a_gene_id.intersect <- intersect(annot.sub$a_gene_id, colnames(dat$pheno))
annot.sub <- annot.sub[annot.sub$a_gene_id %in% a_gene_id.intersect  , ]
dim(annot.sub)
rownames(annot.sub) <- c(1:N)

annot.sub$chr <- factor(annot.sub$chr, levels = c(1:19,"X"))
annot.sub.sorted <- annot.sub[order(annot.sub$chr,annot.sub$pos.cM),]
head(annot.sub.sorted)
tail(annot.sub.sorted)



data.full <- data.frame(vars.temp = c(),
                        pos.temp = c(),
                        probe.pos.temp =c(),
                        lod.temp = c(),
                        colvars.temp = c())
for(c in levels(annot.sub.sorted$chr)){
    annot.temp <-subset(annot.sub.sorted,annot.sub.sorted$chr==c )
    print(max(annot.temp$pos.cM))
    print(min(annot.temp$pos.cM))
    #pheno.idx <- as.numeric(rownames(annot.temp))
    pheno.names <- colnames(dat$pheno)[-c(1,37798,37799)]
    pheno.idx <- match(annot.temp$a_gene_id, pheno.names) 
    pheno.idx <- pheno.idx +1
    lod.temp <- c()
    pos.temp <- c()
    for(pheno_col in pheno.idx ){
        lod.temp <- c(lod.temp, lod.cov.res[[pheno_col-1]])
        pos.temp <- c(pos.temp, pos.cov.res[[pheno_col-1]]   )
    }
    chr.temp <- c(1:19,"X")
    chr.temp <- rep(chr.temp,times=length(pheno.idx))
    vars.temp <- factor(chr.temp, levels = c(1:19,"X"))
    data.temp = data.frame(vars.temp = vars.temp,
                           pos.temp = pos.temp,
                           probe.pos.temp =rep(annot.temp$pos.cM,each=20),
                           lod.temp = lod.temp,
                           colvars.temp = rep(c,length(lod.temp)))
    data.full <- rbind(data.full,data.temp)
}


data.full$colvars.temp_f <- factor(data.full$colvars.temp, levels = c("X",19:1))


threshold.val <- 10
normalized.lod.temp  <-(lod.temp-min(lod.temp) )/(max(lod.temp)-min(lod.temp))-0.1
p1 <- ggplot(subset(data.full,data.full$lod.temp >= threshold.val), 
             aes(pos.temp, probe.pos.temp, color=lod.temp)) + 
    geom_point(size = 0.5, alpha =normalized.lod.temp[data.full$lod.temp >= threshold.val]  ) 
p1 <- p1  +  theme_classic() +theme( axis.text.x=element_blank(),axis.text.y=element_blank(), 
                                     strip.background = element_blank(), 
                                   panel.margin.y = unit(0, "lines"), panel.margin.x = unit(0,"lines")) + 
    labs(x = "QTL position (cM)", y = "Probe position (cM)", colour = "LOD scores")
# # Specify the colors for low and high ends of gradient
p1 <- p1 + scale_color_gradient(low ="grey" , high = "black") + coord_fixed(ratio = 1)
#p1 <- p1 + scale_colour_gradientn(colours = c("red","blue","green","black"), values = c(10,25,50,100)) 
p1 <- p1 + facet_grid( cols = vars(vars.temp) , rows = vars(colvars.temp_f))
p1


##############################################
# subset data that LOD score >= 10, 100
data.sub.10 <- subset(data.full, lod.temp>= 10)
data.sub.100 <- subset(data.full, lod.temp>= 100)
dim(data.sub.10)
dim(data.sub.100)
subset(data.sub.100, vars.temp == 6)

#Find marker closest to a specified position
mar6 <- find.marker(dat, 6, 73.806)

phe6 <- subset(annot.sub.sorted, (chr == 6)& (pos.cM <91.862464)& (pos.cM > 91.862462))
phe.col <- (1:length(colnames(dat$pheno)) )[colnames(dat$pheno) ==  phe6$a_gene_id]
#effectplot(dat, pheno.col=phe.col, mname1=mar6, mname2="Sex", mark2=dat$pheno$Sex, geno2=c("F","M"))
plotPXG(dat, mar6,pheno.col=phe.col, ylab=sub("X","",phe6$a_gene_id))
out.covar <- scanone(dat, chr = c(1:19,"X") ,pheno.col=phe.col, method="hk", addcovar=sex, intcovar=sex)
plot(out.covar, col= "blue",main = sub("X","",phe6$a_gene_id) ,ylab="LOD score")

###############################################################
data.sub.25 <- subset(data.full, lod.temp>= 25)
dim(data.sub.25)
dim(subset(data.sub.25, vars.temp == 6 & colvars.temp != 6 ))


#out.np <- scanone(dat, pheno.col=pheno_col, model="np")
#plot(out.em,out.hk, out.np,chr = c(1:19,"X") ,col = c("red", " blue", "yellow"),main = paste("phe",pheno_col-1 ) ,ylab="LOD score")
# Add a legend
#legend(list(x = 1500,y = 50), legend=c("EM", "HK"), col=c("red", "blue"),lwd=c(2.5,2.5),)


#lod.group.max <- tapply(out.hk$lod, out.hk$chr, max)
#out.hk$pos[out.hk$lod == max(lod.group.max)]
#out.hk$pos[out.hk$lod == max(lod.group.max)]
#out.hk$pos[out.hk$lod >=5]




permo <- scanone(dat, pheno.col=c(2,20,200,2000,20000), n.perm=1000)
summary(permo, alpha=c(0.20, 0.05, 0.01))
# X497628 X497662 X498037 X501623 X10002918255
# 20%    3.15    3.11    3.13    3.22         3.17
# 5%     3.76    3.79    3.84    3.92         3.85
# 1%     4.39    4.45    4.69    4.63         4.82
