# Loading packages
library(sva)
library(Biobase)
library(limma)
library(pheatmap)
library(ggplot2)
library(reshape2)
library(plyr)
library(data.table)
library(ggpubr)
library(RColorBrewer)
theme_set(theme_pubr())


# Loading Data Matrices of Normalized Expression Values
setwd("/Users/sayeh/Documents/uni/pRact/ThesisCode/data/NormalizedDataRMA/")
Data1 <- read.delim("GSE6575NormalMatrix.txt", sep = " ")  #GPL570
Data2 <- read.delim("GSE89594NormalMatrix.txt", sep = " ") #GPL16699** , one color, GeneSymbol in End
Data3 <- read.delim("GSE26415NormalMatrix.txt", sep = " ")  #GPL6480 , one color, GeneSymbol in End
Data4 <- read.delim("GSE18123NormalMatrix.txt", sep = " ") #GPL6244
Datas <- list(Data1[ ,-c(1,2)], Data2[,-95], Data3[,-85], Data4[,-c(1,2)])
# plotting boxes
pdf("/Users/sayeh/Documents/uni/pRact/ThesisCode/data/NormalizedDataRMA/BoxDatas.pdf")
par(mfrow=c(2,2))
boxplotDatas <- lapply(Datas, boxplot)
dev.off()
## Data 1 & 4 both are affymetrix with same summary report of quantiles as for 2,3 agilent arrays
colnames(Data1)[2] <- "GeneSymbol"
colnames(Data2)[95] <- "GeneSymbol"
colnames(Data3)[85] <- "GeneSymbol"
colnames(Data4)[2] <- "GeneSymbol"
# removing gene symbols with more than one names "///" $ excess cols
Data1$`GeneSymbol` <- sub("/.*", "", Data1$`GeneSymbol`)
Data4$`GeneSymbol` <- sub("/.*", "", Data4$`GeneSymbol`)
Data1 <- Data1[,-1]
Data4 <- Data4[,-1]
# Aggregate gene symbol values finding means for each gene symbol with many probe IDs
Data1 <- aggregate(.~GeneSymbol, Data1, mean)
Data2 <- aggregate(.~GeneSymbol, Data2, mean)
Data3 <- aggregate(.~GeneSymbol, Data3, mean)
Data4 <- aggregate(.~GeneSymbol, Data4, mean)
# Merge all datasets 
all <- Reduce(merge,list(Data1, Data2, Data3, Data4, by = "GeneSymbol"))

all <- subset(all, GeneSymbol !="")
rownames(all) <- all$GeneSymbol
all <- subset(all, select=-GeneSymbol)
all <- all[,-length(all)] # an excess col I dont know where it came from!

# removing Developmentally delayed cases


setwd("/Users/sayeh/Documents/uni/pRact/ThesisCode/")
grp <- read.delim("all.txt") 

colnames(all) <-grp$samples
grp <- subset(grp, meta.data != "")
all <- all[,grp$samples]

#write.csv(all, file = "/Users/sayeh/Documents/uni/pRact/ThesisCode/allDDremoved.csv")

#pdf("/Users/sayeh/Documents/uni/pRact/ThesisCode/data/Plots/allBoxplot.pdf")
#boxAll <- boxplot(all)
#dev.off()

############################################### Batch Effect Removal ##########################################################################

batch <- factor(c(rep(1,47), rep(2,62), rep(3,42), rep(4,138)))
# batch effect removal
allc <- ComBat(all, batch)
allc <- as.data.frame(allc)


pdf("/Users/sayeh/Documents/uni/pRact/ThesisCode/data/Plots/allBefore&AfterBatchBoxplot.pdf")
par(mfrow= c(2,1))
ggplot(stack(all), aes(x = ind , y = values)) +
  geom_boxplot(fill = batch) + labs(x = "Samples", y = "log2Expression Values") + scale_fill_brewer(palette="BuPu")

ggplot(stack(allc), aes(x = ind , y = values)) +
  geom_boxplot(fill=batch) + labs(x = "Samples", y = "log2Expression Values", subtitle="After Batch Removed") + scale_fill_brewer(palette="Dark2")
dev.off()

################################################### PCA analysis ##################################################################################

Group <- (grp[,2])
Group2 <- (grp[,3])
# PCA for before batch removal
allm0 <- all - rowMeans(all)
pc0 <- prcomp(allm0) 
pcr0 <- data.frame(pc0$rotation[,1:3], batch,Group)
# PCA after batch removal
allm <- allc - rowMeans(allc) #subtracts mean of every row(gene expression mean) from every sample(columns) to only witness differences
pc <- prcomp(allm)
pcr <- data.frame(pc$rotation[,1:3], batch,Group)

# Compare PCA before and after batch removal 
pdf("/Users/sayeh/Documents/uni/pRact/ThesisCode/data/Plots/PCAallBefore&afterBatch.pdf")
par(mfrow = c(2,2))
PCA <- ggplot(pcr, aes(PC1, PC2, color=Group)) + geom_point(size=3) + theme_bw()
PCA0 <- ggplot(pcr0, aes(PC1, PC2, color=Group)) + geom_point(size=3) + theme_bw()
dev.off()
PCAplot <- ggarrange(PCA, PCA0,
                    labels = c("After Batch Removal", "Before Batch Removal"),
                    ncol = 1, nrow = 2)


#pheatmap(cor(allc))
#labels_row = Group, labels_col = Group, row

########################################### Differential Expression Analysis ################################################################

cl = factor(Group)  
groups <- levels(cl)
designMatrix=model.matrix(~cl + 0, allc)
colnames(designMatrix) <- levels(cl)
fit <- lmFit(allc, designMatrix)  # fit linear model
cts <-  c("ASDfemale-ASDmale", "ASDmale-TDmale", "TDfemale-TDmale", "ASDfemale-TDfemale") 
cont.matrix <- makeContrasts(contrasts=cts, levels=designMatrix)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)
dT <- decideTests(fit2, adjust.method="fdr", p.value=0.05)
DecideResults <- as.data.frame(summary(dT))
write.fit(fit2, results=dT, adjust="BH", method="global", file="/Users/sayeh/Documents/uni/pRact/ThesisCode/IntegratedMatricesDEresults.txt")

########################################### Up-regulated/Down-regulated Genes ###############################################

ASDf.m.Up=tT[which(tT$ASDfemale.ASDmale > 1 & tT$P.Value < 0.05),]
ASDTDDown=tT[which(tT$ASDfemale.ASDmale < -1 & tT$P.Value< 0.05),]

TD.f.m.Up = tT[which(tT$TDfemale.TDmale > 1.5 & tT$P.Value < 0.05),]
TD.f.m.Down =tT[which(tT$TDfemale.TDmale < -1.5 & tT$P.Value < 0.05),]
write.table(ASD.f.m.Down, file = "/Users/sayeh/Documents/uni/pRact/ThesisCode/Results/ASD.f.m.Down.csv", sep = ",", row.names = TRUE, col.names = TRUE)


###################################################### Visualization ##############################################

colnames(fit2) # list contrast names
par(mfrow=c(2,2))
ct <- 4        # choose contrast of interest
volcanoplot(fit2, coef=ct, main=colnames(fit2)[ct], pch=20,
            highlight=length(which(dT[,ct]!=0)), names=rep('+', nrow(fit2)))
abline(v=c(-1.5,1.5))

# MD plot (log fold change vs mean log expression)
plotMD(fit2, column=ct, status=dT[,ct], legend=F, pch=20, cex=1)




