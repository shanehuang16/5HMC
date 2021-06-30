library(data.table)
library(readxl)
library(DESeq2)
library(apeglm)
library(tidyverse)
library(VennDiagram)

setwd("C:/Liang Wang/gene level data")
Final <- read_excel("temporal samples.xlsx", sheet = "Final")

##########################
## Unpaired DESeq Tests ##
##########################

## Time 1 (A) ----

# Load in data and remove bad samples
data <- fread("Qianxia_raw_counts_by_Genes_FINAL.csv")
DESeq.data.1 <- data
DESeq.data.1 <- DESeq.data.1[, -c(2:6)]
DESeq.data.1 <- DESeq.data.1[, c(1:56)]

# Set up metadata/design matrix
metaData.1 <- data.frame(Final[1:56,c(2,4)])
metaData.1 <- metaData.1[-41,]  # A46 at row 41
colnames(metaData.1)[2] <- "Pheno"

# Create DESeq Data format
dds.1 <- DESeqDataSetFromMatrix(countData=DESeq.data.1, 
                              colData=metaData.1, 
                              design=~Pheno, tidy = TRUE)

# Run DESeq
dds.1 <- DESeq(dds.1)

# Pull up names for shrinkage
#resultsNames(dds.1)

# Run shrinkage to get final LFC results
resLFC.1 <- lfcShrink(dds.1, coef="Pheno_Y_vs_N", type="apeglm") # adaptive t prior shrinkage estimator

# LFC results
head(resLFC.1)
summary(resLFC.1)

# Write out results
write.csv(resLFC.1,"Additional Files/Time 1 DESeq2 Results.csv")

# Look at top genes with lowest adj. p-values
res.temp.1 <- resLFC.1[order(resLFC.1$padj),]
head(res.temp.1)

# Look at specific top genes
plotCounts(dds.1, gene="CYSLTR2", intgroup="Pheno")
plotCounts(dds.1, gene="CCDC85C", intgroup="Pheno")
plotCounts(dds.1, gene="CASZ1", intgroup="Pheno")

# Variance stabilization before PCA plot
vsdata.1 <- vst(dds.1, blind=FALSE)

# PCA plot
a.1=plotPCA(vsdata.1, intgroup="Pheno")
a.1
a.1$data
hist(resLFC.1$pvalue)

plot(a.1$data[,1],a.1$data[,2])
identify(a.1$data[,1],a.1$data[,2], a.1$data[,5], cex=0.8)

# Normalize data to write out
normal_RC.1 <- as.data.frame(counts(dds.1, normalized = T))
write.csv(normal_RC.1, "Additional Files/Time 1 Normalized.csv")




## Time 2 (B) ----

DESeq.data.2 <- data

metaData.2 <- data.frame(Final[57:111,c(2,4)])
colnames(metaData.2)[2] <- "Pheno"
month3.55 <- match(metaData.2[,1], colnames(DESeq.data.2))
LST <- c(1, month3.55)
DESeq.data.2 <- data.frame(DESeq.data.2)[, LST]

dds.2 <- DESeqDataSetFromMatrix(countData=DESeq.data.2[,-31],   # remove B35: poor enrichment 
                              colData=metaData.2[-30, ], 
                              design=~Pheno, tidy = TRUE)

dds.2 <- DESeq(dds.2)

resLFC.2 <- lfcShrink(dds.2, coef="Pheno_Y_vs_N", type="apeglm") # adaptive t prior shrinkage estimator
head(resLFC.2)
summary(resLFC.2)

write.csv(resLFC.2,"Additional Files/Time 2 DESeq2 Results.csv")

res.temp.2 <- resLFC.2[order(resLFC.2$padj),]
head(res.temp.2)

vsdata.2 <- vst(dds.2, blind=FALSE)
a.2=plotPCA(vsdata.2, intgroup="Pheno")
a.2
a.2$data
hist(resLFC.2$pvalue)

normal_RC.2 <- as.data.frame(counts(dds.2, normalized = T))
write.csv(normal_RC.2, "Additional Files/Time 2 Normalized.csv")




## Time 3 (C/D) ----

DESeq.data.3 <- data
DESeq.data.3 <- DESeq.data.3[, c(1,117:146)]

Dead <- c( rep("C", 7), "D", rep("C", 5), rep("D", 14), rep("C", 3))
metaData.3 <- data.frame(Sample = colnames(DESeq.data.3)[-c(1, 18)] , Pheno = as.factor(Dead[-17]))

dds.3 <- DESeqDataSetFromMatrix(countData=DESeq.data.3[, -18], # remove D14
                              colData=metaData.3, 
                              design=~Pheno, tidy = TRUE)

dds.3 <- DESeq(dds.3)

resLFC.3 <- lfcShrink(dds.3, coef="Pheno_D_vs_C", type="apeglm") # adaptive t prior shrinkage estimator
head(resLFC.3)
summary(resLFC.3)

write.csv(resLFC.3,"Additional Files/Time 3 DESeq Results.csv")

res.temp.3 <- resLFC.3[order(resLFC.3$padj),]
head(res.temp.3)

vsdata.3 <- vst(dds.3, blind=FALSE)
a.3=plotPCA(vsdata.3, intgroup="Pheno")
a.3
a.3$data
hist(resLFC.3$pvalue)

plot(a.3$data[,1],a.3$data[,2])
identify(a.3$data[,1],a.3$data[,2], a.3$data[,5], cex=0.8)

normal_RC.3 <- as.data.frame(counts(dds.3, normalized = T))
write.csv(normal_RC.3, "Additional Files/Time 3 Normalized.csv")




############################
## Time Specific Analysis ##
############################

## Time 1 (A) recounts ----

levels(as.factor(metaData.1$Pheno))
# *"N" is baseline

recountsA <- colSums(DESeq.data.1[,2:56])

a.1$data$recounts <- recountsA
View(a.1$data)

plot(a.1$data$PC1, a.1$data$recounts)
#plot(a.1$data$recounts, a.1$data$PC2)
plot(a.1$data$recounts, a.1$data$PC2, ylim=c(-6,6))

# PCA Plot
# - Bottom-right "Y":
#   A16, A26, A30, A41, A04 (do not have low recounts)

## Time 1 vs Time 3 ----

# Identify genes with significant adj. p-value from Time 1
lfc.sig.1 <- resLFC.1[!is.na(resLFC.1$padj),]
lfc.sig.1 <- lfc.sig.1[lfc.sig.1$padj < 0.1,]

# Split into LFC up and LFC down genes
lfc.sig.1$up <- lfc.sig.1$log2FoldChange > 0
up.gene.1 <- names(which(lfc.sig.1$up))
down.gene.1 <- names(which(!lfc.sig.1$up))

# Repeat for Time 3
lfc.sig.3 <- resLFC.3[!is.na(resLFC.3$padj),]
lfc.sig.3 <- lfc.sig.3[lfc.sig.3$padj < 0.1,]
lfc.sig.3$up <- lfc.sig.3$log2FoldChange > 0
up.gene.3 <- names(which(lfc.sig.3$up))
down.gene.3 <- names(which(!lfc.sig.3$up))

## Venn Diagram
summary(resLFC.1)
summary(resLFC.3)

library(VennDiagram)
venn.diagram(list(up.gene.1, up.gene.3),
             category.names = c("Up T1", "Up T3"),
             filename='Additional Files/up_venn.png',
             cat.default.pos='text')

venn.diagram(list(down.gene.1, down.gene.3),
             category.names = c("Down T1", "Down T3"),
             filename='Additional Files/down_venn.png',
             cat.default.pos='text',
             inverted=T)

venn.diagram(list(c(up.gene.1,down.gene.1), c(up.gene.3,down.gene.3)),
             category.names = c("T1", "T3"),
             filename='Additional Files/t1vt3.png',
             cat.default.pos='text')

## Chi-sq Test
t1 <- resLFC.1$padj < 0.1
t3 <- resLFC.3$padj < 0.1

table(t1, t3)
chisq.test(table(t1, t3))

## p-value distributions

#hist(resLFC.1$padj)

# Genes at Time 3 that were significant at Time 1
hist(resLFC.3$pvalue[resLFC.1$padj<0.1], nclass=20)

#hist(resLFC.3$padj)

# Genes at Time 1 that were significant at Time 3
hist(resLFC.1$pvalue[resLFC.3$padj<0.1], nclass=20)


## Time A vs C LFC plot
AvsC <- data.frame('lfcA'=resLFC.1$log2FoldChange,
                   'lfcC'=resLFC.3$log2FoldChange,
                   'filter'=(resLFC.1$padj<0.1 | resLFC.3$padj<0.1),
                   'significant'=(resLFC.1$padj<0.1 & resLFC.3$padj<0.1)) %>%
  drop_na() %>% filter(filter==T)

library(grid)
ggplot(data=AvsC) + geom_point(aes(x=lfcC, y=lfcA, color=significant, size=significant)) +
  scale_size_manual(values=c(0.8,1)) +
  labs(x='LFC at Time 1', y='LFC at Time 3', color='Significant\nin T1 & T3', size='Significant\nin T1 & T3') +
  geom_segment(aes(x=-1.5,y=0,xend=1.5,yend=0), arrow=arrow(length = unit(0.3, "cm"))) +
  geom_segment(aes(x=0,y=-1.1,xend=0,yend=1.1), arrow=arrow(length = unit(0.3, "cm"))) +
  xlim(-1.5, 1.5) + ylim(-1.1, 1.1) +
  theme_bw()


#####################
## Subset Analysis ##
#####################

## Set up ----

##
## Subsets are based on dendrogram/heatmap analysis
## Subset analysis is comparing notable subsets of patients with the other patients
##

# Split into subsets
metaData.sub.1 <- metaData.1
rownames(metaData.sub.1) <- NULL
metaData.sub.1$pr1 <- FALSE
metaData.sub.1$pr1[c(4,12,13,23,27,28,36,43)] <- TRUE #assign first subgroup of progressing subjects
metaData.sub.1$pr2 <- FALSE
metaData.sub.1$pr2[c(1,5,11,16,17,22,30,32,34,38,44,50)] <- TRUE #assign second subgroup


## Pr1 vs all non-progressing ----

# Pr1 is the most distinct group of progressing patients
sub.data.1.1 <- DESeq.data.1[,c(T,(metaData.sub.1$pr1==T | metaData.sub.1$Pheno=='N') & metaData.sub.1$pr2==F)]
sub.meta.1.1 <- metaData.sub.1[(metaData.sub.1$pr1==T | metaData.sub.1$Pheno=='N') & metaData.sub.1$pr2==F,c(1,3)]

dds.sub.1.1 <- DESeqDataSetFromMatrix(sub.data.1.1, 
                                      colData=sub.meta.1.1, 
                                      design=~pr1, tidy = TRUE)

dds.sub.1.1 <- DESeq(dds.sub.1.1)

resLFC.sub.1.1 <- lfcShrink(dds.sub.1.1, coef="pr1TRUE", type="apeglm") # adaptive t prior shrinkage estimator
head(resLFC.sub.1.1)
summary(resLFC.sub.1.1)

write.csv(resLFC.sub.1.1, "Additional Files/Pr1 vs nonprogression DESeq.csv")


# Pr2 vs all non-progressing ----

# Pr2 is the other subset of progressing patients
sub.data.1.2 <- DESeq.data.1[,c(T,(metaData.sub.1$pr2==T | metaData.sub.1$Pheno=='N') & metaData.sub.1$pr1==F)]
sub.meta.1.2 <- metaData.sub.1[(metaData.sub.1$pr2==T | metaData.sub.1$Pheno=='N') & metaData.sub.1$pr1==F,c(1,4)]

dds.sub.1.2 <- DESeqDataSetFromMatrix(sub.data.1.2, 
                                      colData=sub.meta.1.2, 
                                      design=~pr2, tidy = TRUE)

dds.sub.1.2 <- DESeq(dds.sub.1.2)

resLFC.sub.1.2 <- lfcShrink(dds.sub.1.2, coef="pr2TRUE", type="apeglm") # adaptive t prior shrinkage estimator
head(resLFC.sub.1.2)
summary(resLFC.sub.1.2)

write.csv(resLFC.sub.1.2, "Additional Files/Pr2 vs nonprogression DESeq.csv")


## Venn Diagrams ----

# Significant genes in Subset 1
a <- rownames(resLFC.sub.1.1[resLFC.sub.1.1$padj<0.1 & !is.na(resLFC.sub.1.1$padj),])
# Significant genes in Subset 2
b <- rownames(resLFC.sub.1.2[resLFC.sub.1.2$padj<0.1 & !is.na(resLFC.sub.1.2$padj),])
# All significant genes from Time 1
c <- rownames(lfc.sig.1)

venn.diagram(list(a,c),
             category.names = c("Subset1", "Original"),
             filename='Additional Files/subset1vsoriginal.png',
             cat.default.pos='text')


venn.diagram(list(a,b),
             category.names = c("Subset 1", "Subset 2"),
             filename='Additional Files/subset1vssubset2.png',
             cat.default.pos='text')





##################
## Paired Tests ##
##################

## A vs B ----
DESeq.paired.1.p <- data
DESeq.paired.1.np <- data

# Select only samples with progressing/non-progressing data
paired.samps.1.p <- c("A01","A04","A05","A14","A16","A19","A20","A25","A31","A36","A38","A41","A51","A57",
                    "B01","B04","B05","B14","B16","B19","B20","B25","B31","B36","B38","B41","B51","B57")
paired.samps.1.np <- c("A07","A10","A11","A12","A17","A18","A35","A37","A39","A42","A44","A49",
                      "B07","B10","B11","B12","B17","B18","B35","B37","B39","B42","B44","B49")

DESeq.paired.1.p <- DESeq.paired.1.p %>% select(all_of(c("Geneid",paired.samps.1.p)))
DESeq.paired.1.np <- DESeq.paired.1.np %>% select(all_of(c("Geneid",paired.samps.1.np)))

# Create design matrix
meta.paired.1.p <- data.frame("ID"=paired.samps.1.p,"Patient"=as.factor(c(1:14,1:14)), "Time"=as.factor(c(rep("A",14),rep("B",14))))
meta.paired.1.np <- data.frame("ID"=paired.samps.1.np,"Patient"=as.factor(c(1:12,1:12)), "Time"=as.factor(c(rep("A",12),rep("B",12))))

dds.paired.1.p <- DESeqDataSetFromMatrix(countData=DESeq.paired.1.p, 
                                         colData=meta.paired.1.p, 
                                         design=~ Patient + Time, tidy = TRUE)

dds.paired.1.np <- DESeqDataSetFromMatrix(countData=DESeq.paired.1.np, 
                                         colData=meta.paired.1.np, 
                                         design=~ Patient + Time, tidy = TRUE)

# Run DESeq
dds.paired.1.p <- DESeq(dds.paired.1.p)
dds.paired.1.np <- DESeq(dds.paired.1.np)

lfc.paired.1.p <- lfcShrink(dds.paired.1.p, coef="Time_B_vs_A", type="apeglm")
lfc.paired.1.np <- lfcShrink(dds.paired.1.np, coef="Time_B_vs_A", type="apeglm")

head(lfc.paired.1.p)
head(lfc.paired.1.np)

summary(lfc.paired.1.p)
summary(lfc.paired.1.np)

hist(lfc.paired.1.p$pvalue)
hist(lfc.paired.1.np$pvalue)

write.csv(lfc.paired.1.p,"Additional Files/AvsB Progression DESeq.csv")
write.csv(lfc.paired.1.np,"Additional Files/AvsB Nonprogression DESeq.csv")

# Normalize data
normal_RC.paired.1.p <- as.data.frame(counts(dds.paired.1.p, normalized = T))
normal_RC.paired.1.np <- as.data.frame(counts(dds.paired.1.np, normalized = T))

# Look at pairwise differences for specific genes
diff <- normal_RC.paired.1.p[,1:14] - normal_RC.paired.1.p[,15:28]
diff[9615,]




## B vs C/D ----
DESeq.paired.2.p <- data
DESeq.paired.2.np <- data

# Select only samples with progressing/non-progressing data
paired.samps.2.p <- c("B01","B04","B05","B14","B16","B19","B20","B25","B31","B36","B38","B41","B51","B57",
                    "D01","D04","D05","D14","D16","D19","D20","D25","D31","C36","D38","D41","D51","D57")
paired.samps.2.np <- c("B07","B10","B11","B12","B17","B18","B35","B37","B39","B42","B44","B49",
                     "C07","C10","C11","C12","C17","C18","C35","C37","C39","C42","C44","C49")

DESeq.paired.2.p <- DESeq.paired.2.p %>% select(all_of(c("Geneid",paired.samps.2.p)))
DESeq.paired.2.np <- DESeq.paired.2.np %>% select(all_of(c("Geneid",paired.samps.2.np)))

# Create design matrix
meta.paired.2.p <- data.frame("ID"=paired.samps.2.p,"Patient"=as.factor(c(1:14,1:14)), "Time"=as.factor(c(rep("B",14),rep("D",14))))
meta.paired.2.np <- data.frame("ID"=paired.samps.2.np,"Patient"=as.factor(c(1:12,1:12)), "Time"=as.factor(c(rep("B",12),rep("C",12))))

dds.paired.2.p <- DESeqDataSetFromMatrix(countData=DESeq.paired.2.p, 
                                         colData=meta.paired.2.p, 
                                         design=~ Patient + Time, tidy = TRUE)

dds.paired.2.np <- DESeqDataSetFromMatrix(countData=DESeq.paired.2.np, 
                                          colData=meta.paired.2.np, 
                                          design=~ Patient + Time, tidy = TRUE)

# Run DESeq
dds.paired.2.p <- DESeq(dds.paired.2.p)
dds.paired.2.np <- DESeq(dds.paired.2.np)

resultsNames(dds.paired.2.np)

lfc.paired.2.p <- lfcShrink(dds.paired.2.p, coef="Time_D_vs_B", type="apeglm")
lfc.paired.2.np <- lfcShrink(dds.paired.2.np, coef="Time_C_vs_B", type="apeglm")

head(lfc.paired.2.p)
head(lfc.paired.2.np)

summary(lfc.paired.2.p)
summary(lfc.paired.2.np)

hist(lfc.paired.2.p$pvalue)
hist(lfc.paired.2.np$pvalue)

write.csv(lfc.paired.2.p,"Additional Files/BvsD Progression DESeq.csv")
write.csv(lfc.paired.2.np,"Additional Files/BvsC Nonprogression DESeq.csv")

# Normalize data
normal_RC.paired.2.p <- as.data.frame(counts(dds.paired.2.p, normalized = T))
normal_RC.paired.2.np <- as.data.frame(counts(dds.paired.2.np, normalized = T))

# Look at pairwise differences for specific genes
diff <- normal_RC.paired.2.p[,1:14] - normal_RC.paired.2.p[,15:28]
diff[9615,]


## Comparing different LFC results ----

# BvD vs AvB
plot(lfc.paired.2.p$log2FoldChange, lfc.paired.1.p$log2FoldChange,
     xlab="BvsD LFC", ylab="AvsB LFC")

# PvNP vs AvB
plot(resLFC.1$log2FoldChange, lfc.paired.1.p$log2FoldChange, pch=16, cex=0.6,
     xlab="Original LFC", ylab="AvsB LFC")
plot(resLFC.1$log2FoldChange, lfc.paired.1.np$log2FoldChange, pch=16, cex=0.6,
     xlab="Original LFC", ylab="AvsB LFC")





## Separating Subsets (Time AvB) ----

DESeq.paired.1.pa <- data
DESeq.paired.1.pb <- data


# Select only samples with progressing/non-progressing data (pulled from Subset Analysis Section)
paired.samps.1.pa <- c("A04", "A16", "A31", "A41", # No time 3 data: "A26", "A30", "A50"
                       "B04", "B16", "B31", "B41")
paired.samps.1.pb <- c("A01", "A05", "A19", "A20", "A25", "A36", "A38", "A51", "A57", # No time 3 data: "A14", "A33", "A43"
                       "B01", "B05", "B19", "B20", "B25", "B36", "B38", "B51", "B57")


DESeq.paired.1.pa <- DESeq.paired.1.pa %>% select(all_of(c("Geneid",paired.samps.1.pa)))
DESeq.paired.1.pb <- DESeq.paired.1.pb %>% select(all_of(c("Geneid",paired.samps.1.pb)))

# Create design matrix
meta.paired.1.pa <- data.frame("ID"=paired.samps.1.pa,"Patient"=as.factor(c(1:4,1:4)), "Time"=as.factor(c(rep("A",4),rep("B",4))))
meta.paired.1.pb <- data.frame("ID"=paired.samps.1.pb,"Patient"=as.factor(c(1:9,1:9)), "Time"=as.factor(c(rep("A",9),rep("B",9))))

dds.paired.1.pa <- DESeqDataSetFromMatrix(countData=DESeq.paired.1.pa, 
                                         colData=meta.paired.1.pa, 
                                         design=~ Patient + Time, tidy = TRUE)

dds.paired.1.pb <- DESeqDataSetFromMatrix(countData=DESeq.paired.1.pb, 
                                         colData=meta.paired.1.pb, 
                                         design=~ Patient + Time, tidy = TRUE)


# Run DESeq
dds.paired.1.pa <- DESeq(dds.paired.1.pa)
dds.paired.1.pb <- DESeq(dds.paired.1.pb)

lfc.paired.1.pa <- lfcShrink(dds.paired.1.pa, coef="Time_B_vs_A", type="apeglm")
lfc.paired.1.pb <- lfcShrink(dds.paired.1.pb, coef="Time_B_vs_A", type="apeglm")

summary(lfc.paired.1.pa)
summary(lfc.paired.1.pb)

hist(lfc.paired.1.pa$pvalue)
hist(lfc.paired.1.pb$pvalue)

write.csv(lfc.paired.1.pa,"Additional Files/AvsB Progression A DESeq.csv")
write.csv(lfc.paired.1.pb,"Additional Files/AvsB Progression B DESeq.csv")

# Normalize data
normal_RC.paired.1.pa <- as.data.frame(counts(dds.paired.1.pa, normalized = T))
normal_RC.paired.1.pb <- as.data.frame(counts(dds.paired.1.pb, normalized = T))

# Venn Diagrams

# Split into LFC up and LFC down genes
lfc.paired.1.pa$up <- lfc.paired.1.pa$log2FoldChange > 0
lfc.paired.1.pb$up <- lfc.paired.1.pb$log2FoldChange > 0


# Up
a <- data$Geneid[lfc.paired.1.pa$padj<0.1 & !is.na(lfc.paired.1.pa$padj) & lfc.paired.1.pa$up]
b <- data$Geneid[lfc.paired.1.pb$padj<0.1 & !is.na(lfc.paired.1.pb$padj) & lfc.paired.1.pb$up]

venn.diagram(list(a,b),
             category.names = c("Subset 1", "Subset 2"),
             filename='Additional Files/subset1vssubset2 up.png',
             cat.default.pos='text')

# Down
a <- data$Geneid[lfc.paired.1.pa$padj<0.1 & !is.na(lfc.paired.1.pa$padj) & !lfc.paired.1.pa$up]
b <- data$Geneid[lfc.paired.1.pb$padj<0.1 & !is.na(lfc.paired.1.pb$padj) & !lfc.paired.1.pb$up]

venn.diagram(list(a,b),
             category.names = c("Subset 1", "Subset 2"),
             filename='Additional Files/subset1vssubset2 down.png',
             cat.default.pos='text')



## With samples that don't have time 3 data
DESeq.paired.1b.pa <- data
DESeq.paired.1b.pb <- data

paired.samps.1b.pa <- c("A04", "A16", "A31", "A41", "A26", "A30", "A50", # No time 3 data: "A26", "A30", "A50"
                       "B04", "B16", "B31", "B41", "B26", "B30", "B50")
paired.samps.1b.pb <- c("A01", "A05", "A19", "A20", "A25", "A36", "A38", "A51", "A57", "A14", "A33", "A43", # No time 3 data: "A14", "A33", "A43"
                       "B01", "B05", "B19", "B20", "B25", "B36", "B38", "B51", "B57", "B14", "B33", "B43")

DESeq.paired.1b.pa <- DESeq.paired.1b.pa %>% select(all_of(c("Geneid",paired.samps.1b.pa)))
DESeq.paired.1b.pb <- DESeq.paired.1b.pb %>% select(all_of(c("Geneid",paired.samps.1b.pb)))

meta.paired.1b.pa <- data.frame("ID"=paired.samps.1b.pa,"Patient"=as.factor(c(1:7,1:7)), "Time"=as.factor(c(rep("A",7),rep("B",7))))
meta.paired.1b.pb <- data.frame("ID"=paired.samps.1b.pb,"Patient"=as.factor(c(1:12,1:12)), "Time"=as.factor(c(rep("A",12),rep("B",12))))

dds.paired.1b.pa <- DESeqDataSetFromMatrix(countData=DESeq.paired.1b.pa, 
                                         colData=meta.paired.1b.pa, 
                                         design=~ Patient + Time, tidy = TRUE)

dds.paired.1b.pb <- DESeqDataSetFromMatrix(countData=DESeq.paired.1b.pb, 
                                         colData=meta.paired.1b.pb, 
                                         design=~ Patient + Time, tidy = TRUE)

dds.paired.1b.pa <- DESeq(dds.paired.1b.pa)
dds.paired.1b.pb <- DESeq(dds.paired.1b.pb)

lfc.paired.1b.pa <- lfcShrink(dds.paired.1b.pa, coef="Time_B_vs_A", type="apeglm")
lfc.paired.1b.pb <- lfcShrink(dds.paired.1b.pb, coef="Time_B_vs_A", type="apeglm")

summary(lfc.paired.1b.pa)
summary(lfc.paired.1b.pb)

write.csv(lfc.paired.1b.pa,"Additional Files/AvsB2 Progression A DESeq.csv")
write.csv(lfc.paired.1b.pb,"Additional Files/AvsB2 Progression B DESeq.csv")

# Venn Diagrams
a <- data$Geneid[lfc.paired.1b.pa$padj<0.1 & !is.na(lfc.paired.1b.pa$padj)]
b <- data$Geneid[lfc.paired.1b.pb$padj<0.1 & !is.na(lfc.paired.1b.pb$padj)]

venn.diagram(list(a,b),
             category.names = c("Subset 1", "Subset 2"),
             filename='Additional Files/subset1vssubset2 v2.png',
             cat.default.pos='text')


## Separating Subsets (Time BvC/D) ----

DESeq.paired.2.pa <- data
DESeq.paired.2.pb <- data


# Select only samples with progressing/non-progressing data (pulled from Subset Analysis Section)
paired.samps.2.pa <- c("B04", "B16", "B38", "B41", "B57",
                       "D04", "D16", "D38", "D41", "D57")
paired.samps.2.pb <- c("B01", "B05", "B19", "B20", "B25", "B31", "B36", "B51",
                       "D01", "D05", "D19", "D20", "D25", "D31", "C36", "D51")

DESeq.paired.2.pa <- DESeq.paired.2.pa %>% select(all_of(c("Geneid",paired.samps.2.pa)))
DESeq.paired.2.pb <- DESeq.paired.2.pb %>% select(all_of(c("Geneid",paired.samps.2.pb)))

# Create design matrix
meta.paired.2.pa <- data.frame("ID"=paired.samps.2.pa,"Patient"=as.factor(c(1:5,1:5)), "Time"=as.factor(c(rep("B",5),rep("D",5))))
meta.paired.2.pb <- data.frame("ID"=paired.samps.2.pb,"Patient"=as.factor(c(1:8,1:8)), "Time"=as.factor(c(rep("B",8),rep("D",8))))

dds.paired.2.pa <- DESeqDataSetFromMatrix(countData=DESeq.paired.2.pa, 
                                         colData=meta.paired.2.pa, 
                                         design=~ Patient + Time, tidy = TRUE)

dds.paired.2.pb <- DESeqDataSetFromMatrix(countData=DESeq.paired.2.pb, 
                                         colData=meta.paired.2.pb, 
                                         design=~ Patient + Time, tidy = TRUE)


# Run DESeq
dds.paired.2.pa <- DESeq(dds.paired.2.pa)
dds.paired.2.pb <- DESeq(dds.paired.2.pb)

lfc.paired.2.pa <- lfcShrink(dds.paired.2.pa, coef="Time_D_vs_B", type="apeglm")
lfc.paired.2.pb <- lfcShrink(dds.paired.2.pb, coef="Time_D_vs_B", type="apeglm")

head(lfc.paired.2.pa)
head(lfc.paired.2.pb)

summary(lfc.paired.2.pa)
summary(lfc.paired.2.pb)

hist(lfc.paired.2.pa$pvalue)
hist(lfc.paired.2.pb$pvalue)

write.csv(lfc.paired.2.pa,"Additional Files/BvsD Progression A DESeq.csv")
write.csv(lfc.paired.2.pb,"Additional Files/BvsD Progression B DESeq.csv")

# Normalize data
normal_RC.paired.2.pa <- as.data.frame(counts(dds.paired.2.pa, normalized = T))
normal_RC.paired.2.pb <- as.data.frame(counts(dds.paired.2.pb, normalized = T))




###################################
## Misc. Cleaning and Outputting ##
###################################

## Prep for heatmap/dendogram ----

# All time 1 samples
sig.genes.1 <- cbind(data$Geneid, normal_RC.1)
sig.genes.1 <- sig.genes.1[resLFC.1$padj<0.1,] %>% drop_na()
write_csv(sig.genes.1, 'Additional Files/time1_heatmap_data.csv')

# Just time 1 progression samples
sig.genes.p.1 <- sig.genes.1 %>% select(c(`data$Geneid`,metaData.1$ID[metaData.1$Pheno=='Y']))
write_csv(sig.genes.p.1, 'Additional Files/time1_heatmap_data_p.csv')

# All time 3 samples
sig.genes.3 <- cbind(data$Geneid, normal_RC.3)
sig.genes.3 <- sig.genes.3[resLFC.3$padj<0.1,] %>% drop_na()
write_csv(sig.genes.3, 'Additional Files/time3_heatmap_data.csv')

# Just time 3 progression samples
sig.genes.p.3 <- sig.genes.3 %>% select(c(`data$Geneid`,metaData.3$Sample[metaData.3$Pheno=='D']))
write_csv(sig.genes.p.3, 'Additional Files/time3_heatmap_data_p.csv')

## Clean clinical table ----
samps <- c(colnames(DESeq.data.1)[2:56], colnames(DESeq.data.2)[2:56], colnames(DESeq.data.3)[2:31])
Final.2 <- Final %>% filter(ID %in% samps)

write.csv(Final.2, "Additional Files/clinical table.csv")

## Update data table to drop out bad samples
data_new <- data %>% select(-c("A08", "A09", "A13", "A34", "A40", "A46", "A47", "A58", "B08", "B09", "B34", "B40", "B47", "B58", "C09", "D13", "D34"))

## Filter out genes with a max recount less than 8 ----
max <- apply(data_new[,7:162], 1, max)
sum(max<8)
hist(max)

data_clean <- data_new[max>8,]

## Write out final data ----
fwrite(data_clean, "Qianxia_raw_counts_by_Genes_FINAL.csv")
nrow(data_clean)

## Write out normalized data for Time 1 ----
fwrite(normal_RC.1[max>5,], "Additional Files/GSEA/Time 1 normalized data.csv")
