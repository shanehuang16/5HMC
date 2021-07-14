#######################################
## Time Comparison Plotting Function ##
#######################################

library(tidyverse)
library(data.table)

data.normal <- fread("Normalized Data All.csv")

a.p.all <- c("A01","A04","A05","A14","A15","A16","A19","A20","A25","A26","A30","A31","A33","A36","A38","A41","A43","A50","A51","A57")
b.p.all <- c("B01","B04","B05","B14","B16","B19","B20","B25","B26","B30","B31","B33","B36","B38","B41","B43","B46","B50","B51","B57")
c.p.all <- c("D01","D04","D05","D14","D15","D16","D19","D20","D25","D31","C36","D38","D41","D51","D57")
a.np.all <- c("A02","A03","A06","A07","A10","A11","A12","A17","A18","A21","A22","A23","A24","A27","A28","A29","A32","A35","A37","A39","A42","A44","A45","A48","A49","A52","A53","A54","A55","A56","A59","A60","A61","A62","A63")
b.np.all <- c("B02","B03","B06","B07","B10","B11","B12","B17","B18","B21","B22","B23","B24","B27","B28","B29","B32","B35","B37","B39","B42","B44","B45","B48","B49","B52","B53","B54","B55","B56","B59","B60","B61","B62","B63")
c.np.all <- c("E02","E03","E06","C07","C10","C11","C12","C17","C18","C35","C37","C39","C42","C44","C49")

a.p.comp <- c("A01","A04","A05","A14","A16","A19","A20","A25","A31","A36","A38","A41","A51","A57")
b.p.comp <- c("B01","B04","B05","B14","B16","B19","B20","B25","B31","B36","B38","B41","B51","B57")
c.p.comp <- c("D01","D04","D05","D14","D16","D19","D20","D25","D31","C36","D38","D41","D51","D57")
a.np.comp <- c("A02","A03","A06","A07","A10","A11","A12","A17","A18","A35","A37","A39","A42","A44","A49")
b.np.comp <- c("B02","B03","B06","B07","B10","B11","B12","B17","B18","B35","B37","B39","B42","B44","B49")
c.np.comp <- c("E02","E03","E06","C07","C10","C11","C12","C17","C18","C35","C37","C39","C42","C44","C49")

a.p.sub1.all <- c("A04","A15","A16","A26","A30","A31","A41","A50")
b.p.sub1.all <- c("B04","B16","B26","B30","B31","B41","B50")
c.p.sub1.all <- c("D04","D16","D38","D41","D57")
a.p.sub2.all <- c("A01","A05","A14","A19","A20","A25","A33","A36","A38","A43","A51","A57")
b.p.sub2.all <- c("B01","B05","B14","B19","B20","B25","B33","B36","B38","B43","B51","B57")
c.p.sub2.all <- c("D01","D05","D15","D19","D20","D25","D31","C36","D51")

a.p.sub1.comp <- c("A04","A16","A31","A41")
b.p.sub1.comp <- c("B04","B16","B31","B41")
c.p.sub1.comp <- c("D04","D16","D38","D41","D57")
a.p.sub2.comp <- c("A01","A05","A19","A20","A25","A36","A38","A51","A57")
b.p.sub2.comp <- c("B01","B05","B19","B20","B25","B36","B38","B51","B57")
c.p.sub2.comp <- c("D01","D05","D19","D20","D25","D31","C36","D51")

sets <- list(list(a.p.all,b.p.all,c.p.all,a.np.all,b.np.all,c.np.all),
             list(a.p.comp,b.p.comp,c.p.comp,a.np.comp,b.np.comp,c.np.comp),
             list(a.p.sub1.all,b.p.sub1.all,c.p.sub1.all,a.np.all,b.np.all,c.np.all),
             list(a.p.sub2.all,b.p.sub2.all,c.p.sub2.all,a.np.all,b.np.all,c.np.all),
             list(a.p.sub1.comp,b.p.sub1.comp,c.p.sub1.comp,a.np.comp,b.np.comp,c.np.comp),
             list(a.p.sub2.comp,b.p.sub2.comp,c.p.sub2.comp,a.np.comp,b.np.comp,c.np.comp),
             # Custom scenarios:
             list(a.p.sub1.all, b.p.all, c.p.sub2.all, a.np.all,b.np.all,c.np.all),
             list(a.p.sub2.all, b.p.all, c.p.sub1.all, a.np.all,b.np.all,c.np.all)
             )




compplot <- function(gene, scenario){
  # Constant settings
  mean <- rep(NA,6)
  sd.u <- rep(NA,6)
  sd.l <- rep(NA,6)
  time <- as.factor(c('A','B','C','A','B','C'))
  pheno <- as.factor(c(rep('Progression',3), rep('Non-Progression',3)))
  
  # Select scenario set
  set <- sets[[scenario]]
  
  # Select specfic gene data
  dat <- data.normal %>% filter(Geneid == gene)
  
  # Calculate means and sds
  for (i in 1:6){
    select <- as.numeric(as.vector(dat %>% select(all_of(set[[i]]))))
    mean[i] <- mean(select)
    sd.u[i] <- mean(select) + sd(select)
    sd.l[i] <- mean(select) - sd(select)
  }
  
  # Create dataset
  df <- data.frame("mean"=mean, "sd.u"=sd.u, "sd.l"=sd.l, "time"=time, "pheno"=pheno)
  
  # Create plot
  ggplot(data=df, aes(x=time, y=mean, color=pheno, group=pheno)) +
    geom_point() +
    geom_errorbar(aes(ymax=sd.u, ymin=sd.l), width=.25) +
    stat_summary(fun=sum, geom="line") +
    labs(title=gene, y="Mean Normalized Count", x="Time Point", color="Condition")
}



## Analysis ----
# ***Run DESeq Analysis Final first

# Finding genes to plot
gene.paired.1.pa <- data$Geneid[lfc.paired.1.pa$padj<0.1 & !is.na(lfc.paired.1.pa$padj)]
gene.sub.1.1 <- data$Geneid[lfc.sub.1.1$padj<0.1 & !is.na(lfc.sub.1.1$padj)]
gene.paired.2.pb <- data$Geneid[lfc.paired.2.pb$padj<0.1 & !is.na(lfc.paired.2.pb$padj)]
gene.sub.3.2 <- data$Geneid[lfc.sub.3.2$padj<0.1 & !is.na(lfc.sub.3.2$padj) & lfc.sub.3.2$log2FoldChange>0]

gene.plot <- Reduce(intersect, list(gene.paired.1.pa, gene.sub.1.1, gene.paired.2.pb, gene.sub.3.2))
gene.plot

# Plotting
compplot("ATP8B2", 7) # Sub1, All, Sub2



head(data$Geneid[order(lfc.1$padj)])
compplot("PDGFA", 7) # Sub1, All, Sub2
"PDGFA" %in% gene.plot

