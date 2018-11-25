###
#
#i.n. vs i.m. Fluidigm anaysis script
#Written by Terry Neeman and Cameron Jack, ANU, 2014
#Adapted by Irwan Jaeson, ANU, 2018
#
#This is a script to analyse data between i.n. FPV and Ad5, and i.m. FPV and Ad5
#
###

library(lme4)
library(lmerTest)
library(texreg) # nice table layout
library(ggplot2)
library(reshape2)
library(plyr)

data <- read.csv('IJ-Fluidigm-18-12-17-FPV-Ad5.csv')

print("ILC analysis")

# these are effectively "missing data". We compensate by using 
# non-parametric stats - also 40 is the highest number possible from
# the fluidigm machine anyway and corresponds to "no expression".

data[data==999.0] = 40.0
values <- data[ ,7:45]

kept_genes <- NULL
for (i in names(values)) {
  num40 <- sum(values[[i]]==40.0)
  totalnum <- nrow(values[i])
  if (((totalnum - num40)/totalnum) > 0.05) {
    if (is.null(kept_genes)) {
      kept_genes <- as.data.frame(values[i])
    } else {
      kept_genes <- cbind(kept_genes, values[i])
    }
  }
  num40 <- 0
  totalnum <- 0
}

values <- as.matrix(kept_genes)

# normalise scores by L32 control

values <- ifelse(values != 40, values - data[,4], values)

# calculate spearman rank correlation matrix
x <- cor(values, method='spearman')

pca1 = prcomp(x, scale.=TRUE, center=TRUE, tol=0)
summary(pca1)


#Calculate and plot Proportion of variance
pca_proportions <- pca1$sdev^2/sum(pca1$sdev^2)*100 

pca_proportions <- as.data.frame(pca_proportions[1:4])

colnames(pca_proportions) <- "PoV"

pca_proportions$PC <- c("PC1", "PC2", "PC3", "PC4") 

tiff("Proportion of Variation Routes.tiff ", units = "in", width = 10, height = 5, res = 600)

ggplot(pca_proportions, aes(pca_proportions$PC,pca_proportions$PoV)) +
  geom_bar(stat = 'identity', fill = "blue", colour = "black") +
  labs(x = "", y = "Proportion of Variation (%)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title.y = element_text(size = 20), axis.text.y = element_text(color = "black", size = 20), axis.text.x = element_text(color = "black", size = 20))
  
dev.off()

# Take first 4 PCs as this should cover most of the variance
y <- as.matrix(pca1$rotation[,1:4])
scores <- (40-values)%*%(y)
#summary(scores)

df <- as.data.frame(values)
df[, c("Vaccine", "Route", "Experiment")] <- data[, c("Vaccine", "Route", "Experiment")]

print("PC1") 
fm1 <- lmer(scores[,1] ~ Vaccine * Route + (1|Experiment), data=df)
screenreg(fm1, single.row = TRUE, naive = TRUE)
anova(fm1)


print("PC2")
fm2 <- lmer(scores[,2] ~ Vaccine * Route + (1|Experiment), data=df)
screenreg(fm2, single.row = TRUE, naive = TRUE)
anova(fm2)
#anova(lmer(scores[,2] ~ as.factor(Group) + (1|Experiment), data=df))

print("PC3")
fm3 <- lmer(scores[,3] ~ Vaccine * Route + (1|Experiment), data=df)
screenreg(fm3, single.row = TRUE, naive = TRUE)
anova(fm3)
#anova(lmer(scores[,3] ~ as.factor(Group) + (1|Experiment), data=df))

print("PC4")
fm4 <- lmer(scores[,4] ~ Vaccine * Route + (1|Experiment), data=df)
screenreg(fm4, single.row = TRUE, naive = TRUE)
anova(fm4)
#anova(lmer(scores[,3] ~ as.factor(Group) + (1|Experiment), data=df))

# Now make heat maps of the correlation coefficients

# to choose ordering, cluster it with affinity propagation
library(apcluster)
apclus <- apcluster(negDistMat(r=1.5), x)
heatmap(apclus)
apclus
order <- append(apclus[[4]], c(apclus[[6]], apclus[[2]], apclus[[3]], apclus[[1]], apclus[[5]]))

x1<-cor(df[,row.names(as.data.frame(order,optional=TRUE))], method='spearman')
colnames(x1)

names_x1<-colnames(x1)
qplot(x=Var1, y=Var2, data=melt(x1), fill=value, geom="tile")

x1plot<-melt(x1)

x1plot$Var1 <- factor(x1plot$Var1, levels=names_x1)
x1plot$Var2 <- factor(x1plot$Var2, levels=rev(names_x1))

tiff('ILC2 correlation heatmap route.tiff', units="in", width=15, height=12, res=600)

ggplot(data =  x1plot, aes(x = Var1, y = Var2)) +
  labs(x = "", y = "") +
  geom_tile(aes(fill = value), colour = "white") +
  theme(axis.text.x = element_text(size = rel(2), angle = 270, hjust = 0, colour = "grey50", face = "italic")) +
  theme(axis.text.y = element_text(size = rel(2), vjust = 0, colour = "grey50", face = "italic")) +
  geom_text(aes(label = sprintf("%0.2f",value)), cex=6,vjust = 0.5) +
  scale_fill_gradient2(low="darkblue", mid="white", high="red", limits = c(-1,1)) +
  theme(legend.text = element_text(size = 8), legend.text.align = 0.5) +
  theme(legend.justification = "top", legend.text = element_text(size = rel(1.5))) + 
  theme(legend.title = element_text(colour = "white"))
 

dev.off()

#Plot PCA vectors

loadings <- as.data.frame(pca1$rotation[, 1:4])
loadings$var <- colnames(x)
loadings <- loadings[order(loadings$PC1),]

library(tidyr)

load <- melt(loadings)
ord <- rep(18:1, 4)

tiff("Route loadings.tiff", units = "in", width = 8, height = 5, res = 600)

ggplot(data = melt(loadings)) +
  geom_bar(aes(x=reorder(var, ord), y=value), stat="identity", fill=3)+
  theme(axis.text.x = element_text(size = rel(0.7), angle = 300, hjust = 0, colour = "black", face = "italic"))+
  facet_wrap(~variable)+
  xlab("Gene")+ ylab("Value") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12,face="bold"))

dev.off()
