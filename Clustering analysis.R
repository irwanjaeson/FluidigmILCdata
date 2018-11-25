###
#
#Stat3, Stat6, Tgfb1 and Ifngr1 Fluidigm analysis script
#Written by Terry Neeman and Cameron Jack, ANU, 2014
#Adapted by Irwan Jaeson, ANU, 2018
#
###

library(lme4)
library(lmerTest)
library(texreg) # nice table layout
library(ggplot2)
library(reshape2)
library(tidyr)
library(plyr)

data <- read.csv('IJ-Fluidigm-10-07-18-all.csv')

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

values <- ifelse(values != 40, values - data[,4], values)

# calculate spearman rank correlation matrix
x <- cor(values, method='spearman')

pca1 <- prcomp(x, scale. = TRUE, center = TRUE)
summary(pca1)

y <- as.matrix(pca1$rotation[,1:4])

scores <- (40-values)%*%y
#summary(scores)

df <- as.data.frame(values)
df[, c("Vaccine", "Route", "Experiment")] <- data[, c("Vaccine", "Route", "Experiment")]

#Plot Stat6, Tgfb1, Stat3, Ifngr1, Il7R, Il2ra, Gata3, Rora, Il5, Il9r, Icos and Cd86

set1 <- values[, c("Stat6", "Tgfb1", "Stat3", "Ifngr1", "Il7r", "Il2ra", "Gata3", "Rora", "Il5", "Il9r", "Icos", "Cd86")]
cor_set1 <- cor(set1, method = 'spearman')
pca_set1 = prcomp(cor_set1, scale. = TRUE, center = TRUE, tol = 0)
summary(pca_set1)

#PC1 and PC2 should capture most of the variance
matrix_set1 <- as.matrix(pca_set1$rotation[,1:2])
scores_set1 <- (40-set1)%*%matrix_set1
scores_set1 <- jitter(scores_set1, amount = 3)
scores_set1_df <- as.data.frame(scores_set1)
scores_set1_df[, c("Vaccine", "Route", "Experiment")] <- data[, c("Vaccine", "Route", "Experiment")]

fm1 <- lmer(scores_set1[,1] ~ Vaccine * Route + (1|Experiment), data=df)
screenreg(fm1, single.row = TRUE, naive = TRUE)
anova(fm1)

df_set1 <- as.data.frame(set1)
df_set1[, c("Vaccine", "Route", "Experiment")] <- data[, c("Vaccine", "Route", "Experiment")]

#Determine the appropriate number of clusters for this set

library(NbClust)
library(factoextra)
fviz_nbclust(NbClust(scores_set1_df[,1:2], distance = "euclidean", min.nc = 2,
                     max.nc = 10, method = "kmeans"))

cl1 <- kmeans(scores_set1_df[,1:2], 2, nstart = 50)
scores_set1_df <- cbind(scores_set1_df, as.character(cl1$cluster))
scores_set1_df <- rename(scores_set1_df, c("as.character(cl1$cluster)" = "cluster"))

scores_set1_df <- cbind(scores_set1_df, as.vector(1:207))
scores_set1_df <- rename(scores_set1_df,c("as.vector(1:207)" = "ID"))

lmer(scores_set1[,1] ~ Vaccine * Route + (1|Experiment), data=df_set1) #PC1
anova(lmer(scores_set1[,1] ~ Vaccine * Route + (1|Experiment), data=df_set1))
lmer(scores_set1[,2] ~ Vaccine * Route + (1|Experiment), data=df_set1) #PC2
anova(lmer(scores_set1[,2] ~ Vaccine * Route + (1|Experiment), data=df_set1))

mean_scores_set1 <- aggregate(scores_set1_df[, 1:2], list(scores_set1_df$Vaccine, scores_set1_df$Route), mean)

names(mean_scores_set1)[names(mean_scores_set1) == "Group.1"] <- "Vaccine"
names(mean_scores_set1)[names(mean_scores_set1) == "Group.2"] <- "Route"

VRint_set1 <- interaction(scores_set1_df$Route, scores_set1_df$Vaccine)
VRint_mean_set1 <- interaction(mean_scores_set1$Route, mean_scores_set1$Vaccine)

tiff("clusterall_ILC2genes.tiff", width = 12, height = 8, units = "in", res = 600)

ggplot(data = scores_set1_df, aes(x = scores_set1_df$PC1, y = scores_set1_df$PC2, colour = VRint_set1, shape = VRint_set1)) + 
  #geom_point(size = 2.5) +
  #stat_ellipse(aes(group = scores_set1_df$cluster)) +
  theme_bw() + 
  geom_point(data = scores_set1_df, aes(x = scores_set1_df$PC1, y = scores_set1_df$PC2, colour = VRint_set1, shape = VRint_set1), size = 4, colour = "black") +
  geom_point(data = scores_set1_df, aes(x = scores_set1_df$PC1, y = scores_set1_df$PC2, colour = VRint_set1, shape = VRint_set1), size = 3) +
  #geom_text(aes(label = scores_set1_df$ID)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.margin = margin(1,1,1,1, "cm")) +
  geom_vline(xintercept = 41, size= 0.5) +
  geom_hline(yintercept = 25, size = 0.5) +
  labs(x = "PC1", y = "PC2") +
  theme(axis.title = element_text(size = 18, face = "bold")) +
  theme(axis.text = element_text(size = 16)) +
  scale_shape_manual(values = c(16,16,18,17,17), name = "Vaccine/Route", labels = (c("FPV/i.m.", "FPV/i.n.", "MVA/i.m.", "Ad5/i.m.", "Ad5/i.n."))) +
  scale_colour_manual("Vaccine/Route", values = c(2, 3, "yellow", 4, "grey"), labels = (c("FPV/i.m.", "FPV/i.n.", "MVA/i.m.", "Ad5/i.m.", "Ad5/i.n."))) +
  theme(legend.position = "bottom", legend.direction = "horizontal") +
  theme(legend.key.size = unit(1.5, "cm")) + 
  theme(legend.background = element_rect(colour = "black", fill= "white", size = 0.5), legend.key = element_rect(colour = "white", fill = "white")) +
  theme(legend.text = element_text(size = "16"), legend.title = element_text(size = "16"))

dev.off()

#plot TGFb1, STAT3, STAT6, IFNgRI

set5 <- values[, c("Tgfb1", "Stat3","Stat6", "Ifngr1")]
cor_set5 <- cor(set5, method = 'spearman')
pca_set5 = prcomp(cor_set5, scale. = TRUE, center = TRUE, tol = 0)
summary(pca_set5)

#PC1 and PC2 should capture most of the variance
matrix_set5 <- as.matrix(pca_set5$rotation[,1:2])
scores_set5 <- (40-set5)%*%matrix_set5
scores_set5 <- jitter(scores_set5, amount = 2)
df_set5 <- as.data.frame(set5)
df_set5[, c("Vaccine", "Route", "Experiment")] <- data[, c("Vaccine", "Route", "Experiment")]

scores_set5_df <- as.data.frame(scores_set5)
scores_set5_df[, c("Vaccine", "Route", "Experiment")] <- data[, c("Vaccine", "Route", "Experiment")]

cl5 <- kmeans(scores_set5_df[,1:2], 12, nstart = 1000)
scores_set5_df <- cbind(scores_set5_df, as.character(cl5$cluster))
scores_set5_df <- rename(scores_set5_df, c("as.character(cl5$cluster)" = "cluster"))

scores_set5_df <- cbind(scores_set5_df, as.vector(1:207))
scores_set5_df <- rename(scores_set5_df,c("as.vector(1:207)" = "ID"))

cl5_mean <- aggregate(scores_set5_df[,1:2], list(scores_set5_df$cluster), mean)

VRint_set5 <- interaction(scores_set5_df$Route, scores_set5_df$Vaccine)

tiff("clusterall.tiff", width = 9, height = 8, units = "in", res = 600)

ggplot(data = scores_set5_df, aes(x = scores_set5_df$PC1, y = scores_set5_df$PC2, colour = VRint_set5, shape = VRint_set5)) + 
  #geom_point(size = 2.5) +
  stat_ellipse(aes(group = scores_set5_df$cluster), type = "norm", show.legend = FALSE, na.rm = TRUE) +
  theme_bw() +
  geom_point(data = scores_set5_df, aes(x = scores_set5_df$PC1, y = scores_set5_df$PC2, colour = VRint_set5, shape = VRint_set5), size = 4, colour = "black") +
  geom_point(data = scores_set5_df, aes(x = scores_set5_df$PC1, y = scores_set5_df$PC2, colour = VRint_set5, shape = VRint_set5), size = 3) +
  #geom_text(aes(label = scores_set5_df$cluster), colour = "black", position = position_jitter(width = 2, height = 2)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.margin = margin(1,1,1,1, "cm")) +
  geom_vline(xintercept = 3, size= 0.5) +
  geom_hline(yintercept = -1, size = 0.5) +
  labs(x = "PC1", y = "PC2") +
  theme(axis.title = element_text(size = 18, face = "bold")) +
  theme(axis.text = element_text(size = 16)) +
  theme(axis.title.x = element_text(vjust = -1)) +
  scale_shape_manual(values = c(16,16,18,17,17), name = "Vaccine/Route", labels = (c("FPV/i.m.", "FPV/i.n.", "MVA/i.m.", "Ad5/i.m.", "Ad5/i.n."))) +
  scale_colour_manual("Vaccine/Route", values = c(2, 3, "yellow", 4, "grey"), labels = (c("FPV/i.m.", "FPV/i.n.", "MVA/i.m.", "Ad5/i.m.", "Ad5/i.n."))) +
  #geom_text(aes(label = scores_set5_df$cluster), color = "black", size = 10, check_overlap = TRUE) +
  theme(legend.position = "bottom", legend.direction = "horizontal") +
  theme(legend.key.size = unit(1.5, "cm")) + 
  theme(legend.background = element_rect(colour = "black", fill= "white", size = 1), legend.key = element_rect(colour = "white", fill = "white")) +
  theme(legend.text = element_text(size = "16"), legend.title = element_text(size = "16"))

dev.off()

cluster_all <- cbind(scores_set5_df, df$STAT3, df$STAT6, df$TGFb1, df$IFNgR1)

#Generate expression combinations for Stat3, Stat6, Tgfb1 and Ifngr1

b.data.inFPV <- data[data$Vaccine == "A" & data$Route == "N",c("Stat3", "Stat6", "Tgfb1", "Ifngr1")]
b.data.inFPV <- as.matrix(b.data.inFPV)
b.data.inFPV[b.data.inFPV != 40] <- 1
b.data.inFPV[b.data.inFPV == 40] <- 0
b.data.inFPV <- as.data.frame(b.data.inFPV)
b.data.inFPV <- b.data.inFPV[!(b.data.inFPV$STAT3 == 0 & b.data.inFPV$STAT6 == 0 & b.data.inFPV$TGFb1 == 0 & b.data.inFPV$IFNgR1 == 0),]

R <- NULL
b.data.inFPV$R <- NULL
for (i in row.names.data.frame(b.data.inFPV)) {
  R[i] <- as.vector(paste(b.data.inFPV[i,], collapse=""))
  
}
b.data.inFPV <- cbind(b.data.inFPV, R)

freq.inFPV <- count(b.data.inFPV, 'R')

prop <- NULL
freq.inFPV$prop <- NULL
for(i in row.names.data.frame(freq.inFPV)) {
  prop[i] <- as.numeric((freq.inFPV[i,2])/sum(freq.inFPV$freq)*100)
}

freq.inFPV <- cbind(freq.inFPV,prop)

b.data.imFPV <- data[data$Vaccine == "A" & data$Route == "M",c("STAT3", "STAT6", "TGFb1", "IFNgR1")]
b.data.imFPV <- as.matrix(b.data.imFPV)
b.data.imFPV[b.data.imFPV != 40] <- 1
b.data.imFPV[b.data.imFPV == 40] <- 0
b.data.imFPV <- as.data.frame(b.data.imFPV)
b.data.imFPV <- b.data.imFPV[!(b.data.imFPV$STAT3 == 0 & b.data.imFPV$STAT6 == 0 & b.data.imFPV$TGFb1 == 0 & b.data.imFPV$IFNgR1 == 0),]

R <- NULL
b.data.imFPV$R <- NULL
for (i in row.names.data.frame(b.data.imFPV)) {
  R[i] <- as.vector(paste(b.data.imFPV[i,], collapse=""))
  
}
b.data.imFPV <- cbind(b.data.imFPV, R)

freq.imFPV <- count(b.data.imFPV, 'R')

prop <- NULL
freq.imFPV$prop <- NULL
for(i in row.names.data.frame(freq.imFPV)) {
  prop[i] <- as.numeric((freq.imFPV[i,2])/sum(freq.imFPV$freq)*100)
}

freq.imFPV <- cbind(freq.imFPV,prop)

b.data.inAd5 <- data[data$Vaccine == "C" & data$Route == "N",c("STAT3", "STAT6", "TGFb1", "IFNgR1")]
b.data.inAd5 <- as.matrix(b.data.inAd5)
b.data.inAd5[b.data.inAd5 != 40] <- 1
b.data.inAd5[b.data.inAd5 == 40] <- 0
b.data.inAd5 <- as.data.frame(b.data.inAd5)
b.data.inAd5 <- b.data.inAd5[!(b.data.inAd5$STAT3 == 0 & b.data.inAd5$STAT6 == 0 & b.data.inAd5$TGFb1 == 0 & b.data.inAd5$IFNgR1 == 0),]

R <- NULL
b.data.inAd5$R <- NULL
for (i in row.names.data.frame(b.data.inAd5)) {
  R[i] <- as.vector(paste(b.data.inAd5[i,], collapse=""))
  
}
b.data.inAd5 <- cbind(b.data.inAd5, R)

freq.inAd5 <- count(b.data.inAd5, 'R')

prop <- NULL
freq.inAd5$prop <- NULL
for(i in row.names.data.frame(freq.inAd5)) {
  prop[i] <- as.numeric((freq.inAd5[i,2])/sum(freq.inAd5$freq)*100)
}

freq.inAd5 <- cbind(freq.inAd5,prop)

b.data.imAd5 <- data[data$Vaccine == "C" & data$Route == "M",c("STAT3", "STAT6", "TGFb1", "IFNgR1")]
b.data.imAd5 <- as.matrix(b.data.imAd5)
b.data.imAd5[b.data.imAd5 != 40] <- 1
b.data.imAd5[b.data.imAd5 == 40] <- 0
b.data.imAd5 <- as.data.frame(b.data.imAd5)
b.data.imAd5 <- b.data.imAd5[!(b.data.imAd5$STAT3 == 0 & b.data.imAd5$STAT6 == 0 & b.data.imAd5$TGFb1 == 0 & b.data.imAd5$IFNgR1 == 0),]

R <- NULL
b.data.imAd5$R <- NULL
for (i in row.names.data.frame(b.data.imAd5)) {
  R[i] <- as.vector(paste(b.data.imAd5[i,], collapse=""))
  
}
b.data.imAd5 <- cbind(b.data.imAd5, R)

freq.imAd5 <- count(b.data.imAd5, 'R')

prop <- NULL
freq.imAd5$prop <- NULL
for(i in row.names.data.frame(freq.imAd5)) {
  prop[i] <- as.numeric((freq.imAd5[i,2])/sum(freq.imAd5$freq)*100)
}

freq.imAd5 <- cbind(freq.imAd5,prop)

b.data.imMVA <- data[data$Vaccine == "B" & data$Route == "M",c("STAT3", "STAT6", "TGFb1", "IFNgR1")]
b.data.imMVA <- as.matrix(b.data.imMVA)
b.data.imMVA[b.data.imMVA != 40] <- 1
b.data.imMVA[b.data.imMVA == 40] <- 0
b.data.imMVA <- as.data.frame(b.data.imMVA)
b.data.imMVA <- b.data.imMVA[!(b.data.imMVA$STAT3 == 0 & b.data.imMVA$STAT6 == 0 & b.data.imMVA$TGFb1 == 0 & b.data.imMVA$IFNgR1 == 0),]

R <- NULL
b.data.imMVA$R <- NULL
for (i in row.names.data.frame(b.data.imMVA)) {
  R[i] <- as.vector(paste(b.data.imMVA[i,], collapse=""))
  
}
b.data.imMVA <- cbind(b.data.imMVA, R)

freq.imMVA <- count(b.data.imMVA, 'R')

prop <- NULL
freq.imMVA$prop <- NULL
for(i in row.names.data.frame(freq.imMVA)) {
  prop[i] <- as.numeric((freq.imMVA[i,2])/sum(freq.imMVA$freq)*100)
}

freq.imMVA <- cbind(freq.imMVA,prop)

binarydata <- df
binarydata[binarydata < 40] = 1
binarydata[binarydata == 40] = 0
binarydata$Experiment <- NULL

write.csv(cluster_all, file = "Cluster all.csv")

freq.inFPV
freq.inAd5
freq.imFPV
freq.imMVA
freq.imAd5
