##This is specifically to plot the data using ggplot2 for visualising differences between routes and vaccine vector

library(multcomp)
library(plyr)

scores.df <- as.data.frame(scores)

scores.df[, c("Vaccine", "Route", "Experiment")] <- data[, c("Vaccine", "Route", "Experiment")]

library(ggplot2)
 
mean_scores <- aggregate(scores.df[, 1:3], list(scores.df$Vaccine, scores.df$Route), mean)
mean_scores <- rename(mean_scores, c("Group.1" = "Vaccine", "Group.2" = "Route"))

#Group variables by interaction between Vaccine and Route

VRint <- interaction(scores.df$Route, scores.df$Vaccine)
VRint_mean <- interaction(mean_scores$Route, mean_scores$Vaccine)

p1a <- ggplot(data = scores.df, aes(x = scores.df$PC1, y = scores.df$PC2, colour = VRint, shape = VRint)) + 
  geom_point(data = scores.df, aes(x = scores.df$PC1, y = scores.df$PC2, colour = VRint, shape = VRint), colour = "black", size = 4) +
  geom_point(data = scores.df, aes(x = scores.df$PC1, y = scores.df$PC2, colour = VRint, shape = VRint), size = 3) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_vline(xintercept = -55, size= 0.5) +
  geom_hline(yintercept = 18, size = 0.5) +
  labs(x = "PC1", y = "PC2") +
  theme(axis.title = element_text(size = 20, face = "bold")) +
  theme(axis.text = element_text(size = 18)) +
  geom_point(data = mean_scores, aes(x = mean_scores$PC1, y = mean_scores$PC2, colour = VRint_mean, shape = VRint_mean), colour = "black", size = 9) +
  geom_point(data = mean_scores, aes(x = mean_scores$PC1, y = mean_scores$PC2, colour = VRint_mean, shape = VRint_mean), size = 8) +
  scale_shape_manual(values = c(16,16,17,17), name = "Vaccine/Route", labels = (c("FPV/i.m.", "FPV/i.n.", "Ad5/i.m.", "Ad5/i.n."))) +
  scale_colour_manual("Vaccine/Route", values = c(2, 3, 4, "grey"), labels = (c("FPV/i.m.", "FPV/i.n.", "Ad5/i.m.", "Ad5/i.n."))) +
  theme(legend.position = "bottom") +
  theme(legend.key.size = unit(0.3, "cm")) + 
  theme(legend.background = element_rect(colour = "white", fill= "white", size = 0.1), legend.key = element_rect(colour = "white", fill = "white"), legend.text = element_text(face = "bold", size = 16), legend.title = element_text(face = "bold", size = 16))

p2a <- ggplot(data = scores.df, aes(x = scores.df$PC1, y = scores.df$PC3, colour = VRint, shape = VRint)) + 
  geom_point(data = scores.df, aes(x = scores.df$PC1, y = scores.df$PC3, colour = VRint, shape = VRint), colour = "black", size = 4) +
  geom_point(data = scores.df, aes(x = scores.df$PC1, y = scores.df$PC3, colour = VRint, shape = VRint), size = 3) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_vline(xintercept = -55, size= 0.5) +
  geom_hline(yintercept = 5, size = 0.5) +
  labs(x = "PC1", y = "PC3") +
  theme(axis.title = element_text(size = 20, face = "bold")) +
  theme(axis.text = element_text(size = 18)) +
  geom_point(data = mean_scores, aes(x = mean_scores$PC1, y = mean_scores$PC3, colour = VRint_mean, shape = VRint_mean), colour = "black", size = 9) +
  geom_point(data = mean_scores, aes(x = mean_scores$PC1, y = mean_scores$PC3, colour = VRint_mean, shape = VRint_mean), size = 8) +
  scale_shape_manual(values = c(16,16,17,17), name = "Vaccine/Route", labels = (c("FPV/i.m.", "FPV/i.n.", "Ad5/i.m.", "Ad5/i.n."))) +
  scale_colour_manual("Vaccine/Route", values = c(2, 3, 4, "grey"), labels = (c("FPV/i.m.", "FPV/i.n.", "Ad5/i.m.", "Ad5/i.n."))) +
  theme(legend.position = "bottom") +
  theme(legend.key.size = unit(0.3, "cm")) + 
  theme(legend.background = element_rect(colour = "white", fill= "white", size = 0.1), legend.key = element_rect(colour = "white", fill = "white"), legend.text = element_text(face = "bold", size = 16), legend.title = element_text(face = "bold", size = 16))

tiff('ILC2 group ggplot Route PC1 vs PC2.tiff', units="in", width=7.5, height=7.5, res=600)

p1a

dev.off()

tiff('ILC2 group ggplot Route PC1 vs PC3.tiff', units="in", width=7.5, height=7.5, res=600)

p2a

dev.off()
