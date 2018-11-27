##This is specifically to plot the data using ggplot2 for visualising differences between vectors following i.m. vaccination
scores.df <- as.data.frame(scores)

scores.df[, c("Vaccine", "Route", "Experiment")] <- data[, c("Vaccine", "Route", "Experiment")]

library(ggplot2)
 
mean_scores <- aggregate(scores.df[, 1:4], list(scores.df$Vaccine, scores.df$Route), mean)

names(mean_scores)[names(mean_scores) == "Group.1"] <- "Vaccine"
names(mean_scores)[names(mean_scores) == "Group.2"] <- "Route"

p1 <- ggplot(data = scores.df, aes(x = scores.df$PC1, y = scores.df$PC4, colour = Vaccine, shape = Vaccine)) + 
  geom_point(data = scores.df, aes(x = scores.df$PC1, y = scores.df$PC4, colour = Vaccine, shape = Vaccine), colour = "black", size = 4) +
  geom_point(data = scores.df, aes(x = scores.df$PC1, y = scores.df$PC4, colour = Vaccine, shape = Vaccine), size = 3) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_vline(xintercept = 30, size= 0.5) +
  geom_hline(yintercept = 5, size = 0.5) +
  scale_shape_manual(values = c(16,18,17), name = "Vaccine", labels = (c("FPV/i.m.", "MVA/i.m.", "Ad5/i.m."))) +
  scale_colour_manual("Vaccine", values = c(2,'yellow',3), labels = c("FPV/i.m.", "MVA/i.m.","Ad5/i.m.")) +
  labs(x = "PC1", y = "PC4") +
  theme(axis.title = element_text(size = 20, face = "bold")) +
  theme(axis.text = element_text(size = 18)) +
  geom_point(data = mean_scores, aes(x = mean_scores$PC1, y = mean_scores$PC4, colour = Vaccine), colour = "black", size = 9) +
  geom_point(data = mean_scores, aes(x = mean_scores$PC1, y = mean_scores$PC4, colour = Vaccine), size = 8) +
  theme(legend.position = "bottom") +
  theme(legend.key.size = unit(0.3, "cm")) + 
  theme(legend.background = element_rect(colour = "white", fill= "white", size = 0.1), legend.key = element_rect(colour = "white", fill = "white"), legend.text = element_text(face = "bold", size = 16), legend.title = element_text(face = "bold", size = 16))

tiff('ILC2 grouped ggplot im only.tiff', units="in", width=7.5, height=7.5, res=600)

p1

dev.off()

