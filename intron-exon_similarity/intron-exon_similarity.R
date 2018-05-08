library(ggplot2)
library(forcats)
library(svglite)

mt<-read.csv("intron_exon_similarity.csv",header=TRUE, sep=",")

p <- ggplot(mt,aes(fct_rev(locus),conservation, colour=purple, size=10)) +
  guides(size=FALSE) +
  geom_jitter(aes(shape=triangle), size=12, width=0.2, height=0.2) +
  scale_shape_manual(values=c(1,1)) +
  coord_flip() +
  theme_classic() +
  theme(legend.text=element_text(size=15)) +
  theme(legend.title=element_blank()) +
  guides(colour=guide_legend(override.aes=list(size=5)))

p + scale_colour_manual(values=c("grey25","grey70")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  expand_limits(y=c(50,100)) +
  xlab("Diatom-to-Chattonella percent similarity") +
  theme(axis.text.x = element_text(colour="grey1", size=15)) +
  theme(axis.title.x = element_text(colour="grey1", size=20)) +
  scale_x_discrete(name="")
