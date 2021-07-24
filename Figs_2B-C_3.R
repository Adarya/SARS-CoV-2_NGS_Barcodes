library(ggplot2)
library(reshape2)

###-----
#Figs 2B, 2C

path <- "C:/Users/AdarY/OneDrive - mail.huji.ac.il/SARS-CoV-2/96pos/"

pos96 <- read.csv(paste0(path, "summary_R1_001.csv")) 
pos96 = pos96[pos96$sample != "unmatched",]

#Fig 2B
pdf("plots/pos96_boxplot.pdf")
ggplot(pos96)+
  geom_boxplot(aes(x=Name, y=log2(NumReads), fill = Name))+
  scale_y_continuous(limits = c(0,15))+
  scale_x_discrete(name = "Reads across 96 barcodes")+
  theme_classic()+
  theme(
    axis.title = element_text(size = 20, face = "bold"),
    axis.title.x = element_blank(),
    axis.text = element_text(size = 15, face= "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size =15, face = "bold"),
    legend.position = "none",
  )
dev.off()


#Fig 2C
pdf("plots/pos96_bars.pdf", width = 14)
ggplot(pos96)+
  geom_col(aes(x=sample, y=log2(NumReads), fill = Name))+
  theme(
    axis.text.x = element_text(size = 10, face = "bold", angle = 90, vjust = 0.5),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 15, face = "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size =12, face = "bold")
  )
dev.off()


####----
# Fig 3
mult <- read.csv("data/mult2.csv",header = T)
colnames(mult)[1] <- "Sample"
mult<-melt(mult[mult$status !=0.5,], measure.vars = c(5,4,3,2), value.name = "Reads")

RNASE <- "#C77CFF"
E <- "#00BFC4"
N1 <- "#7CAE00"
N2 <- '#F8766D'
ORF <- "#3C5488B2"

mult2_plot <- ggplot(mult)+
  geom_col(aes(x=reorder(Sample, status),y=log2(Reads+1), fill = variable))+ 
  #ggtitle("S10 L001") +
  scale_x_discrete(name = NULL)+
  scale_fill_manual(values = c(N2,N1,ORF,RNASE))+
  scale_y_continuous(limits = c(0,52))+
  theme(
    axis.text.x = element_text(size = 10, face = "bold", angle = 90, vjust = 0.5),
    title = element_text(size = 10, face="bold"),
    axis.title.x = element_text(size = 11, face = "bold"),
    axis.title.y = element_text(size = 15, face = "bold"),
    axis.text.y = element_text(size = 9),
    legend.title = element_blank(),
    legend.position = c(0.13,0.8),
    legend.key.size = unit(0.8, 'cm'),
    legend.text = element_text(size =12, face = "bold")
  )

pdf("plots/mult2 log2.pdf", width = 12)
mult2_plot
dev.off()
