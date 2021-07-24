library(ggplot2)
library(gridExtra)
library(stringr)
library(reshape2)
library(ggpubr)
library(readxl)

# Figs 4A, 4B

#A Summary file example can be found at the data dir 
path <- "C:/Users/AdarY/OneDrive - mail.huji.ac.il/SARS-CoV-2/960/Summary_files/"
summary_files <- list.files(path, pattern = paste0("*.csv"))
all_plates <- NULL
for(f in summary_files){
  temp <- read.csv(paste0(path, f), header = T)
  temp$lane <- substr(f,9,16)
  temp$plate <- str_split(f,"_")[[1]][2]
  temp <- temp[temp$sample != 'unmatched',]
  assign(substr(f,9,16), temp)
  all_plates <- rbind(all_plates, temp)
}
colnames(all_plates) <- colnames(get(substr(f,9,16)))

grids <- list()
for(i in 1:length(summary_files)){
  grids[[i]] <- ggplot(get(substr(summary_files[i],9,16))[get(substr(summary_files[i],9,16))$Name != 'RNASE_P',])+
    geom_col(aes(x=sample,y=NumReads, fill = sample))+ 
    ggtitle(paste(str_split(summary_files[i],"_")[[1]][2],
                  str_split(summary_files[i],"_")[[1]][3])) +
    scale_x_discrete(name = "Barcode")+
    scale_y_continuous(limits = c(0,95000))+
    theme(
      axis.text.x = element_blank(),
      title = element_text(size = 18, face="bold"),
      axis.title.x = element_text(size = 15, face = "bold"),
      axis.title.y = element_text(size = 15, face = "bold"),
      axis.text.y = element_text(size = 10),
      legend.position = "None"
    )
}

png("plots/plates and lanes viral reads.png", 1600, 900)
do.call(grid.arrange, c(grids,ncol = 8))
dev.off()

#combine lanes in each plate
for(i in 1:length(unique(all_plates$plate))){
  if(unique(all_plates$plate)[i] == "S2"){ 
    temp <- get(paste0("S2","_L001_"))
    temp$NumReads <- temp$NumReads + 
      get(paste0("S2","_L003_"))$NumReads+
      get(paste0("S2","_L004_"))$NumReads
    
    assign(paste0("S2","_combined"), temp)
    
  }
  else if(unique(all_plates$plate)[i] == "S10"){
    temp <- get(paste0("S10","_L001"))
    temp$NumReads <- temp$NumReads + 
      get(paste0("S10","_L002"))$NumReads+
      get(paste0("S10","_L003"))$NumReads+
      get(paste0("S10","_L004"))$NumReads
    
    assign(paste0(unique(all_plates$plate)[i]),"_combined", temp)
  }
  else{
    temp <- get(paste0(unique(all_plates$plate)[i],"_L001_"))
    temp$NumReads <- temp$NumReads + 
      get(paste0(unique(all_plates$plate)[i],"_L002_"))$NumReads+
      get(paste0(unique(all_plates$plate)[i],"_L003_"))$NumReads+
      get(paste0(unique(all_plates$plate)[i],"_L004_"))$NumReads
    
    assign(paste0(unique(all_plates$plate)[i],"_combined"), temp)
  }
}

#One pool example - Fig 4A
pdf("plots/S7 log.pdf", width = 12)
ggplot(S7_combined[S7_combined$Name %in% c('RNASE_P','E_gene', 'N1_gene', 'ORF1a'),])+
  geom_col(aes(x=sample,y=log(NumReads+1), fill = Name))+ 
  ggtitle("S7") +
  scale_x_discrete(name = "Barcode")+
  #scale_y_continuous(limits = c(0,65000))+
  theme(
    axis.text.x = element_text(size = 11, face = "bold", angle = 90, vjust = 0.5),
    title = element_text(size = 12, face="bold"),
    axis.title.x = element_text(size = 15, face = "bold"),
    axis.title.y = element_text(size = 15, face = "bold"),
    axis.text.y = element_text(size = 12),
    legend.position = c(0.9,0.85),
    legend.title = element_blank(),
    legend.key.size = unit(0.8, 'cm'),
    legend.text = element_text(size =12, face = "bold")
  )
dev.off()

###----
#Fig 4B
all_plates_matrix <- NULL
for(i in 1:10){
  all_plates_matrix <- rbind(all_plates_matrix, 
                             read_xlsx("data/All_pools_results_matrix.xlsx", sheet = i)[-6])
}

all_plates_matrix$log_viral <- log10(all_plates_matrix$E_gene + 
                                       all_plates_matrix$N1_gene +
                                       all_plates_matrix$ORF1a)
  
all_plates_matrix$ct_40 <- all_plates_matrix$N 
for(i in 1:nrow(all_plates_matrix)){
  if(is.na(all_plates_matrix$ct_40[i])){
    all_plates_matrix$ct_40[i] <- 40
  }
}

a<-ggplot(all_plates_matrix[all_plates_matrix$RNA_sample == 1,],aes(x= N, y = log_viral))+
  geom_point(size = 2, alpha = 0.75)+ 
  geom_smooth(method='lm', formula= y~x)+
  #geom_vline(xintercept = 30, type = "dashed")+
  scale_y_continuous(name = "Viral genes reads (log)", limits = c(0,6))+
  scale_x_continuous(name = "Positive samples") + 
  theme_classic()+
  theme(
    axis.title.x = element_text(size = 15, face = "bold"),
    axis.text = element_text(size = 15, face = "bold"),
    axis.title.y = element_text(size = 15, face = "bold")
  )

b<-ggplot(all_plates_matrix[all_plates_matrix$RNA_sample == 0,], 
       aes(x=ct_40,y=log_viral))+
  geom_jitter(size = 2, alpha = 0.75)+
  scale_x_continuous(name = "Negative samples", expand = c(0, 0))+
  scale_y_continuous(limits = c(0,6))+theme_classic()+
  theme(axis.title.y = element_blank(),
        axis.title.x = element_text(size = 15, face = "bold"),
        axis.text.y = element_text(size = 15, face = "bold"),
        axis.text.x = element_text(size = 15, face = "bold", colour = "white"))

pdf("plots/Ct vs reads.pdf", width = 6, height = 4)
gridExtra::grid.arrange(a,b,ncol = 2)
dev.off()

