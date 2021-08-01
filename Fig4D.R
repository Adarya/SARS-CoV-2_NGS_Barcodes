library(ggplot2)
library(stringr)
library(vcfR)

#bcftools outputs are as .vcf files
path <- "C:/Users/Yoga920/OneDrive - mail.huji.ac.il/Barcode/28_01_2021/vcf/"
vcf_files <- list.files(path, pattern = paste0("*.vcf$"))
                  
vcf_reader <- function(path, vcf_files){
  vcf_df <- data.frame()
  for(i in 1:length(vcf_files)){
    t = read.vcfR(paste0(path,vcf_files[i]))
    if(length(getPOS(t)) > 0){
      vcf_df <- rbind(vcf_df,
                      cbind(ID = substr(vcf_files[i],1,10),
                            CHROM = getCHROM(t), POS = getPOS(t), REF = getREF(t), 
                            ALT = getALT(t), QUAL = getQUAL(t), INFO2df(t)))
    }
  }
  return(as.data.frame(vcf_df))
}

vcf_df <- vcf_reader(path, vcf_files)
vcf_df <- vcf_df[vcf_df$DP >= 50, ]

write.csv(vcf_df, "data/unannotated_filtered_calls.csv")

#un-annotated data were annotated and curated manually for consequence and sample ID SARS-CoV-2 status 
mut_analysis <- read.csv("data/bcftools_variants_annotated.csv",header = T)[1:88,]

p <- ggplot(mut_analysis[mut_analysis$CONSEQU != "none",], (aes(x=as.factor(POS), y=stat(count), col = STATUS, 
                          shape = CONSEQU)))+
  stat_count(geom = "point", size = 5)+
  theme_bw()+
  scale_x_discrete(name = NULL, labels = unique(paste0(mut_analysis$POS[mut_analysis$CONSEQU != "Artefact"], " ", mut_analysis$REF, ">", mut_analysis$ALT)))+
  scale_shape(name = "Mutation type")+
  scale_color_discrete(name = "Status", labels = c("Negative","Positive"))+
  scale_y_continuous(name = "Number of mutations", breaks = seq(0,30,3))+
  theme(
    axis.text.x = element_text(size =15,face = "bold", angle = 90, 
                               vjust = 0.5, colour = c("red","skyblue2",
                                                       "skyblue2", rep("black",10))),
    axis.text.y = element_text(size=15,face = "bold"),
    axis.title = element_text(size = 20, face = "bold"),
    legend.title = element_text(size = 17, face = "bold"),
    legend.text = element_text(size = 17, face = "bold")
  )

pdf("plots/mutations.pdf", width = 7, height = 5)
p
dev.off()
