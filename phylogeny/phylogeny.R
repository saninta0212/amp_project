library(ggplot2)
library(tidyverse)
library(ape)
library(ggtree)

# Read the nwk file (where the tree is given as parenthetic format) and the csv file with the TPM value
tree1 <- read.tree("/projects/btl_scratch/saninta/amp_pipeline/phylogeny/amp.mfft.anura.nwk",text = NULL, tree.names = NULL, skip = 0,
                  comment.char = "", keep.multi = FALSE)

TPM_data <- read.csv("/projects/btl_scratch/saninta/amp_pipeline/phylogeny/filtered.csv",sep=",",header=TRUE)

TPM_data <- cbind(TPM_data, lnTPM=sort(log(TPM_data$TPM)))
TPM_data <- cbind(TPM_data, Label=1:nrow(TPM_data))
tip_dir <- as.character(TPM_data$Name)
tip_dir <- cbind(tip_dir,lnTPM=TPM_data$lnTPM)
tip_dir <- as.data.frame(tip_dir)
tip_dir <- cbind(tip_dir,Transcription_Levels= lnTPM)

tip_dir$lnTPM <- as.numeric(lnTPM)
tip_dir$Transcription_Levels <- "Low"
tip_dir$Transcription_Levels[tip_dir$lnTPM >2 & tip_dir$lnTPM<4] <- "Medium"
tip_dir$Transcription_Levels[tip_dir$lnTPM > 4] <- "High"

# Plot the tree

# Make a bar plot with the TPM values 
make_bar <- function() {
  bar <- ggplot(TPM_data,aes(x=Label,y=lnTPM)) + geom_col(fill="lightblue") + xlab("Labels of sequences") + 
    ylab("lnTPM (Transcriptis Per Kilobase Million)") + ggtitle("Expression levels of AMPs") +
    theme_light() + labs(caption="*The sequence names are represented as discrete labels") +
    theme(axis.text.x = element_text(angle=90, hjust=0.95,vjust=0.2,size=12),
          text = element_text(size = 12),
          plot.title = element_text(hjust = 0.5,size=14), plot.caption = element_text(face = "italic"))  
  ggsave("/projects/btl_scratch/saninta/amp_pipeline/phylogeny/TPM.bars.png",
         bar, units="in", height=4, width=9, dpi=300)
}

amp_tree <- ggtree(tree1) %<+% tip_dir + # Plot default tree
  geom_tiplab(aes(fill=factor(Transcription_Levels, levels = c("Low", "Medium", "High"))), color = "black", geom = "label", size=1, label.padding = unit(0.15, "lines"),label.size=0) + 
  ggtitle("Expression characterization of AMPs") + theme_tree2() + labs(fill="Transcription levels") +
  theme(plot.title = element_text(hjust = 0.5,face="bold",size=14)) 

  ggsave("/projects/btl_scratch/saninta/amp_pipeline/phylogeny/amp_tree.png", amp_tree, device = "png", width=50, height=100, units = "cm", dpi = 300)
  ggsave("/projects/btl_scratch/saninta/amp_pipeline/phylogeny/amp_tree.pdf", amp_tree, device = "pdf", width=50, height=100, units = "cm", dpi = 300)