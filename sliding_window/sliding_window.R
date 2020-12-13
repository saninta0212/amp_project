# Read the TSV files

slid_seq30 <- read.table("/projects/btl_scratch/saninta/amp_pipeline/anura/apetersi/skin_liver/amplify/slid_seq/AMPlify_results.tsv",sep = "\t", header = TRUE)
slid_seq40 <- read.table("/projects/btl_scratch/saninta/amp_pipeline/anura/apetersi/skin_liver/amplify/slid_seq/40mer/AMPlify_results.tsv",sep = "\t", header = TRUE)
slid_seq50 <- read.table("/projects/btl_scratch/saninta/amp_pipeline/anura/apetersi/skin_liver/amplify/slid_seq/50mer/AMPlify_results.tsv",sep = "\t", header = TRUE)
protein  <- read.table("/projects/btl_scratch/saninta/amp_pipeline/fasta.tsv",sep="\t",header=TRUE)

# Plot functions of amplify score vs window positions

make_plots <- function(sequence_id1,sequence_id2,kmer,dataset) {
  sequence1 <- filter(dataset, str_detect(Sequence_ID,sequence_id1))
  sequence1 <- cbind(sequence1, ID=sequence_id1)
  sequence1 <- cbind(sequence1, Starting_window=1:nrow(sequence1)) 
  sequence2 <- filter(dataset, str_detect(Sequence_ID,sequence_id2))
  sequence2 <- cbind(sequence2, ID=sequence_id2)
  sequence2 <- cbind(sequence2, Starting_window=1:nrow(sequence2)) 
  sequence <- rbind(sequence1,sequence2)
  score <- factor(sequence$Score) 
  amp <- as.numeric(as.character(score)) 
  len1 <- filter(protein, str_detect(Sequence_ID,sequence_id1))
  len1 <- as.character(len1$Length)
  len2 <- filter(protein, str_detect(Sequence_ID,sequence_id2))
  len2 <- as.character(len2$Length)
  len <- paste("Lengths of the peptide precursors are : ",len1,len2,sep=" ")
  kk <- paste("Length of window size: ",kmer,sep="")
  x <- sequence$Starting_window
  plot <- ggplot(sequence,aes(x=Starting_window,y=amp,color=Charge)) + xlab("Starting position of the sliding window") + 
    ylab("AMPlify score") + geom_point(aes(fill=Charge),color="black",shape=21,size=2) + facet_grid(~ID,scales="free_x") +
    geom_hline(yintercept=0.5,linetype="dotdash")+ geom_hline(yintercept=0.99,linetype="dotdash") + theme_bw() +labs(subtitle=kk, caption=len) + 
    scale_fill_gradient2(low="#f43a4a", mid = "white",high="#a346f4") + theme(strip.text = element_text(size=10)) + 
    theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))
  ggsave("/projects/btl_scratch/saninta/amp_pipeline/plot.png",
         plot, units="in", height=4, width=9, dpi=300)
  return(plot) 
  
}

# Makes line plots of AMPs
make_line <- function(sequence_id1,sequence_id2,kmer,dataset) {
  sequence1 <- filter(dataset, str_detect(Sequence_ID,sequence_id1))
  sequence1 <- cbind(sequence1, ID=sequence_id1)
  sequence1 <- cbind(sequence1, Starting_window=1:nrow(sequence1)) 
  sequence2 <- filter(dataset, str_detect(Sequence_ID,sequence_id2))
  sequence2 <- cbind(sequence2, ID=sequence_id2)
  sequence2 <- cbind(sequence2, Starting_window=1:nrow(sequence2)) 
  sequence <- rbind(sequence1,sequence2)
  score <- factor(sequence$Score) 
  amp <- as.numeric(as.character(score)) 
  len1 <- filter(protein, str_detect(Sequence_ID,sequence_id1))
  len1 <- as.character(len1$Length)
  len2 <- filter(protein, str_detect(Sequence_ID,sequence_id2))
  len2 <- as.character(len2$Length)
  len <- paste("Lengths of the peptide precursors are :",len1,len2,sep=" ")
  kk <- paste("Length of window size: ",kmer,sep="")
  plot <- ggplot(sequence,aes(x=Starting_window,y=amp)) + xlab("Starting position of the sliding window") + 
    ylab("AMPlify score") + geom_line(color="darkgreen") + facet_grid(~ID,scales="free_x") +
    geom_hline(yintercept=0.5,linetype="dotdash") + theme_light() + labs(subtitle=kk, caption=len)
  ggsave("/projects/btl_scratch/saninta/amp_pipeline/plot-line.png",
         plot, units="in", height=4, width=9, dpi=300)
  return(plot) 
  
}

# Plot a single AMP
make_splot <- function(sequence_id,kmer,dataset) {
  seq <- filter(dataset, str_detect(Sequence_ID,sequence_id))
  seq <- cbind(seq, Starting_window=1:nrow(seq)) 
  score <- factor(seq$Score) 
  amp <- as.numeric(as.character(score)) 
  kk <- paste("Length of window size: ",kmer,sep="")
  len1 <- filter(protein, str_detect(Sequence_ID,sequence_id))
  len1 <- as.character(len1$Length)
  len <- paste("Length of the peptide precursor is :",len1,sep=" ")
  plot <- ggplot(seq,aes(x=Starting_window,y=amp,fill=Charge)) + xlab("Starting position of the sliding window") + 
    ylab("AMPlify score") + geom_point(aes(fill=Charge),color="black",shape=21,size=2) + 
    geom_hline(yintercept=0.5,linetype="dotdash")+geom_hline(yintercept=0.99,linetype="dotdash")+ geom_hline(yintercept=0.99,linetype="dotdash") + theme_bw() + labs(title=sequence_id,subtitle=kk, caption=len) +
    scale_fill_gradient2(low="#f43a4a",mid = "white",high="#a346f4") + theme(strip.text = element_text(size=12)) + theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))
  ggsave("/projects/btl_scratch/saninta/amp_pipeline/plot_single.png",
         plot, units="in", height=4, width=9, dpi=300)
  return(plot) 
}

# Plot acf 

calculate_acf <- function(sequence_id,dataset) {
  sequence <- filter(dataset, str_detect(Sequence_ID,sequence_id))
  AMPlify_score <- sequence$Score
  acf_p <- ggAcf(AMPlify_score,lag.max=100,type=c("correlation"),plot=TRUE) + theme_bw() +   
   theme(axis.text=element_text(size=12)) + theme(axis.title=element_text(size=14,face="bold")) 
  ggsave("/projects/btl_scratch/saninta/amp_pipeline/kamp-acf.png",
         acf_p, units="in", height=4, width=9, dpi=300)
  
  return(acf_p)
}

# Makes a grid with the sequence AMPlify plot with its acf plot 

arrange_acf <- function(pp,acfp) {
    legend <- get_legend(pp)   # get the legend of the first one plot
    
    # here the plots in a grid
    prow <- plot_grid( pp + theme(legend.position="none"),
                       # here you add the percentage
                       acfp + theme(legend.position="none"),
                       align = 'v',
                       hjust = -1,
                       nrow = 2)
    
    # here you add the legend
    p <- plot_grid(prow, legend, rel_widths = c(3, .3))
    ggsave("/projects/btl_scratch/saninta/amp_pipeline/kamp-acf+plot.png",
           plot=p, units="in", height=4, width=9, dpi=300)
    
}


# Plot the known AMPs

# Reading the files with mature length cleaving


  HP2_m <- read.table("/projects/btl_scratch/saninta/amp_pipeline/sliding_window/known_amp/amplify/29mer/AMPlify_results.tsv",sep="\t",header=TRUE)
  HP2_m <- cbind(HP2_m,actual_name="HP2")
  HP2_m <- cbind(HP2_m, Starting_window=1:nrow(HP2_m)) 
  HP2_m <- cbind(HP2_m,kmer="Mature_seq: 29mer")
  HP4_m <- read.table("/projects/btl_scratch/saninta/amp_pipeline/sliding_window/known_amp/amplify/25mer/AMPlify_results.tsv",sep="\t",header=TRUE)
  HP4_m <- cbind(HP4_m,actual_name="HP4")
  HP4_m <- cbind(HP4_m, Starting_window=1:nrow(HP4_m)) 
  HP4_m <- cbind(HP4_m,kmer="Mature_seq: 25mer")
  HP9_m <- read.table("/projects/btl_scratch/saninta/amp_pipeline/sliding_window/known_amp/amplify/31mer/AMPlify_results.tsv",sep="\t",header=TRUE)
  HP9_m <- cbind(HP9_m,actual_name="HP9")
  HP9_m <- cbind(HP9_m, Starting_window=1:nrow(HP9_m)) 
  HP9_m <- cbind(HP9_m,kmer="Mature_seq: 31mer")


# Filter specific sequences and make a table with its actual name and starting window as a column

filter_all_k <- function(dataset,seq_id,name,kmer) {
  data_k <- filter(dataset, str_detect(Sequence_ID,seq_id))
  data_k <- cbind(data_k,actual_name=name)
  data_k <- cbind(data_k, Starting_window=1:nrow(data_k))
  data_k <- cbind(data_k,kmer)
  return(data_k) 
  
}

# Combine all 6 known amp cleaving sizes in one table (helps to facet by kmer size later)
combine_all_k <- function(data1,data2,data3,data4,data5,seq_id,name,mature) {
  k_10 <- filter_all_k(data1,seq_id,name,"10mer window")
  k_20 <- filter_all_k(data2,seq_id,name,"20mer window")
  k_30 <- filter_all_k(data3,seq_id,name,"30mer window")
  k_40 <- filter_all_k(data4,seq_id,name,"40mer window")
  k_50 <- filter_all_k(data5,seq_id,name,"50mer window")
  amp  <- rbind(k_10,k_20,k_30,k_40,k_50,mature)
  return(amp)
  
}

# Make a table with all 6 of its cleaved pieces
# hp2 <- combine_all_k(k1,k2,k3,k4,k5,"HP2","HP2: AFR43665.1-Ranatuerin-2PRc",HP2_m)
# hp4 <- combine_all_k(k1,k2,k3,k4,k5,"HP4","HP4: ACR46972.1-Ranatuerin-1",HP4_m)
# hp9 <- combine_all_k(k1,k2,k3,k4,k5,"HP9","HP9: ACR84085.1-Palustrin-Ca",HP9_m)

# Plot each known amp faceted by the kmer sizes

make_kplot <- function(kamp,leng) {
  score <- factor(kamp$Score) 
  amp <- as.numeric(as.character(score))
  len <- paste("Length of the peptide precursor is: ",leng,sep="")
  plot <- ggplot(kamp,aes(x=Starting_window,y=amp,fill=Charge)) + xlab("Starting position of the sliding window") + 
    ylab("AMPlify score") + geom_point(aes(fill=Charge),color="black",shape=21,size=2) + 
    geom_hline(yintercept=0.5,linetype="dotdash") + theme_bw() + facet_wrap(~kmer) + 
    scale_fill_gradient2(low="#f43a4a",mid = "white",high="#a346f4") + labs(title=kamp$actual_name,caption=len) +
    theme(strip.text = element_text(size=12)) + theme(axis.text=element_text(size=12),
                                                      axis.title=element_text(size=14,face="bold"))
  
  ggsave("/projects/btl_scratch/saninta/amp_pipeline/kamp-known.png",
         plot, units="in", height=4, width=9, dpi=300)
  
  return(plot) 
}

# Plotting the known AMPs with mature seq cleaving

make_kplot1 <- function() {
  seq <- rbind(HP2_m,HP4_m,HP9_m)
  score <- factor(seq$Score) 
  amp <- as.numeric(as.character(score))
  kplot <- ggplot(seq,aes(x=Starting_window,y=amp,fill=Charge)) + xlab("Starting position of the sliding window") + 
    ylab("AMPlify score") + geom_point(aes(fill=Charge),color="black",shape=21,size=2) + 
    geom_hline(yintercept=0.99,linetype="dotdash") + geom_hline(yintercept=0.5,linetype="dotdash") + theme_bw() + 
    scale_fill_gradient2(low="#f43a4a",mid = "white",high="#a346f4") + facet_wrap(~actual_name,nrow=3) +
    theme(strip.text = element_text(size=12)) + theme(axis.text=element_text(size=12),
                                                      axis.title=element_text(size=14,face="bold"))
  
  ggsave("/projects/btl_scratch/saninta/amp_pipeline/known_amp-mature.png",
         kplot, units="in", height=4, width=9, dpi=300)

}

# Plot for a single known AMP

make_kplot2 <- function(dataset) {
  score <- factor(HP4_m$Score) 
  amp <- as.numeric(as.character(score))
  kplot <- ggplot(HP4_m,aes(x=Starting_window,y=amp,fill=Charge)) + xlab("Starting position of the sliding window") + 
    ylab("AMPlify score") + geom_point(aes(fill=Charge),color="black",shape=21,size=2) + 
    geom_hline(yintercept=0.99,linetype="dotdash") + geom_hline(yintercept=0.5,linetype="dotdash") + theme_bw() + 
    scale_fill_gradient2(low="#f43a4a",mid = "white",high="#a346f4") + 
    theme(axis.title=element_text(size=14,face="bold"))
  ggsave("/projects/btl_scratch/saninta/amp_pipeline/plot-kamp.png",
         kplot, units="in", height=4, width=9, dpi=300)
  
}



# Negative control peptides

neg <- read.table("/projects/btl_scratch/saninta/amp_pipeline/negative_control/amplify/AMPlify_results.tsv",sep="\t",header=TRUE)

make_neg_plot <- function(kamp) {
  kamp<-neg
  albumin <- filter(kamp,str_detect(Sequence_ID,"ALBUB_XENLA"))
  albumin <- cbind(albumin,id="sp|P14872|ALBUB_XENLA")
  albumin <- cbind(albumin, Starting_window=1:nrow(albumin)) 
  collagen <- filter(kamp,str_detect(Sequence_ID,"CO2A1_XENLA"))
  collagen <- cbind(collagen,id="sp|Q91717|CO2A1_XENLA")
  collagen <- cbind(collagen, Starting_window=1:nrow(collagen)) 
  keratin <- filter(kamp,str_detect(Sequence_ID,"K118A_XENLA"))
  keratin <- cbind(keratin,id="sp|P08802|K118A_XENLA")
  keratin <- cbind(keratin, Starting_window=1:nrow(keratin)) 
  neg <- rbind(albumin,collagen,keratin)
  score <- factor(neg$Score) 
  amp <- as.numeric(as.character(score))
  
  plot <- ggplot(neg,aes(x=Starting_window,y=amp,fill=Charge)) + xlab("Starting position of the sliding window") + 
    ylab("AMPlify score") + geom_point(aes(fill=Charge),color="black",shape=21,size=2) + 
    geom_hline(yintercept=0.5,linetype="dotdash") + theme_bw() + facet_wrap(~id,scales="free_x") + 
    scale_fill_gradient2(low="#f43a4a",mid = "white",high="#a346f4") + 
    labs(title="Keratin, Albumin and Collagen from Frog skin") +  theme(strip.text = element_text(size=10)) + 
    theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))
  
  ggsave("/projects/btl_scratch/saninta/amp_pipeline/neg-amp.png",
         plot, units="in", height=4, width=9, dpi=300)
  
  return(plot) 
}

# Making the plots

# make_neg_plot(neg)
# 
# make_plots("apetersi_skin-liver_15980.p1","apetersi_skin-liver_12383.p1",30,slid_seq30)
# make_plots("apetersi_skin-liver_15980.p1","apetersi_skin-liver_24469.p2",40,slid_seq40)
# make_plots("apetersi_skin-liver_49400.p1","apetersi_skin-liver_100789.p1",30,slid_seq30)
# make_plots("apetersi_skin-liver_49400.p1","apetersi_skin-liver_100789.p1",40,slid_seq40)
# make_plots("apetersi_skin-liver_50068.p1","apetersi_skin-liver_15845.p1",30,slid_seq30)
# make_plots("apetersi_skin-liver_50068.p1","apetersi_skin-liver_15845.p1",40,slid_seq40)
# 
# 
# make_splot("apetersi_skin-liver_106.p1",30,slid_seq30)
# make_splot("apetersi_skin-liver_52327.p2",30,slid_seq30)

# Make all the plots from 30mer
# make_splot("apetersi_skin-liver_100789.p1",30,slid_seq30)
# make_splot("apetersi_skin-liver_106.p1",30,slid_seq30)
# 
# 
# make_line("apetersi_skin-liver_15980.p1","apetersi_skin-liver_24469.p2",30,slid_seq30)
# make_line("apetersi_skin-liver_49400.p1","apetersi_skin-liver_100789.p1",30,slid_seq30)
# make_line("apetersi_skin-liver_50068.p1","apetersi_skin-liver_15845.p1",30,slid_seq30)

# Calling functions to make

# hp2_30 <- make_kplot(k_amp,"HP2","AFR43665.1-Ranatuerin-2PRc",43,71,30)
# hp4_30 <- make_kplot(k_amp,"HP4","ACR46972.1-Ranatuerin-1",37,61,30)
# hp9_30 <- make_kplot(k_amp,"HP9","ACR84085.1-Palustrin-Ca",21,51,30)
# make_kplot(hp2,71)
# make_kplot(hp4,61)
# make_kplot(hp9,54)
# 
# arrange_k <- function(p1,p2,p3) {
#   legend <- get_legend(p1)   # get the legend of the first one plot
#   
#   # here the plots in a grid
#   prow <- plot_grid( p1 + theme(legend.position="none"),
#                      p2 + theme(legend.position="none"),
#                      p3 + theme(legend.position="none"),
#                      align = 'v',
#                      hjust = -1,
#                      ncol = 3)
#   
#   # here you add the legend
#   p <- plot_grid(prow, legend, rel_widths = c(3, .3))
#   return (p)
# }
# 
# arrange_k(hp2_30,hp4_30,hp9_30)
# 
# 
# # Start plotting only the sequences that have length < 100
# 
# seq_50 <- read.table("/projects/btl_scratch/saninta/amp_pipeline/anura/apetersi/seq_100/amplify/50/AMPlify_results.tsv",sep="\t",header=TRUE)
# make_plots("apetersi_skin-liver_13087.p1","apetersi_skin-liver_12383.p1",50,seq_50)
# make_plots("apetersi_skin-liver_14142.p1","apetersi_skin-liver_15236.p1",50,seq_50)
# make_plots("apetersi_skin-liver_15980.p1","apetersi_skin-liver_24469.p2",50,seq_50)
# make_plots("apetersi_skin-liver_16922.p1","apetersi_skin-liver_17196.p1",50,seq_50)
# make_plots("apetersi_skin-liver_18134.p1","apetersi_skin-liver_18295.p1",50,seq_50)
# make_plots("apetersi_skin-liver_18346.p1","apetersi_skin-liver_19016.p1",50,seq_50)
# make_plots("apetersi_skin-liver_190.p1","apetersi_skin-liver_1945.p1",50,seq_50)
# make_plots("apetersi_skin-liver_21541.p2","apetersi_skin-liver_22835.p1",50,seq_50)
# make_plots("apetersi_skin-liver_22867.p1","apetersi_skin-liver_23019.p1",50,seq_50)
# make_plots("apetersi_skin-liver_23586.p2","apetersi_skin-liver_24469.p2",50,seq_50)
# make_plots("apetersi_skin-liver_24525.p1","apetersi_skin-liver_2468.p1",50,seq_50)
# make_plots("apetersi_skin-liver_28128.p1","apetersi_skin-liver_28614.p1",50,seq_50)
# make_plots("apetersi_skin-liver_29228.p1","apetersi_skin-liver_30304.p1",50,seq_50)
# make_plots("apetersi_skin-liver_30467.p2","apetersi_skin-liver_31007.p2",50,seq_50)
# make_plots("apetersi_skin-liver_31029.p1","apetersi_skin-liver_31044.p2",50,seq_50)
# make_plots("apetersi_skin-liver_33082.p2","apetersi_skin-liver_33581.p1",50,seq_50)
# make_plots("apetersi_skin-liver_34515.p1","apetersi_skin-liver_36369.p1",50,seq_50)
# make_plots("apetersi_skin-liver_37328.p1","apetersi_skin-liver_37581.p1",50,seq_50)
# make_plots("apetersi_skin-liver_38051.p1","apetersi_skin-liver_39055.p1",50,seq_50)
# make_plots("apetersi_skin-liver_39282.p1","apetersi_skin-liver_39586.p2",50,seq_50)
# make_plots("apetersi_skin-liver_40550.p1","apetersi_skin-liver_40731.p1",50,seq_50)
# make_plots("apetersi_skin-liver_97665.p1","apetersi_skin-liver_.p1",50,seq_50)
# 
# make_plots("apetersi_skin-liver_15849","apetersi_skin-liver_15845",50,seq_50)
# make_plots("apetersi_skin-liver_23586.p2","apetersi_skin-liver_24469.p2",50,seq_50)
# make_plots("apetersi_skin-liver_23586.p2","apetersi_skin-liver_24469.p2",50,seq_50)
# 
# 
# seq_60 <- read.table("/projects/btl_scratch/saninta/amp_pipeline/anura/apetersi/seq_100/amplify/60/AMPlify_results.tsv",sep="\t",header=TRUE)
# 
# 
# 
# make_kplot(HP9_m,54)
