rm(list=ls())
library(Biostrings ,verbose=FALSE ,warn.conflicts=FALSE, quietly = TRUE)
library(ggplot2 ,verbose=FALSE ,warn.conflicts=FALSE, quietly = TRUE)
library(Cairo ,verbose=FALSE ,warn.conflicts=FALSE, quietly = TRUE)


ChIP_file <- "../results/bs_ChIP"
DAP_file <- "../results/bs_DAP"
out <-  "../results/boxplot.png"
    
ChIP <- read.table(ChIP_file,sep="\t")
DAP <- read.table(DAP_file,sep="\t")

tab <- data.frame(rbind(cbind(ChIP$V8,rep("ChIP-Seq",dim(ChIP)[1])),cbind(DAP$V8,rep("DAP-Seq",dim(DAP)[1]))))
colnames(tab) <- c("Scores","exp")
tab$Scores <- as.numeric(as.character(tab$Scores))

student <- t.test(ChIP$V8,DAP$V8,alternative="less")$p.val

Cairo(width = 600, height = 500, file=out, type="png", pointsize=12, bg = "transparent", canvas = "white", units = "px", dpi = "auto")
g <- ggplot(tab, aes(x=exp, y=Scores,fill=exp)) + geom_boxplot() + 
    theme_classic() +
    scale_fill_manual(values=c("#0D7294", "#40A5C7")) +
annotate(size=8,geom="text",x=1.5,y=1.1,label=paste(" - log(p) ~ ",as.character(-floor(log(student,10))))) +
    scale_y_continuous(breaks=c(0,0.2,0.4,0.6,0.8,1)) +
theme(axis.title.x=element_blank(),axis.title.y=element_text(size=20),axis.text=element_text(size=18),legend.position="none")
g
dev.off()

