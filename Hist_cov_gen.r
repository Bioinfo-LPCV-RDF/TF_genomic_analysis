rm(list=ls())
library(Biostrings ,verbose=FALSE ,warn.conflicts=FALSE, quietly = TRUE)
library(ggplot2 ,verbose=FALSE ,warn.conflicts=FALSE, quietly = TRUE)
library(Cairo ,verbose=FALSE ,warn.conflicts=FALSE, quietly = TRUE)

args=commandArgs(trailingOnly=TRUE)
dir=args[1]
name1=args[2]
name2=args[3]
tableName=args[4]
rationame1=as.numeric(args[5]) # >1
rationame2=as.numeric(args[6]) # <1

# Rscript Hist_cov_gen.r dir name1 name2 tablename rationame1 rationame2

if (rationame1 < rationame2){ # if ratios are not in good order, just exchange them
	rationame1 = rationame1 + rationame2
	rationame2 = rationame1 - rationame2
	rationame1 = rationame1 - rationame2
}
print(rationame1)
print(rationame2)
out_dir <- dir
#in_dir <- "../../results/test_pipeline_all_size"

dir.create(file.path(out_dir),recursive=TRUE ,showWarnings = FALSE)

table <- read.table(tableName,sep="\t",header=TRUE)
len <- dim(table)[1]
## nbR <- len %/% frac
bed <- table[,(1:3)]
#vec <- unname(unlist(table["SEP3delAG"]/table["SEP3AG"]))
vec <- unname(unlist(table[name1]/table[name2]))


bed$ratio = vec
bed$group = "Common"
bed$group[bed$ratio > rationame1] = paste(name1,"specific",sep=" ")
bed$group[bed$ratio < rationame2] =  paste(name2,"specific",sep=" ")
bed$group <- as.factor(bed$group)
print(paste(name1," vs ", name2,sep=""))
summary(bed$group)

Cairo(width = 600, height = 500, file=paste(out_dir,"/plots/ratio_coverages_",name1,"_",name2,"_ALLDAPpeaks.png",sep=""), type="png", pointsize=12, bg = "transparent", canvas = "white", units = "px", dpi = "auto")
ggplot(bed, aes(x=bed$ratio, color=bed$group)) + 
geom_histogram(fill="white",binwidth=0.05) + 
scale_x_log10() + 
xlab(paste("ratio Cov ",name1," / Cov ",name2))
dev.off()


Cairo(width = 1500, height = 1000, file=paste(out_dir,"/plots/covplot_",name1,"_",name2,".png",sep=""), type="png", pointsize=12, bg = "transparent", canvas = "white", units = "px", dpi = "auto")
g <- ggplot(data=table,aes(y=table[,name1],x=table[,name2]))
g + geom_point(alpha=0.4,aes(color=as.character(name))) + scale_x_log10() + scale_y_log10() + labs(color="peaks") + theme(axis.title=element_text(size=25),legend.text=element_text(size=25)) + guides(colour = guide_legend(override.aes = list(size=10,alpha=1))) +  geom_abline(intercept = 0) + xlab (paste("cov ",name2,sep="")) + ylab (paste("cov ",name1,sep=""))
dev.off()


Name1glob_temp=table[table$name==paste(name1,sep=" ")|table$name==paste("both",sep=" "),]
Name1glob<-Name1glob_temp[order(Name1glob_temp[,name1],decreasing=TRUE),]
write.table(Name1glob,file.path(out_dir,paste(name1,".bed",sep="")),quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
Name2glob_temp=table[table$name==paste(name2,sep=" ")|table$name==paste("both",sep=" "),]
Name2glob<-Name2glob_temp[order(Name2glob_temp[,name2],decreasing=TRUE),]
write.table(Name2glob,file.path(out_dir,paste(name2,".bed",sep="")),quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)



table$ratio = vec
table$group = "Common"
table$group[table$ratio > rationame1] = paste(name1,"specific",sep=" ")
table$group[table$ratio < rationame2] =  paste(name2,"specific",sep=" ")
table$group <- as.factor(table$group)

Cairo(width = 1500, height = 1000, file=paste(out_dir,"/plots/covplot_",name1,"_",name2,"_NewGroups.png",sep=""), type="png", pointsize=12, bg = "transparent", canvas = "white", units = "px", dpi = "auto")
g <- ggplot(data=table,aes(y=table[,name1],x=table[,name2]))
g + geom_point(alpha=0.4,aes(color=group)) + scale_x_log10() + scale_y_log10() + labs(color="peaks") + theme(axis.title=element_text(size=25),legend.text=element_text(size=25)) + guides(colour = guide_legend(override.aes = list(size=10,alpha=1))) +  geom_abline(intercept = 0) + xlab (paste("cov ",name2,sep="")) + ylab (paste("cov ",name1,sep=""))
dev.off()


Name1Spe_temp=table[table$group==paste(name1,"specific",sep=" "),]
Name1Spe<-Name1Spe_temp[order(Name1Spe_temp[,name1],decreasing=TRUE),]
Name2Spe_temp=table[table$group==paste(name2,"specific",sep=" "),]
Name2Spe<-Name2Spe_temp[order(Name2Spe_temp[,name2],decreasing=TRUE),]
Common_temp=table[table$group=="Common",]
Common_temp$mean=(Common_temp[,name1]+Common_temp[,name2])/2
Common<-Common_temp[order(Common_temp$mean,decreasing=TRUE),]
Common$mean<-NULL
write.table(Name1Spe,file.path(out_dir,paste(name1,"_spePeaks_",rationame1,".bed",sep="")),quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
write.table(Name2Spe,file.path(out_dir,paste(name2,"_spePeaks_",rationame2,".bed",sep="")),quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
write.table(Common,file.path(out_dir,paste(name1,"_",name2,"_commonPeaks_",rationame1,"_",rationame2,".bed",sep="")),quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)


