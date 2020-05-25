rm(list=ls())
#library(Biostrings ,verbose=FALSE ,warn.conflicts=FALSE, quietly = TRUE)
library(ggplot2 ,verbose=FALSE ,warn.conflicts=FALSE, quietly = TRUE)
library(Cairo ,verbose=FALSE ,warn.conflicts=FALSE, quietly = TRUE)
library(latex2exp ,verbose=FALSE ,warn.conflicts=FALSE, quietly = TRUE)

#library(gridExtra)

library("optparse")
option_list = list(
  make_option(c("-f", "--files"), type="character", default=NULL, 
              help="tables files ", metavar="character"),
  make_option(c("-n", "--names"), type="character", default=NULL, 
              help=" ", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

myfilelist <- strsplit(opt$files, ",")
name1=strsplit(opt$names, ",")[[1]][1]
name2=strsplit(opt$names, ",")[[1]][2]
#table$MSPC <- "common"
#table$MSPC[table$CFR<=0.125] <- paste(name2,"specific",sep=" ")
#table$MSPC[table$CFR>=8] <- paste(name1,"specific",sep=" ")

#blankPlot <- ggplot()+geom_blank(aes(1,1))+
#  theme(plot.background = element_blank(), 
#   panel.grid.major = element_blank(),
#   panel.grid.minor = element_blank(), 
#   panel.border = element_blank(),
#   panel.background = element_blank(),
#   axis.title.x = element_blank(),
#   axis.title.y = element_blank(),
#   axis.text.x = element_blank(), 
#   axis.text.y = element_blank(),
#   axis.ticks = element_blank()
#     )
standard_error <- function(x) {sd(x) / sqrt(length(x))}
################################################################
meancol <- function(colorder,colscores,table,nbgroups,nbfile){
	ordered<-table[order(table[,colscores]),]
	leng <- round(dim(ordered)[1] / nbgroups,0)
	colinterest <- ordered[,colscores]
	colinterest2 <- ordered[,colorder]
	meanCOV<-c()
	SpacingMean<-c()
	SpacingRaw<-c()
	nbmax<-c()
	Group<-c()
	i<-floor(leng/2)+1
	j=1
	continue=TRUE
	while (continue){
		deb=i-floor(leng/2)
		fin=i+floor(leng/2)-1
		if (fin > length(colinterest)){
			fin=length(colinterest)
			continue=FALSE
		}
		if (length(colinterest2[deb:fin])>=(0.74*leng)){
		meanCOV[j]=mean(colinterest[deb:fin])
#		stde[j]=sd(colinterest[deb:fin]) # standard deviation
#		stde[j]=standard_error(colinterest[deb:fin]) # standard error
		SpacingMean[j]=(mean(colinterest2[deb:fin])*100)
		SpacingRaw[j]=(sum(colinterest2[deb:fin]))
#		nbmax[j]=length(colinterest2[deb:fin])
		Group[j]=j
#		if(Spacing[j]==0){
#			j=j-1
#		}
		i=i+leng
		j=j+1
		}#else{
#		print(length(colinterest2[deb:fin]))
#		}
	}
#	if(Spacing[j-1]==0){
#		Spacing=head(Spacing,-1)
#		meanCOV=head(meanCOV,-1)
#	}
	df <- data.frame(meanCOV,SpacingMean,Group,SpacingRaw)
	colnames(df)<-c(paste("meanCOV",nbfile,sep=""),paste("SpacingMean",nbfile,sep=""),paste("Group",nbfile,sep=""),paste("SpacingRaw",nbfile,sep=""))
#	df$stde[df$stde>df$meanCOV]=df$meanCOV-0.01
	return(df)
}


ordered<-data.frame(c(1:10))
for (i in 1:length(myfilelist[[1]])){
	tableName=myfilelist[[1]][i]
	table <- read.table(tableName,header=TRUE)
	table$CFR <- table[,name1] /table[,name2]
	ordered<-cbind(ordered,meancol( "Spacing", "CFR", table,10,i))
#	
	
}




#if (sd(ordered$Spacing)==0){
#	print(0.0)
#}else{
corr2<-cor.test(log10(ordered$meanCOV5),ordered$SpacingMean5)
print(as.numeric((corr2$estimate)*(corr2$estimate)))

print(ordered)

#}
#g<-
#for (i in 1:length(myfilelist[[1]])){
#g<-g + 
#	geom_point(aes(x=ordered[,paste("Group",i,sep="")],y=ordered[,paste("SpacingMean",i,sep="")])) 
#}
cbPalette=c('#40A5C7', '#307C95', '#205364')
Cairo(width = 600, height = 500, file=paste("/home/304.3-STRUCTPDEV/MADS/mapping_DAP/pipepline_rep/results/reanalysis/spaced2","/plot_groupsSpacing_all.svg",sep=""), type="svg", pointsize=12, bg = "transparent", canvas = "white", units = "px", dpi = 75)
ggplot(ordered,aes(x=ordered$Group1))+
geom_line(aes(y=ordered$SpacingMean1,color="0.3"),size=1.25)+
geom_point(aes(y=ordered$SpacingMean1,color="0.3"),size=3)+
#geom_line(aes(y=ordered$SpacingMean2,color="0.35"))+
#geom_point(aes(y=ordered$SpacingMean2,color="0.35"))+
geom_line(aes(y=ordered$SpacingMean3,color="0.4"),size=1.25)+
geom_point(aes(y=ordered$SpacingMean3,color="0.4"),size=3)+
#geom_line(aes(y=ordered$SpacingMean4,color="0.45"))+
#geom_point(aes(y=ordered$SpacingMean4,color="0.45"))+
geom_line(aes(y=ordered$SpacingMean5,color="0.5"),size=1.25)+
geom_point(aes(y=ordered$SpacingMean5,color="0.5"),size=3)+
ylim(0,65)+
theme_classic(base_size=15)+ scale_color_manual(values=cbPalette)+
ylab("percentages of sequences with prefered spacing")+
xlab(paste("Groups of sequences ordered by CFR log(","SEP3","\U0394","tet-AG","/SEP3-AG)"))
#ggplot(ordered, aes(x=ordered$Group, y=ordered$SpacingMean)) + 
#geom_line() + geom_point() + geom_bar(aes(x=ordered) + ylim(0,100)
#scale_color_gradientn(colours = rainbow(5))
##geom_histogram(fill="white",binwidth=0.05) + 
#scale_x_log10()
##xlab(paste("ratio Cov ",name1," / Cov ",name2))
dev.off()

