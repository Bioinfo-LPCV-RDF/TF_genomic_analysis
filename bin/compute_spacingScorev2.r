rm(list=ls())
library(ggplot2 ,verbose=FALSE ,warn.conflicts=FALSE, quietly = TRUE)
library(Cairo ,verbose=FALSE ,warn.conflicts=FALSE, quietly = TRUE)
library(latex2exp ,verbose=FALSE ,warn.conflicts=FALSE, quietly = TRUE)

library("optparse")
option_list = list(
  make_option(c("-f", "--files"), type="character", default=NULL, 
              help="tables files ", metavar="character"),
  make_option(c("-n", "--names"), type="character", default=NULL, 
              help=" ", metavar="character"),
  make_option(c("-d", "--out_dir"), type="character", default=NULL, 
              help=" ", metavar="character"),
  make_option(c("-c", "--colors"), type="character", default=NULL, 
              help=" ", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

myfilelist=strsplit(opt$files , ",");
name1=strsplit(opt$names, ",")[[1]][1]
name2=strsplit(opt$names, ",")[[1]][2]
dirout=opt$out_dir

standard_error <- function(x) {sd(x) / sqrt(length(x))}

meancol <- function(colorder,colscores,table,nbgroups,nbfile,Name){
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
		Group[j]=j
		i=i+leng
		j=j+1
		}
	}
	df <- data.frame(meanCOV,SpacingMean,Group,SpacingRaw)
	
	corr2<-cor.test(log10(df$meanCOV),df$SpacingMean)
	print(paste("R2 for decile ",nbfile,": ",as.numeric((corr2$estimate)*(corr2$estimate)),sep=""))
	df$file<-Name
	colnames(df)<-c("meanCOV","SpacingMean","Decile","SpacingRaw","file")
	return(df)
}
################################################################
ordered <- data.frame(meanCOV=double(), SpacingMean=double(), Decile=character(), SpacingRaw=double(), file=character(), stringsAsFactors=FALSE) 


for (i in 1:length(myfilelist[[1]])){
	tableName=myfilelist[[1]][i]
	table <- read.table(tableName,row.names = NULL,header=TRUE)
	Name=strsplit(myfilelist[[1]][i],"/")[[1]][length(strsplit(myfilelist[[1]][i],"/")[[1]])]
	table$CFR <- table[,name1] /table[,name2]
	ordered<-rbind(ordered,meancol( "Spacing", "CFR", table,10,i,Name))
#	
	
}

colors=strsplit(opt$colors, ",")[[1]]
print(colors)
cbPalette=colors#c('#40A5C7', '#307C95', '#205364')
print(ordered)

Cairo(width = 650, height = 500, file=paste(dirout,"/plot_groupsSpacing_all.svg",sep=""), type="svg", pointsize=12, bg = "transparent", canvas = "white", units = "px", dpi = 75)

ggplot(ordered,aes(x=ordered$Decile))+
geom_line(aes(y=ordered$SpacingMean,group=factor(ordered$file),color=factor(ordered$file)),size=1.25)+
geom_point(aes(y=ordered$SpacingMean,group=factor(ordered$file),color=factor(ordered$file)),size=3)+
ylim(0,70)+ # TODO change MAX
theme_classic(base_size=15)+ scale_color_manual(values=cbPalette)+
ylab("percentages of sequences with prefered spacing")+
xlab(paste("Groups of sequences ordered by CFR log(",name1,"/",name2,")",sep=""))

dev.off()
