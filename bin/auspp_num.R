#!/usr/bin/env Rscript
# Version: 1.0
# Modified date: 2018.03.15
Args <- commandArgs()
#Args[1];Args[2];Args[3];Args[4];Args[5];Args[6];
#Args[1];Args[2];Args[3];Args[4];Args[5];Args[6];
#Args[1]= "/usr/lib/R/bin/exec/R"
#Args[2]= "--slave"
#Args[3]= "--no-restore"
#Args[4]= "--file=./lenDist.R"      your program name
#Args[5]= "--args" or NA
#Args[6]= your first parameter start from 6
#Args[6]
if (length(Args) > 6) {
	filename=Args[7]
} else {
	filename="auspp"
}
x=length(strsplit(Args[6],";")[[1]])
#strsplit(Args[6],";")[[1]]
#x
y=length( strsplit(strsplit(Args[6],";")[[1]][1],",")[[1]] )
#y
k=matrix(nrow=x-1,ncol=y)
colnames(k)=strsplit(strsplit(Args[6],";")[[1]][1],",")[[1]]
#k
for(j in 2:x){
	k[j-1,]=strsplit(strsplit(Args[6],";")[[1]][j],",")[[1]]
}
#k

#x=c(as.numeric(strsplit(Args[6],"\n")[[1]]))
#x=c(Args[6])
#barplot(k[,2:y],beside=TRUE)
pdf(file=paste(filename,".auspp.info.pdf",sep=""))
if (y > 2) {
	par(mfrow=c(y-1,1),plt=c(0.2,0.9,0.2,0.9))
	z=matrix(nrow=x-1,ncol=y-1)
	for(j in 2:y) {
		z[,j-1]=as.numeric(k[,j])
	}
#	z=as.matrix(k[,-1])
	rownames(z)=k[,1]
	colnames(z)=colnames(k)[2:y]
	#barplot(z[,1],beside=TRUE,ylim=c(0,round(max(z)/10+0.5)*10) )
	for(j in 1:(y-1)){#2:(y-1)){
		barplot(as.matrix(z[,j]),beside=TRUE,ylim=c(0,round(max(z)*1.05)),horiz=F,names.arg=rownames(z),las=2)#,axes=T)
		title(main = colnames(z)[j], ylab = "Raw read number",mgp=c(4,1,0)) #xlab = "Length (nt)", 
	}
} else {
	par(mfrow=c(y-1,1),plt=c(0.2,0.9,0.2,0.9)) #par(mar=c(8,8,0.5,0.5))
	z=matrix(nrow=x-1,ncol=y-1)
	z[,1]=as.numeric(k[,2])
#	z=as.matrix(k[,-1])
	rownames(z)=k[,1]
	colnames(z)=colnames(k)[2]
#	z
	#barplot(z[,1],beside=TRUE,ylim=c(0,round(max(z)/10+0.5)*10) )
	for(j in 1:(y-1)){#2:(y-1)){
		barplot(as.matrix(z[,j]),beside=TRUE,ylim=c(0,round(max(z)*1.05)),horiz=F,names.arg=rownames(z),las=2)#,axes=T)
		title(main = colnames(z)[j], ylab = "Raw read number",mgp=c(4,1,0))
	}
}

dev.off()
#aplot<-function(x){
#x<-scan()
#mydata<-data.frame(age=numeric(0),gender=character(0),weight=numeric(0))
#mydata<-edit(mydata)
#plot(x)
#}#
