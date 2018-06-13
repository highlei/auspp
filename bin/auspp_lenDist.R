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
	k[j-1,]=as.numeric(strsplit(strsplit(Args[6],";")[[1]][j],",")[[1]])
}
#k
#x=c(as.numeric(strsplit(Args[6],"\n")[[1]]))
#x=c(Args[6])
#barplot(k[,2:y],beside=TRUE)
pdf(file=paste(filename,".auspp.lenDist.pdf",sep=""))
if (y > 2) {
	par(mfrow=c(y-1,1))
	z=k[,-1]
	rownames(z)=k[,1]
	for(j in 1:(y-1)){
		z[,j]=z[,j]/sum(z[,j])
	}
	z=z*100
	#z
	#barplot(z[,1],beside=TRUE,ylim=c(0,round(max(z)/10+0.5)*10) )
	for(j in 1:(y-1)){#2:(y-1)){
		barplot(z[,j],beside=TRUE,ylim=c(0,round(max(z)/10+0.5)*10))#,axes=F)
		title(main = colnames(z)[j], xlab = "Length (nt)", ylab = "Percentage %")
	}
} else {
#	par(mfrow=c(y-1,1))
	z=as.matrix(k[,-1])
	rownames(z)=k[,1]
	colnames(z)=colnames(k)[2]
	for(j in 1:(y-1)){
		z[,j]=z[,j]/sum(z[,j])
	}
	z=z*100
	#z
	#barplot(z[,1],beside=TRUE,ylim=c(0,round(max(z)/10+0.5)*10) )
	for(j in 1:(y-1)){#2:(y-1)){
		barplot(z[,j],beside=TRUE,ylim=c(0,round(max(z)/10+0.5)*10))#,axes=F)
		title(main = colnames(z)[j], xlab = "Length (nt)", ylab = "Percentage %")
	}
}

dev.off()
#aplot<-function(x){
#x<-scan()
#mydata<-data.frame(age=numeric(0),gender=character(0),weight=numeric(0))
#mydata<-edit(mydata)
#plot(x)
#}#
