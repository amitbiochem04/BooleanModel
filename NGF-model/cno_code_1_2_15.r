
#### CNO model simulation 
library(CellNOptR)
library(CNORdt)
library(CNORfeeder)
rm(list=ls())
# Read prior knowledge network
pknmodel=readSIF("PKN-network.sif")
###plot the pknmodel to look how the interaction between them.  
pdf(file="PKN.pdf",width=21/1.54,height=10/1.54)
plotModel(pknmodel)
dev.off()
### Read the Normalized data (between 0 and 1). I took log2fold change value of qpcer data and applied logistic function to have value between 0 and 1
data1 = CNOlist("final_data_model_5_10.csv")

###Now arrange the data with list format because it required for further analysis
a<-colnames(getCues(data1))
b<-colnames(getStimuli(data1))
c<-colnames(getInhibitors(data1))
d<-colnames(getSignals(data1)[[1]])
e<-c(0,1,2,4,6,8,12,24)
e<-as.integer(e)
###matrixformat
g<-getCues(data1)
h<-getInhibitors(data1)
i<-getStimuli(data1)
###all the in time points
j<-getSignals(data1)[[1]]
k<-getSignals(data1)[[2]]
l<-getSignals(data1)[[3]]
m<-getSignals(data1)[[4]]
n<-getSignals(data1)[[5]]
o<-getSignals(data1)[[6]]
p<-getSignals(data1)[[7]]
q<-getSignals(data1)[[8]]
#r<-getSignals(data1)[[9]]
###list the timepoint
data2<-list("t0"=j,"t1"=k,"t2"=l,"t4"=m,"t6"=n,"t8"=o,"t12"=p,"t24"=q)
cnolist2<-list("namesCues"=a,"namesStimuli"=b,"namesInhibitors"=c,"namesSignals"=d,"timeSignals"=e,"valueCues"=g,"valueInhibitors"=h,"valueStimuli"=i,"valueSignals"=data2)
colnames(cnolist2$valueCues)<-NULL
#colnames(cnolist2$valueSignals[[9]])<-NULL
colnames(cnolist2$valueSignals[[8]])<-NULL
colnames(cnolist2$valueSignals[[7]])<-NULL
colnames(cnolist2$valueSignals[[6]])<-NULL
colnames(cnolist2$valueSignals[[5]])<-NULL
colnames(cnolist2$valueSignals[[4]])<-NULL
colnames(cnolist2$valueSignals[[3]])<-NULL
colnames(cnolist2$valueSignals[[2]])<-NULL
colnames(cnolist2$valueSignals[[1]])<-NULL
save(cnolist2,file="cnolist_ngf_last_dat.mat")

###plot experimental data
quartz()
pdf(file="plotdata.pdf",width=21/1.54,height=10/1.54)
plotCNOlist(cnolist2)
dev.off()
###### pknmodel pre-processing, model was compressed and you can find the index of the reactions of the model
model <- preprocessing(cnolist2, pknmodel,expansion=TRUE,compression=TRUE,verbose=FALSE,cutNONC=TRUE)
indexOrig<-indexFinder(cnolist2, model,verbose=FALSE)
res<-residualError(cnolist2)
init<-rep(1,length(model$reacID))
##### Here you can fix some of the edge to 1 to check which edge you want to fix ""call comand model$reacID
priorBitString = rep(NA, length(model$reacID))
####Here I only show two edge I fixed but there are 6-7 edge I fixed 
priorBitString[1] = 1
priorBitString[19] = 1
priorBitString[20] = 1
priorBitString[22] = 1
priorBitString[36] = 1
priorBitString[38] = 1
priorBitString[39] = 1
priorBitString[77] = 1
priorBitString[63] = 1
fields4Sim <- prep4sim(model)

#####model optimization. I wanted to run 100 time "gaBinaryDT" function because each and every time if I run "gaBinaryDT" then I get different best string (edges). If I run 100 times and if i get more than 90% same edges then it make sense. 
## The choice of repeating the analysis multiple times stems both from the stochastic nature of the optimization procedure and also from the fact that the training data were not sufficient to fully constrain the model (Saez‐Rodriguez et al, 2009).
N =100
opt = list()
optres=list()
for (i in 1:N){
	opt1 <- gaBinaryDT(CNOlist=cnolist2, model=model, initBstring=init,verbose=FALSE, boolUpdates=10, maxTime=100, lowerB=0.5, upperB=10,popSize=50,relTol=0.01,priorBitString=priorBitString)
    opt[[i]] = opt1
	simResults <- simulatorDT(CNOlist=cnolist2,model=model,simList=fields4Sim,indices=indexOrig,boolUpdates=10)
	simResults = convert2array(simResults, dim(cnolist2$valueSignals[[1]])[1],length(model$namesSpecies), 10)
	optimRes <- getFitDT(simResults=simResults,CNOlist=cnolist2,model=model,indexList=indexOrig,boolUpdates=10,lowerB=0.5,upperB=10,nInTot=length(which(model$interMat == -1)))
	optres[[i]]=optimRes
}
save(opt,file="opt_res_with_com.mat")
save(optimres,file="opt_residual_with_com.mat")

####Extracting the best string  
####extract the best string (best edges) from all 100 runs. Here, I have 121 edges in my compressed model. The main aim of 100 times run is to have stochastic nature of the network.

totalbeststring<-rep(0,114)
for (i in opt)
{
	bs<-i$bString
	totalbeststring<-totalbeststring+bs
}
cutofffct<-function(x,cutoff)
{
	if (x > cutoff) return(1)
	return(0)
}
#### one can uses different cutoff to have a right topology of the network but it would be nice to have 90 %plot the data
avgbs<-lapply(totalbeststring,cutofffct,cutoff=70)
avgbs <- unlist(avgbs)

#####plot the bestmodel
plotModel(model, cnolist2, bString=avgbs)
#plotModel(model, cnolist2, bString=opt1$bString)
pdf(file="plot_fit.pdf",width=21/1.54,height=10/1.54)
plotCNOlist(cnolist2)
dev.off()


########################################## plot the optimized model
#plot the model back, it will plot only the edges which obtained after optimization.
#plotModel(model, cnolist2, bString=opt1$bString)
#cutAndPlotResultsDT(model,CNOlist=cnolist2,bString=opt1$bString,plotPDF=FALSE, boolUpdates=15,lowerB=0.5,upperB=10,)
#plotModel(model, cnolist2, bString=avgbs)
#cutAndPlotResultsDT(model,CNOlist=cnolist2,bString=avgbs,plotPDF=FALSE, boolUpdates=15,lowerB=0.5,upperB=10,plotParams=list(cex=0.8, cmap_scale=0.5, margin=0.2,maxrow=25, margin=0.1, width=20, height=20))
#l<-c(1,1, 0, 1 ,0 ,0 ,0 ,1 ,1, 0 ,0 ,0 ,1 ,0 ,0 ,0 ,1 ,1 ,1 ,1 ,1 ,0,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0, 0 ,0 ,0 ,0 ,0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
#0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,0,0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0)
#plotModel(model, cnolist2, bString=l)
#cutAndPlotResultsDT(model,CNOlist=cnolist2,bString=opt1$bString,plotPDF=FALSE, boolUpdates=10,lowerB=0.5,upperB=10,plotParams=list(cex=0.8, width=20, height=20))
#ggplot(data=a,aes(reorder(x= reaction,app)),y=app+ coord_flip() +geom_bar(colour="black", stat="identity",width=.5)
#p<-qplot(data=meltdf, x = variable, y=value, geom="bar", stat = "identity" )+ facet_wrap( "symbol" ))
#ggplot(data=meltdf,aes(x=reaction,y=value)+ coord_flip() +geom_bar(colour="black", stat="identity",width=.3)

#################  
# Extracting best score which is the minimum value
#best <- c()
#as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}
#for (i in 1:N){
#b <- as.data.frame(opt[[i]]$results)
#	best[i] <- min(as.numeric.factor(b$Best_score))
#}



#############histogram of avergae score
#h<-hist(best, breaks=100,col="Lightblue", xlab="Avrage Bestscore",ylim=c(0,25),xlim=c(0.012,0.018))
#xfit<-seq(min(best),max(best),length=100)
#yfit<-dnorm(xfit,mean=mean(best),sd=sd(best))
#yfit <- yfit*diff(h$mids[1:2])*length(best)
#lines(xfit, yfit, col="blue",lwd=2)

####### Here I looked for percentage of best string in 100 runs
#N=100
#matr <- matrix(nrow=106,ncol=N)
#for (i in 1:N){
#	matr[,i] <- opt[[i]]$bString
#}
#percent <- rowSums(matr)*100/N
#save(matr,file="matr.mat")
###
a<-read.csv("test_corelation.csv",header=T,sep=",",)
m<-as.matrix(a[-1,])
row.names(m)<-a[,1]
m1<-cor(t(m),method="spearman")
m1[is.na(m1)] <- 0

#####Plot corelation matrix
row.names(matr)<-model$reacID
m <- cor(t(matr),method="spearman")
m[is.na(m)] <- 0
CairoPDF(file="corelationmatrix",width=35/1.54,height=35/1.54)
pheatmap(m,cluster_rows=FALSE,cluster_col=FALSE, fontsize =7,col=c("black","navajowhite","snow4","lightskyblue"))
dev.off()

options(scipen=999)
d<-read.delim("ECvalue.txt",header=T,sep="\t")
hill<- function(x){
    
    for (i in d$EC50value){
        EC50Data<-i
    }
    
    HillCoef=3
    x^HillCoef/(EC50Data^HillCoef+x^HillCoef)
}


hill<- function(x){
    
    x^HillCoef/(EC50Data^HillCoef+x^HillCoef)
}

EC50Data<-0.9
HillCoef=3



