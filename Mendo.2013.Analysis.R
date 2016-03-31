#clear all objects from R console
  rm(list=ls())
  graphics.off()
  pardefault <- par(no.readonly = T)

#set working directory
  setwd("/Users/mmcmunn/Desktop/GitHub/Mendocino")

#load ggplot2, vegan, and reshape
  library("ggplot2")
  library("vegan")
  library("reshape")
  library("plyr")
  library("treemap")
  library("mvabund")
  # graphic or grid package unavail, problem?

#read in data
  d<-read.csv("Mendo.July.2013.Night.Day.family.11414.csv", header=T)
  clim<-read.csv("mendocino.climate.var.all.summ.csv",header=T)
  trophic<-read.csv("trophic.family.assign.csv", header=T,na.strings="")
  
#error bar function
  error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
        if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
        stop("vectors must be same length")
        arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)}
  
#define standard error that removes NA's
  se<-function(x) sqrt( var(x[ !is.na(x) ] ) / length(x[ !is.na(x) ] ))

#create new vector for unique families
  d$ord.fam<-paste(d$Order,d$Family,sep=".")

#write unique families to csv
#ord.fam<-sort(unique(d$ord.fam))
#write.csv(ord.fam,file="ord.fam.csv")

#all collembolans stripped of family ID for analyses (unreliable ID to family)
  d$ord.fam <- ifelse(d$Order == "Collembola" , "Collembola." , d$ord.fam)

  
#create new vector for insect volume
  d$volume <- d$Length * ( pi * ( .5 * d$Width )^2)

  
#clean up time variables, change to R time objects
  d$dateTimeEnd <- strptime(paste("2013" , d$End.Date, d$End.Time, sep = "-"), format = "%Y-%d-%b-%H:%M")
  clim$dateTimeEnd <- strptime(clim$end.time , format = "%m/%d/%y %H:%M")
  
#trophic assignments
#convert to characters, rather than factors, which R turns into factor-levels
  #do any ord.fams fail to match?
  d[ (which(is.na(match(d$ord.fam,trophic$ord.fam)))) , ]  
  
  trophic$trophic.assign  <- with(trophic,  ifelse(is.na(trophic.assign), as.character(trophic.position), as.character(trophic.assign) ))

  d$trophic.assign <- trophic$trophic.assign[match(d$ord.fam,trophic$ord.fam)]

#sampleID
    d$sampleID <- with (d, paste(Malaise.Pit ,as.numeric(as.factor(as.character(d$dateTimeEnd))), sep = ""))
    d$timeChar <- as.character(d$dateTimeEnd)
        ##############################################
#make community matrices
#average # of arthropods in each sample
  mean(table(d$sampleID))
#abundance within samples by family
  abund.FxS <- table(d$timeChar, d$ord.fam)

#abundance within samples by trophic position
  abund.TxS <-table(d$timeChar,d$trophic.assign)

#biomass within samples by family
  vol.FxS <- with(d, tapply(volume, list(ord.fam, timeChar), sum))
  vol.FxS[is.na(vol.FxS)] <- 0

#biomass within samples by trophic position
  vol.TxS <- with(d, tapply(volume, list(trophic.assign, timeChar), sum))
  vol.TxS[is.na(vol.TxS)] <- 0

###############################
#abundance based analyses
###############################
#multivariate abundance models using mvabund package
#####by family
m.time <- manyglm(abund.FxS ~ clim$end.time, family = "negative.binomial")
m.time.sig <- anova(m.time, p.uni = "adjusted")
m.time.sig
#time of day is significant (.001) as a factor determining community composition

#check residuals vs fitted
plot(m.time)
#check mean variance relationship
meanvar.plot(comm.by.fam[,] ~ time.of.day)

#how many taxa varied in abundance significantly by time of day
sum(m.time.sig$uni.p[2,]<=.05)
which(m.time.sig$uni.p[2,]<=.05)

#which taxa were most abundant
sort(colSums(comm.by.fam),decreasing=TRUE)[1:15]
#some notable taxa that were abundant, but did not vary by time of day - ants, cecidomyiids

#a closer look at taxa driving changes over the course of the day
#get top 15 likelihood statistics

top.fifteen.uni <- names(sort(-m.time$two.loglike, decreasing=TRUE)[1:15])
top.fifteen.lik <- sort(-m.time$two.loglike, decreasing=TRUE)[1:15]
top.fifteen.values <- top.fifteen.uni


d$sample.type.time <- paste(d$time.method,d$End.Date,sep=".")


temp.counts <- with(d, aggregate(Malaise.Pit, by=list(ord.fam,sample.type.time, end.hour),length))
relev.counts <- temp.counts[which(!is.na(match(temp.counts$Group.1 , top.fifteen.uni))) , ]
means.to.plot <- with(relev.counts, aggregate(x , by=list(Group.1,Group.3), mean))
colnames(means.to.plot) <- c("ord.fam", "end.hour", "mean.abund")
means.to.plot<-means.to.plot [order(match(means.to.plot$ord.fam, top.fifteen.uni)),]
to.plot <- cbind( top.fifteen.lik[match(means.to.plot$ord.fam, top.fifteen.uni)] , means.to.plot)

plot(to.plot[,1], to.plot[,4], pch=16, xlab="likelihood statistic", ylab="mean abundance")


#####by trophic
m.time.tro <- manyglm(comm.by.tro ~ time.of.day, family = "negative.binomial")
m.time.tro.sig <- anova(m.time.tro, p.uni = "adjusted")

#how many trophic levels varied in abundance significantly by time of day
sum(m.time.tro.sig$uni.p[2,]<=.05)
which(m.time.tro.sig$uni.p[2,]<=.05)

########################################
#plot NMDS of community family matrix with Malaise and Pitfall bulked
set.seed(1231231)
ord.all= metaMDS(comm=as.matrix(comm.by.fam),distance="bray")

cols<-c(rep("red",5),rep("yellow",5),rep("green",5),rep("black",5),rep("purple",5),rep("gray50",5))

xy.coord<-scores(ord.all)
rownames(xy.coord)
hours <- as.numeric(unlist(lapply(strsplit(rownames(xy.coord),split="\\."), function(x) x[1])))
days <- as.factor(unlist(lapply(strsplit(rownames(xy.coord),split="\\."), function(x) x[2])))
type <- ifelse(hours==22|hours==10, "crepuscular", ifelse(hours==2|hours==6, "night", "day"))

xy.ordering <- data.frame(hours = as.numeric(hours),days = as.numeric(days),type)
xy.coord <- xy.coord[order(xy.ordering[,2],xy.ordering[,1] ) , ]
cols <- cols[order(xy.ordering[,2],xy.ordering[,1] )]

typeOrder <- xy.ordering[order(xy.ordering[,2],xy.ordering[,1] ),"type"]


#for gif
setwd("/Users/mmcmunn/Desktop/Yang_Lab/Insect Temporal Diversity/ordination animation")
name<-paste("ord.plot",31,"png",sep=".")
png(name, width = 600, height = 600)
par(oma = c(0,1,0,0))
plot(x=xy.coord[,1], y=xy.coord[,2],col=cols,pch=16,cex=2.5,cex.lab=1.5,cex.axis=1.5,
     xlim = c(min(xy.coord[,1])+min(xy.coord[,1]*.3) , max(xy.coord[,1]+max(xy.coord[,1]*.3))),
     ylim = c(min(xy.coord[,2])+min(xy.coord[,2]*.3) , max(xy.coord[,2]+max(xy.coord[,2]*.4))),
     xlab = "NMDS 1", ylab = "NMDS 2"
)
legend("topright",legend=c("10pm-2am","2am-6am","6am-10am","10am-2pm","2pm-6pm","6pm-10pm"),pch=rep(16,6),col=c("black","gray50","red","yellow","green","purple"),cex=1.5)

crep.points <- xy.coord[which(typeOrder=="crepuscular"),]
crep.bound <- crep.points[chull(crep.points),]
polygon(crep.bound)
text(mean(crep.points[,1]), mean(crep.points[,2])-.5, "crepuscular", cex = 1.5)


day.points <- xy.coord[which(typeOrder=="day"),]
day.bound <- day.points[chull(day.points),]
polygon(day.bound)
text(mean(day.points[,1]), mean(day.points[,2])-.5, "diurnal", cex = 1.5)


night.points <- xy.coord[which(typeOrder=="night"),]
night.bound <- night.points[chull(night.points),]
polygon(night.bound)
text(mean(night.points[,1]), mean(night.points[,2]-.4), "nocturnal", cex = 1.5)
dev.off()

for(i in 1:30){
	name<-paste("ord.plot",i,"png",sep=".")
	png(name, width = 600, height = 600)
	plot(xy.coord,type="n",col=cols,pch=16,cex=1.5,cex.lab=1.5,cex.axis=1.5,
	xlim = c(min(xy.coord[,1])+min(xy.coord[,1]*.3) , max(xy.coord[,1]+max(xy.coord[,1]*.3))),
	ylim = c(min(xy.coord[,2])+min(xy.coord[,2]*.3) , max(xy.coord[,2]+max(xy.coord[,2]*.3))),
	xlab = "NMDS 1", ylab = "NMDS 2", main = "Family-level abundance by collection time"
	)
points(x=xy.coord[1:i,1], y=xy.coord[1:i,2],col=cols[1:i],pch=16,cex=2.5) 
legend("topright",legend=c("10pm-2am","2am-6am","6am-10am","10am-2pm","2pm-6pm","6pm-10pm"),pch=rep(16,6),col=c("black","gray50","red","yellow","green","purple"),cex=1.5)

for (h in 1:i ){
	segments(x0 = xy.coord[ifelse(h >1, h-1, h) , 1 ], y0= xy.coord[ifelse(h >1, h-1, h) , 2], x1 = xy.coord[h, 1], y1 = xy.coord[h, 2]   )	
}
arrows(x0 = xy.coord[ifelse(i >1, i-1, i) , 1 ], y0= xy.coord[ifelse(i >1, i-1, i) , 2], x1 = xy.coord[i, 1], y1 = xy.coord[i, 2]   , lwd=6 )
dev.off()
}


setwd("/Users/mmcmunn/Desktop/Yang_Lab/Insect Temporal Diversity")



#PERMANOVA for time of day, second for abiotic variables
adonis(comm.by.fam ~ as.factor(sample.info[,2]))

sample.info<-cbind(rownames(comm.by.tro),c(rep(10,5),rep(14,5),rep(18,5),rep(2,5),rep(22,5),rep(6,5)))

adonis(comm.by.fam ~ sample.info$light.shade+sample.info$mean.T+sample.info$wind.max)


#plot NMDS of community trophic matrix with Malaise and Pitfall bulked
ord.all= metaMDS(comm=comm.by.tro)
plot(ord.all,type="t",display="sites")
points(ord.all,col=cols,pch=1,cex=4,lwd=2)
adonis(comm.by.tro ~ as.factor(sample.info[,2]))
legend(.3,.6,legend=c("10pm-2am","2am-6am","6am-10am","10am-2pm","2pm-6pm","6pm-10pm"),pch=rep(16,6),col=c("black","gray50","red","yellow","green","purple"))
mtext("mendocino insect trophic position by time of day",3,line=1,cex=1.5)

#repeat NMDS with Malaise and Pitfall separate
comm.by.fam.type<-table(d$sample.time, d$ord.fam,d$Malaise.Pit)

pitfall <- comm.by.fam.type[,,"P"]
malaise <- comm.by.fam.type[,,"M"]

rownames(pitfall) <- paste(rownames(pitfall), "P", sep = ".")
rownames(malaise) <- paste(rownames(malaise), "M", sep = ".")
comm.by.fam.type <- rbind(pitfall, malaise)

set.seed(1231231)
ord.all= metaMDS(comm=comm.by.fam.type)

xy.coord<-scores(ord.all)
#pitfall are columns 1 and 2. malaise are columns 3 and 4
xy.coord <- cbind(xy.coord[1:30,], xy.coord[31:60,])

hours <- as.numeric(unlist(lapply(strsplit(rownames(xy.coord),split="\\."), function(x) x[1])))
days <- as.factor(unlist(lapply(strsplit(rownames(xy.coord),split="\\."), function(x) x[2])))
type <- ifelse(hours==22|hours==10, "crepuscular", ifelse(hours==2|hours==6, "night", "day"))

xy.ordering <- data.frame(hours = as.numeric(hours),days = as.numeric(days),type)
xy.coord <- xy.coord[order(xy.ordering[,2],xy.ordering[,1] ) , ]
cols<-c(rep("red",5),rep("yellow",5),rep("green",5),rep("black",5),rep("purple",5),rep("gray50",5))
cols <- cols[order(xy.ordering[,2],xy.ordering[,1] )]

typeOrder <- xy.ordering[order(xy.ordering[,2],xy.ordering[,1] ),"type"]


#for gif
setwd("/Users/mmcmunn/Desktop/Yang_Lab/Insect Temporal Diversity/ordination animation 2")
name<-paste("ord.plot",31,"png",sep=".")
png(name, width = 600, height = 600)
plot(x=c(xy.coord[,1],xy.coord[,3]), y=c(xy.coord[,2],xy.coord[,4]),
     col=cols,cex=2.5,cex.lab=1.5,cex.axis=1.5, pch = c(replicate(30, 16), replicate(30, 17)),
     xlab = "NMDS 1", ylab = "NMDS 2", main = "Family-level abundance by collection time and trap type")

legend("bottomleft",ncol = 3,legend=c("10pm-2am","2am-6am","6am-10am","10am-2pm","2pm-6pm","6pm-10pm", "Pitfall", "Malaise"),
       pch=c(rep(16,6), 16, 17),col=c("black","gray50","red","yellow","green","purple", "black", "black"),cex=1)
dev.off()
for(i in 1:30){
  name<-paste("ord.plot",i,"png",sep=".")
  png(name, width = 600, height = 600)
  plot(x=c(xy.coord[,1],xy.coord[,3]), y=c(xy.coord[,2],xy.coord[,4]),type = "n",
       col=cols,cex=2.5,cex.lab=1.5,cex.axis=1.5, pch = c(replicate(30, 16), replicate(30, 17)),
       xlab = "NMDS 1", ylab = "NMDS 2", main = "Family-level abundance by collection time and trap type")
  
  legend("bottomleft",ncol = 3,legend=c("10pm-2am","2am-6am","6am-10am","10am-2pm","2pm-6pm","6pm-10pm", "Pitfall", "Malaise"),
         pch=c(rep(16,6), 16, 17),col=c("black","gray50","red","yellow","green","purple", "black", "black"),cex=1)
  
  
  points(x=xy.coord[1:i,1], y=xy.coord[1:i,2],col=cols[1:i],pch=16,cex=2.5)  
  points(x=xy.coord[1:i,3], y=xy.coord[1:i,4],col=cols[1:i],pch=17,cex=2.5)   
  
  for (h in 1:i ){
    segments(x0 = xy.coord[ifelse(h >1, h-1, h) , 1 ], y0= xy.coord[ifelse(h >1, h-1, h) , 2], x1 = xy.coord[h, 1], y1 = xy.coord[h, 2]   )	
    segments(x0 = xy.coord[ifelse(h >1, h-1, h) , 3 ], y0= xy.coord[ifelse(h >1, h-1, h) , 4], x1 = xy.coord[h, 3], y1 = xy.coord[h, 4]   )  
    
    
  }
  arrows(x0 = xy.coord[ifelse(i >1, i-1, i) , 1 ], y0= xy.coord[ifelse(i >1, i-1, i) , 2], x1 = xy.coord[i, 1], y1 = xy.coord[i, 2]   , lwd=6 )
  arrows(x0 = xy.coord[ifelse(i >1, i-1, i) , 3 ], y0= xy.coord[ifelse(i >1, i-1, i) , 4], x1 = xy.coord[i, 3], y1 = xy.coord[i, 4]   , lwd=6 )
  
  dev.off()
}


setwd("/Users/mmcmunn/Desktop/Yang_Lab/Insect Temporal Diversity")











plot(ord.all,type="t",display="sites")
cols<-c(rep("red",5),rep("yellow",5),rep("green",5),rep("black",5),rep("purple",5),rep("gray50",5))
points(ord.all,col=c(cols),pch=1,cex=3,lwd=2)
legend(-1.05,.95,legend=c("10pm-2am","2am-6am","6am-10am","10am-2pm","2pm-6pm","6pm-10pm"),pch=rep(16,6),col=c("black","gray50","red","yellow","green","purple"))
mtext("mendocino insect families by time of day - MP",3,line=1,cex=1.5)




################################
#biomass based analysis
################################
#plot NMDS of community family matrix with Malaise and Pitfall bulked

ord.all= metaMDS(comm=vol.by.fam,distance="bray")

plot(ord.all,type="t",display="sites")
cols<-c(rep("red",5),rep("yellow",5),rep("green",5),rep("black",5),rep("purple",5),rep("gray50",5))
points(ord.all,col=cols,pch=1,cex=3,lwd=2)
legend(.4,-.4,legend=c("10pm-2am","2am-6am","6am-10am","10am-2pm","2pm-6pm","6pm-10pm"),pch=rep(16,6),col=c("black","gray50","red","yellow","green","purple"))
mtext("mendocino insect family biomass by time of day",3,line=1,cex=1.5)

#PERMANOVA for time of day, second for abiotic variables
adonis(vol.by.fam ~ as.factor(sample.info[,2]))
adonis(vol.by.fam ~ sample.info$light.shade+sample.info$mean.T+sample.info$wind.max)
vol.by.fam[1:31,1]

#plot NMDS of community trophic matrix with Malaise and Pitfall bulked
sample.info<-cbind(rownames(comm.by.tro),c(rep(10,5),rep(14,5),rep(18,5),rep(2,5),rep(22,5),rep(6,5)))
ord.all= metaMDS(comm=comm.by.tro)
plot(ord.all,type="t",display="sites")
points(ord.all,col=cols,pch=1,cex=4,lwd=2)
adonis(comm.by.tro ~ as.factor(sample.info[,2]))
legend(.3,-.3,legend=c("10pm-2am","2am-6am","6am-10am","10am-2pm","2pm-6pm","6pm-10pm"),pch=rep(16,6),col=c("black","gray50","red","yellow","green","purple"))
mtext("mendocino insect trophic position biomass by time of day",3,line=1,cex=1.5)



############################
#sum of squares within family
#############################
#cut out families collected in one sample
comm.by.fam.multiples<-comm.by.fam[,which(colSums(comm.by.fam>0)>1)]

#standardize abundance in each column
comm.by.fam.multiples.st <- sapply(1:ncol(comm.by.fam.multiples), function(x) (comm.by.fam.multiples[,x]-mean(comm.by.fam.multiples[,x]))/
sd(comm.by.fam.multiples[,x]))


#vector of times for each sample	
times<-as.factor(unlist(lapply(strsplit(rownames(comm.by.fam.multiples),split="\\."), function(x) x[1])))

#empty lists for loop
sum.of.sq.fam <-list()
f.stat.fam <- list()

#loop over each family, pull out sum.sq and f-statistic
for(i in 1:ncol(comm.by.fam.multiples.st)){
	
temp <- cbind(comm.by.fam.multiples.st[,i],times)

f.stat.fam[i] <- summary(aov(temp[,1]~temp[,2]))[[1]][,4][1]
sum.of.sq.fam[i] <- summary(aov(temp[,1]~temp[,2]))[[1]][,2][1]
}

#look at loop output
hist(unlist(f.stat.fam))
hist(unlist(sum.of.sq.fam))
cbind(colnames(comm.by.fam.multiples), sum.of.sq.fam)

#############################################
#mean within orders
length.means<-tapply(X=d$Length,INDEX=d$Order,FUN=mean,na.rm=T)
length.means

#se within orders
length.ses<-tapply(X=d$Length,INDEX=d$Order,FUN=se)
length.ses

length.summary<-cbind(length.means,length.ses)


#to order by size
length.to.plot<-length.summary[order(-length.summary[,1]),]
labels<-rownames(length.to.plot)
#the plot
predplot <- barplot(height=length.to.plot[,1], names.arg= c(rownames(length.to.plot)) ,las=2,main="body length by Order" ,ylim = c(0,20),col="grey",xaxt="n",ylab="body length (mm)")
error.bar(predplot,length.to.plot[,1],length.to.plot[,2])
text(cex=.8,x=predplot,y=-1.75,labels,xpd=T,srt=90)

#calculate sizes within each order and time block
sizes<-aggregate(x=d$Length,by=list(d$Order,d$end.hour),FUN=mean)
sizes
#se within order and family
size.se<-aggregate(d$Length,list(d$Order,d$end.hour),se)
size.se

#total sample sizes to check for suspiciously low # samples
aggregate(x=d$Length,by=list(d$End.Date,d$time.method),FUN=length)

#abundances of insects in each of these catagories
abund.temp<-aggregate(d$Length,list(d$Order,d$end.hour,d$End.Date),FUN=length)
abund.mean<-aggregate(abund.temp[,4],list(abund.temp[,1],abund.temp[,2]),FUN=mean)
abund.se<-aggregate(abund.temp[,4],list(abund.temp[,1],abund.temp[,2]),FUN=se)

abundances
#volumes
volumes<-aggregate(x=d$volume,by=list(d$Order,d$end.hour),FUN=mean)
volume.se<-aggregate(x=d$volume,by=list(d$Order,d$end.hour),FUN=se)

#paste together
summary.table<-cbind(sizes,size.se[,3],abund.mean[,3],volumes[,3],volume.se[,3],abund.se[,3])
summary.table

#the data look right, rename columns
colnames(summary.table)<-c("order","end.hour","mean.length","se.length","abundance","volume","volume.se","abund.se")
summary.table

#reorder by order then by time
summary.table<-summary.table[order(summary.table[,1],as.numeric(summary.table[,2])),]


#we want only orders that occur in multiple time steps for comparison of size
#counts of windows of occurence
sum.occur<-tapply(X=summary.table$end.hour,INDEX=summary.table$order,FUN=length)
sum.occur<-subset(sum.occur,sum.occur>1)
length(sum.occur)
ord.names<-rownames(sum.occur)

#open an empty 4x4 matrix of plots of BODY LENGTH
par(mfrow=c(4,4))

#using a loop, fill in each of the 16 plots
for(i in 1:length(sum.occur)){
d.temp<-subset(summary.table,summary.table[,1]==ord.names[i])

predplot <- barplot(height=d.temp[,3], names.arg= d.temp[,2], ylim = c(0,28),col="grey",ylab="mean length (mm)",xlab="collection end time (hours)",main=ord.names[i])
error.bar(predplot,d.temp[,3],d.temp[,4])	
		
}
summary.table

#open an empty 4x4 matrix of plots of BODY VOLUME
par(mfrow=c(4,4))

#using a loop, fill in each of the 12 plots
for(i in 1:length(sum.occur)){
d.temp<-subset(summary.table,summary.table[,1]==ord.names[i])
predplot <- barplot(height=d.temp[,6], names.arg= d.temp[,2],col="grey",ylim=c(0,(max(d.temp[,6])*1.5)),ylab="mean volume (mm^3)",xlab="collection end time (hours)",main=ord.names[i])
error.bar(predplot,d.temp[,6],d.temp[,7])	
		
}



summary.table

#open an empty 3x4 matrix of plots of ABUNDANCE
par(mfrow=c(4,4))

#using a loop, fill in each of the 12 plots
for(i in 1:length(sum.occur)){
d.temp<-subset(summary.table,summary.table[,1]==ord.names[i])
predplot <- barplot(height=d.temp[,5], names.arg= d.temp[,2],col="grey",ylab="mean abundance",xlab="collection end time (hours)",main=ord.names[i])
error.bar(predplot,d.temp[,5],d.temp[,8])	

}



#redo same graphs splitting by day, night, and crepuscular
d$comm.type<-0
d$comm.type<-ifelse(d$end.hour==10|d$end.hour==22,"crepuscular",d$comm.type)
d$comm.type<-ifelse(d$end.hour==14|d$end.hour==18,"day",d$comm.type)
d$comm.type<-ifelse(d$end.hour==2|d$end.hour==6,"night",d$comm.type)

size<-aggregate(d$Length,list(d$comm.type,d$Order),mean)
size.se<-aggregate(d$Length,list(d$comm.type,d$Order),se)
vol<-aggregate(d$volume,list(d$comm.type,d$Order),mean)
vol.se<-aggregate(d$volume,list(d$comm.type,d$Order),se)
abund<-aggregate(d$Length,list(d$comm.type,d$Order),length)
abund.se<-aggregate(d$Length,list(d$comm.type,d$Order),se)

day.night.summ<-cbind(size,size.se[,3],abund[,3],abund.se[,3],vol[,3],vol.se[,3])

##########For ENT SOC

means.comm <- with(d , tapply(Length , comm.type, mean))
ses.comm <- with(d , tapply(Length , comm.type, se))


predplot <- barplot(means.comm ,col=c("grey80","white","grey20"),ylab="",xlab="", ylim=c(0,5.5), border=c(FALSE,TRUE,FALSE),cex.names=1.7,cex.axis=1.5)
error.bar(predplot, means.comm, ses.comm, length=.2,lwd=1)	

abund.comm <- with(d , tapply(Length , comm.type, length))


predplot <- barplot(abund.comm ,col=c("grey80","white","grey20"),ylab="",xlab="", ylim=c(0,2000), border=c(FALSE,TRUE,FALSE),cex.names=1.7,cex.axis=1.5)
error.bar(predplot, means.comm, ses.comm, length=.2,lwd=1)	


#plot all the data

length(unique(day.night.summ[,2]))
ord.names<-unique(day.night.summ[,2])


#night and day BODY LENGTH, labeled with sample size on x axis
par(mfrow=c(5,4))
#using a loop, fill in each of the 17 plots
for(i in 1: length(unique(day.night.summ[,2]))){
	
d.temp<-subset(day.night.summ,day.night.summ[,2]==ord.names[i])
predplot <- barplot(height=d.temp[,3], names.arg= paste(d.temp[,1],"  ","(",d.temp[,5],")",sep=""),col=ifelse(d.temp[,1]=="day","white",ifelse(d.temp[,1]=="night","grey10","grey80")),ylab="mean length (mm)",ylim=c(0,(max(d.temp[,3])*1.5)),xlab="",xlim=c(0,4),main=ord.names[i])

error.bar(predplot,d.temp[,3],d.temp[,4])	

}


#STATISTICS for day/night body length question

#subset to only day and night
d.dn<-subset(d,d$comm.type=="day"|d$comm.type=="night")

#tally number of days and nights each order is represented in. Drop singletons
order.day.night<-aggregate(d.dn$Length,list(d.dn$Order,d.dn$comm.type),FUN=length)
order.day.night <-subset(order.day.night, order.day.night[,3]>1)

#subset to only orders in both time periods, list the orders that fit analysis criteria
sum.occur<-tapply(X= order.day.night[,1],INDEX= order.day.night[,1],FUN=length)
sum.occur<-subset(sum.occur,sum.occur==2)
orders.for.analysis<-rownames(sum.occur)

#function to pull out p-value from lm object
lmp <- function (modelobject) {
    if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
    f <- summary(modelobject)$fstatistic
    p <- pf(f[1],f[2],f[3],lower.tail=F)
    attributes(p) <- NULL
    return(p)
}

#loop to do the analysis for each order, print order, parameter estimate, and p-value
for (i in 1:length(orders.for.analysis)){
	
	d.temp<-subset(d.dn, d.dn$Order==orders.for.analysis[i])
	m.temp<-lm(d.temp$Length~d.temp$comm.type)
	print(orders.for.analysis[i])
	print(m.temp$coefficients)
	print(c("p-value",lmp(m.temp)))
	print("")
	
}

#A first stab at trying to quantify the contribution of family ID vs order*day/night interactions to body length
summary(lm(Length~Order*comm.type+Family,data=d.dn))
summary(lm(Length~trophic.assign,data=d.dn))



#night and day BODY VOLUME, labeled with sample size on x axis
par(mfrow=c(5,4))
#using a loop, fill in each of the 17 plots
for(i in 1: length(unique(day.night.summ[,2]))){
	
d.temp<-subset(day.night.summ,day.night.summ[,2]==ord.names[i])
predplot <- barplot(height=d.temp[,7], names.arg= paste(d.temp[,1],"  ","(",d.temp[,5],")",sep=""),col=ifelse(d.temp[,1]=="day","white",ifelse(d.temp[,1]=="night","grey10","grey80")),ylab="body volume (mm^3)",ylim=c(0,(max(d.temp[,7])*1.5)),xlab="",main=ord.names[i],xlim=c(0,4))
error.bar(predplot,d.temp[,7],d.temp[,8])	
		
}




#plot Hyms and Dips seperately to look at how different families contribute to the total counts
dH<-subset(d,d$Order=="Hymenoptera")

dH.fam.count<-aggregate(dH$Length,list(dH$Family,dH$end.hour),length)
ggplot(dH.fam.count,aes(x=dH.fam.count[,2],y=dH.fam.count[,3],fill=dH.fam.count[,1]))+geom_bar(stat="identity")+labs(title="Hymenoptera abundance")+ scale_x_discrete(breaks=c("2","6","10","14","18","22"), labels=c("2:00", "6:00", "10:00","14:00","18:00","22:00"))+theme(axis.text=element_text(size=14))+labs(fill="Family",x=("Time"),y=("Abundance"))+guides(fill=guide_legend(reverse=TRUE))


dD<-subset(d,d$Order=="Diptera")

dD.fam.count<-aggregate(dD$Length,list(dD$Family,dD$end.hour),length)

ggplot(dD.fam.count,aes(x=dD.fam.count[,2],y=dD.fam.count[,3],fill=dD.fam.count[,1]))+geom_bar(stat="identity")+labs(title="Diptera abundance")+ scale_x_discrete(breaks=c("2","6","10","14","18","22"), labels=c("2:00", "6:00", "10:00","14:00","18:00","22:00"))+theme(axis.text=element_text(size=14))+labs(fill="Family",x=("Time"),y=("Abundance"))+guides(fill=guide_legend(reverse=TRUE))

######################
#when are more families active?
######################

#table of number of families per sampling time
fam.count.time<-ddply(d, .(sample.time), summarize, NumFams = length(unique(ord.fam)))

#pull out number for time of day
fam.count.time$hours<-gsub("(.+).......+","\\1",fam.count.time$sample.time)
fam.count.time<-as.data.frame(fam.count.time)

#calculate mean number of families and se in each time window
num.fams.mean<-tapply(fam.count.time$NumFams,as.factor(fam.count.time$hours),mean)
num.fams.se<-tapply(fam.count.time$NumFams,as.factor(fam.count.time$hours),se)

#stick hours back on as a column, order by time of day for plotting
hours<-as.numeric(rownames(num.fams.mean))
fam.count.summ<-cbind(hours,num.fams.mean, num.fams.se)
fam.count.summ <-fam.count.summ[order(fam.count.summ[,1]),]

#plot mean number of families in each timestep
predplot <- barplot(height= fam.count.summ[,2], names.arg=c("10pm-2am","2am-6am","6am-10am","10am-2pm","2pm-6pm","6pm-10pm"),ylab="families per sampling interval",ylim=c(0,35),xlab="",xlim=c(0,7),main="family level diversity in time")

error.bar(predplot, fam.count.summ[,2], fam.count.summ[,3])	

######################

dP<-subset(d,d$Malaise.Pit=="P")
dM<-subset(d,d$Malaise.Pit=="M")

#total abundance through the day - comparison of pitfall and malaise
par(mfrow=c(3,1))

M.abund<-tapply(dM$Length,dM$end.hour,length)
plot(M.abund,xaxt="n",type="b",pch=16,cex=1.2,xlab="Time of day",ylab="Total Abundance",main="Malaise Abundance",ylim=c(0,500))
axis(side=1,at=c(1,2,3,4,5,6),labels=c("2:00", "6:00", "10:00","14:00","18:00","22:00"))

P.abund<-tapply(dP$Length,dP$end.hour,length)
plot(P.abund,xaxt="n",type="b",pch=16,cex=1.2,xlab="Time of day",ylab="Total Abundance",main="Pitfall Abundance",ylim=c(0,800))
axis(side=1,at=c(1,2,3,4,5,6),labels=c("2:00", "6:00", "10:00","14:00","18:00","22:00"))

abund<-tapply(d$Length,d$end.hour,length)
plot(abund,xaxt="n",type="b",pch=16,cex=1.2,xlab="Time of day",ylab="Total Abundance",main="Total Abundance Abundance",ylim=c(0,1200))
axis(side=1,at=c(1,2,3,4,5,6),labels=c("2:00", "6:00", "10:00","14:00","18:00","22:00"))


#total biomass through the day - comparison of pitfall and malaise
par(mfrow=c(2,1))

M.vol<-tapply(dM$volume,dM$end.hour,sum)
plot(M.vol,xaxt="n",type="b",pch=16,cex=1.2,xlab="Time of day",ylab="Total Volume (mm^3)",main="Malaise Biomass",ylim=c(0,25000))
axis(side=1,at=c(1,2,3,4,5,6),labels=c("2:00", "6:00", "10:00","14:00","18:00","22:00"))

P.vol<-tapply(dP$volume,dP$end.hour,sum)
plot(P.vol,xaxt="n",type="b",pch=16,cex=1.2,xlab="Time of day",ylab="Total Volume (mm^3)",main="Pitfall Biomass",ylim=c(0,25000))
axis(side=1,at=c(1,2,3,4,5,6),labels=c("2:00", "6:00", "10:00","14:00","18:00","22:00"))


#body size through the day - comparison of pitfall and malaise
M.size<-tapply(dM$volume,dM$end.hour,sum)/tapply(dM$Length,dM$end.hour,length)
P.size<-tapply(dP$volume,dP$end.hour,sum)/tapply(dP$Length,dP$end.hour,length)
par(mfrow=c(3,1))

plot(M.size,xaxt="n",type="b",pch=16,cex=1.2,xlab="Time of day",ylab="Average size (mm^3)",main="Malaise body size",ylim=c(0,150))
axis(side=1,at=c(1,2,3,4,5,6),labels=c("2:00", "6:00", "10:00","14:00","18:00","22:00"))

plot(P.size,xaxt="n",type="b",pch=16,cex=1.2,xlab="Time of day",ylab="Average size (mm^3)",main="Pitfall body size",ylim=c(0,50))
axis(side=1,at=c(1,2,3,4,5,6),labels=c("2:00", "6:00", "10:00","14:00","18:00","22:00"))

size<-tapply(d$volume,d$end.hour,sum)/tapply(d$Length,d$end.hour,length)
plot(size,xaxt="n",type="b",pch=16,cex=1.2,xlab="Time of day",ylab="Average size (mm^3)",main="Body size",ylim=c(0,50))
axis(side=1,at=c(1,2,3,4,5,6),labels=c("2:00", "6:00", "10:00","14:00","18:00","22:00"))
#######################
#description of most dominant taxa through the day
#######################
#most abundant taxon in each sample
fam.count<-aggregate(d$Length,list(d$ord.fam,d$Malaise.Pit,d$sample.time),FUN=length)
fam.count$vars<-paste(fam.count$Group.2,fam.count$Group.3,sep=".")
abundant.taxa<-ddply(fam.count,.(vars),subset,x==max(x))

#most massive taxon in each sample
fam.mass<-aggregate(d$volume,list(d$ord.fam,d$Malaise.Pit,d$sample.time),FUN=sum)
fam.mass$vars<-paste(fam.mass$Group.2, fam.mass$Group.3,sep=".")
massive.taxa<-ddply(fam.mass,.(vars),subset,x==max(x))

######################
#family accumulation curve

#accum<-matrix(0,nrow=length(d$ord.fam),ncol=10)
#for(h in 1:10){
#d.rand<-d[sample(nrow(d)),]
#new.row<-1:length(d.rand$ord.fam)
#d.rand<-cbind(d.rand,new.row)
#for(i in 1:length(d.rand$ord.fam)){
#	temp<-subset(d.rand,d.rand$new.row<i+1)
#	accum[i,h]<-length(unique(temp$ord.fam))

#}
#}
#accum.to.fit<-rowMeans(accum)
#accum.to.fit<-cbind(1:length(accum.to.fit),accum.to.fit)
#plot(accum.to.fit,xlab="Individuals collected",ylab="Family count",main="Family accumulation")

##########################
#plot all data together, summing across the 5 days
##########################
trophic.counts<-aggregate(d$Length,list(d$end.hour,d$trophic.assign),FUN=length)


#plot absolute total abundance
ggplot(trophic.counts,aes(x= trophic.counts[,1],y= trophic.counts[,3],fill= trophic.counts[,2]))+geom_bar(stat="identity")+labs(title="Trophic abundance")+ scale_x_discrete(breaks=c("2","6","10","14","18","22"), labels=c("2:00", "6:00", "10:00","14:00","18:00","22:00"))+theme(axis.text=element_text(color="black",size=18))+labs(fill="Trophic position",x=("Time"),y=("Abundance"))+guides(fill=guide_legend(reverse=TRUE))



#plot proportion total abundance
ggplot(trophic.counts, aes(x=trophic.counts[,1],y=trophic.counts[,3],group=trophic.counts[,2],fill=trophic.counts[,2])) + geom_area(position="fill") + scale_x_discrete(breaks=c("2","6","10","14","18","22","26","30","34","38","42","46","50"), labels=c("2:00", "6:00", "10:00","14:00","18:00","22:00","2:00","6:00", "10:00","14:00","18:00","22:00","2:00"))+labs(fill="Trophic position",size=18,x=("Time"),y=("Relative abundance"))+theme(axis.text=element_text(size=16,color="black"))+guides(size=18,fill=guide_legend(reverse=TRUE))



#volume -total within each trophic level
trophic.volumes<-aggregate(d$volume,list(d$end.hour,d$trophic.assign),FUN=sum)
ggplot(trophic.volumes,aes(x= trophic.volumes[,1],y= trophic.volumes[,3],fill= trophic.volumes[,2]))+geom_bar(stat="identity")+labs(title="Trophic biomass distribution")+ scale_x_discrete(breaks=c("2","6","10","14","18","22"), labels=c("2:00", "6:00", "10:00","14:00","18:00","22:00"))+theme(axis.text=element_text(size=14))+labs(fill="Trophic position",x=("Time"),y=("Biomass (mm^3)"))

#volume - proportion within each trophic level
ggplot(trophic.volumes, aes(x=trophic.volumes[,1],y=trophic.volumes[,3],group=trophic.volumes[,2],fill=trophic.volumes[,2])) + geom_area(position="fill")+labs(title="Proportion biomass by trophic position and time of day") + scale_x_discrete(breaks=c("2","6","10","14","18","22","26"), labels=c("2:00", "6:00", "10:00","14:00","18:00","22:00","2:00"))+theme(axis.text=element_text(size=14))+labs(fill="Trophic position",x=("Time"),y=("Proportion biomass"))

###############
#for ent soc presentation
dev.off()

trophic.counts<-aggregate(d$Length,list(d$end.hour,d$trophic.assign),FUN=length)

t.count.dub <- cbind((trophic.counts[,1]+24),trophic.counts[,c(2,3)])
colnames(t.count.dub)<-colnames(trophic.counts)
t.count.dub <-rbind(trophic.counts, t.count.dub)


unique(t.count.dub[,2])

#colors<-c( 
#rgb(161,158,103,maxColor=255) ,
#rgb(34,59,14,maxColor=255)  ,
#rgb(158,212,114,maxColor=255)   ,
#rgb(255,172,101,maxColor=255) ,   
#rgb(114,130,212,maxColor=255),
#rgb(172,97,95,maxColor=255),
#rgb(114,212,186,maxColor=255),
#rgb(161,54,42,maxColor=255)  ,
#rgb(161,42,127,maxColor=255)
#)

#plot proportion total abundance

ggplot(t.count.dub, aes(x= t.count.dub[,1],y= t.count.dub[,3],group= t.count.dub[,2],fill= t.count.dub[,2])) + geom_area(position="fill")+scale_x_discrete(breaks=c("2","6","10","14","18","22","26","30","34","38","42","46","50"), labels=c("2:00", "6:00", "10:00","14:00","18:00","22:00","2:00","6:00", "10:00","14:00","18:00","22:00","2:00"))+labs(fill="Trophic position",size=18,x=("Time"),y=("Relative abundance"))+theme(axis.text=element_text(size=16,color="black"), legend.text = element_text(size = 16, face = "bold"), legend.title = element_text(size=18, face="bold"))+guides(size=18,fill=guide_legend(reverse=TRUE))

theme(legend.text = element_text(size = 16, face = "bold"))


############################
#plot data by sample, to illustrate sample relative homogeneity within timesteps
############################
#abundance
#facet by time of day, stack by trophic, display replicate days along x-axis
trophic.counts<-aggregate(d$Length,list(d$end.hour,d$End.Date,d$trophic.assign),FUN=length)
colnames(trophic.counts)<-c("end.hour","End.Date","trophic.assign","abund")
p<-ggplot(trophic.counts,aes(x= End.Date,y= abund,fill= trophic.assign))+geom_bar(stat="identity")+labs(title="Trophic abundance")+theme(axis.text=element_text(color="black",size=18))+labs(fill="Trophic position",x=("Time"),y=("Abundance"))+guides(fill=guide_legend(reverse=TRUE))
p+facet_wrap(~end.hour,ncol=2)


#volume
#facet by time of day, stack by trophic, display replicate days along x-axis
trophic.volumes<-aggregate(d$volume,list(d$end.hour,d$End.Date,d$trophic.assign),FUN=sum)
colnames(trophic.volumes)<-c("end.hour","End.Date","trophic.assign","volume")
p<-ggplot(trophic.volumes,aes(x= End.Date,y= volume,fill= trophic.assign))+geom_bar(stat="identity")+labs(title="Trophic volume")+theme(axis.text=element_text(color="black",size=18))+labs(fill="Trophic position",x=("Time"),y=("Volume (mm^3)"))+guides(fill=guide_legend(reverse=TRUE))
p+facet_wrap(~end.hour,ncol=2)
###########################
#family-time curve
###########################
length(unique(d$Family))
#for this number and less sample 1:30, then get unique ord.fam
time.list<-sort(unique(d$time))
d$sample.num <- match(d$time,time.list )

fam.accum<-vector(length=length(time.list))
for(i in 1:length(unique(time.list))){
	temp <- d[d$sample.num<=i,]
	fam.accum[i] <- length(unique(temp$ord.fam))
}

plot(fam.accum)

############################
#treemaps for proportions
############################
#try a plot for predators only
#need columns for time of day, trophic position, order, family, abundance
tree.map.trophic <- function(name.of.group){
temp <- d[d$trophic.assign==name.of.group,]

for.tree<-with(temp, aggregate(ord.fam, by=list(End.Time , Order ,Family), FUN=length))

colnames(for.tree) <- c("End.Time" , "Order" ,"Family",name.of.group)

treemap(for.tree, index = c("End.Time","Family"),vSize=name.of.group,type="index", fontsize.labels=c(25,14),fontsize.title=25,bg.labels=0,overlap.labels=1,force.print=TRUE,ymod.labels=c(.05,0))
}

dev.new(width=10, height=14)

tree.map.trophic("Parasitoid")

tree.map.trophic("Detritivore")

d[d$trophic.assign=="Detritivore",]

for.tree<-with(d, aggregate(ord.fam, by=list(End.Time , trophic.assign , ord.fam), FUN=length))

colnames(for.tree) <- c("End.Time" , "trophic.assign" ,"Family","Time Communities")


dev.new(width=50, height=40)
treemap(for.tree, index = c("End.Time","trophic.assign","Family"),vSize="Time Communities",type="index", fontsize.labels=c(100,60,14),fontsize.title=80,bg.labels=0,overlap.labels=1,force.print=TRUE,ymod.labels=c(0,.5,0))

for.tree<-with(d, aggregate(volume, by=list(End.Time , trophic.assign , ord.fam), FUN=sum))

colnames(for.tree) <- c("End.Time" , "trophic.assign" ,"Family","Time Communities")

dev.new(width=50, height=40)
treemap(for.tree, index = c("End.Time","trophic.assign","Family"),vSize="Time Communities",type="index", fontsize.labels=c(100,60,14),fontsize.title=80,bg.labels=0,overlap.labels=1,force.print=TRUE,ymod.labels=c(0,.5,0))





