#setup
        #clear all objects from R console
          rm(list=ls())
          graphics.off()
          pardefault <- par(no.readonly = T)
        
        #set working directory
          setwd("/Users/mmcmunn/Desktop/GitHub/Mendocino")
        
        #load packages
          library("codyn")
          library("ggplot2")
          library("vegan")
          library("reshape")
          library("plyr")
          library("treemap")
          library("mvabund")
          library("gridExtra")
          library("tidyr")
          
          cbbPalette <- c("#009E73", "#e79f00", "#9ad0f3", "#0072B2", "#D55E00", 
                          "#CC79A7", "#F0E442")
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
          se <- function(x) sqrt( var(x[ !is.na(x) ] ) / length(x[ !is.na(x) ] ))
          
#data cleaning
        #create new vector for unique families
          d$ord.fam<-paste(d$Order,d$Family,sep=".")
        
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
            
        #create a categorical variable for day, night, crepuscular
            
            d$comm.type<-ifelse(d$End.Time=="10:00"|d$End.Time=="22:00","crepuscular",NA)
            d$comm.type<-ifelse(d$End.Time=="14:00"|d$End.Time=="18:00","day",d$comm.type)
            d$comm.type<-ifelse(d$End.Time=="2:00"|d$End.Time=="6:00","night",d$comm.type)
            d$comm.type<-as.factor(d$comm.type)
            combinations  <- as.data.frame( xtabs(~ Order + End.Time + comm.type, d))[ , c("Order" , "End.Time", "comm.type")]
            dExpand <- merge(combinations , d, all.x = TRUE)   
            
        #make End.Time and ordered factor so it behaves in plots  
            d <- within(d,  End.Time <- factor(End.Time, levels=c("2:00", "6:00", "10:00", "14:00" , "18:00" , "22:00") ))
            
#summary statistics, tables        
        #average # of arthropods in each sample
             mean(table(d$sampleID))
        
        #which taxa were most abundant
            head(sort(table(d$ord.fam), decreasing = TRUE) , 20)    
    
            
            
########## FIGURE 2 ##############
### Temporal diversity indices ###
#when are more families active?

#table of number of families per sampling time, standard error
divTime <-  with( with(d ,  aggregate(x = ord.fam, by = list(ord.fam = ord.fam,End.Time = End.Time), FUN = length )) ,
                  aggregate( x , list(End.Time) , function(x){c(length = length(x) , se = se(x)) }))
#plot diversity by time  
p1 <- ggplot(divTime , aes(x = Group.1 , y = x[,"length"]), group = 1) + labs(x = "Time", y = "Average family diversity") 
p1 + geom_errorbar(aes(ymin=x[,"length"]-x[,"se"], ymax=x[,"length"]+x[,"se"] , position= Group.1), width=.1) + geom_point()

### Panel B: turnover ###
head(d)
sampleCounts <-  count(d , c("ord.fam" , "dateTimeEnd"))
sampleCounts$timeNum <- as.numeric(sampleCounts$dateTimeEnd)
mendoTurnover <- turnover(df = sampleCounts,  
                         time.var = "timeNum",  
                         species.var = "ord.fam", 
                         abundance.var = "freq",
                         metric = "total")


mendoAppear <- turnover(df = sampleCounts,  
                                            time.var = "timeNum",  
                                            species.var = "ord.fam", 
                                            abundance.var = "freq",
                                            metric = "appearance")

mendoDisappear <- turnover(df = sampleCounts,  
                           time.var = "timeNum",  
                           species.var = "ord.fam", 
                           abundance.var = "freq",
                           metric = "disappearance")

mendoTurnover$time <- as.POSIXct(mendoTurnover$timeNum, origin="1970-01-01")
reps <- rep(1:6, 5)
reps <- reps[-length(reps)]
mendoTurnover$group <- reps
mendoTurnover$appear <- mendoAppear[,1]
mendoTurnover$disappear <- mendoDisappear[,1]


means<- with(mendoTurnover, tapply(total, group, mean))
ses<-with(mendoTurnover, tapply(total, group, se))
plot(means, xaxt = "n" , ylim = c(0,1), type = "b", xlab = "time community transition", ylab = "proportion family turnover", pch = 16)
arrows(x0 = 1:6, y0 = means - ses, y1 = means + ses, angle = 90, code = 3)

means<- with(mendoTurnover, tapply(appear, group, mean))
ses<-with(mendoTurnover, tapply(appear, group, se))
points(means, type = "b", col = "blue", pch = 16)
arrows(x0 = 1:6, y0 = means - ses, y1 = means + ses, angle = 90, code = 3, col = "blue")

means<- with(mendoTurnover, tapply(disappear, group, mean))
ses<-with(mendoTurnover, tapply(disappear, group, se))
points(means, type = "b", col = "red", pch = 16)
arrows(x0 = 1:6, y0 = means - ses, y1 = means + ses, angle = 90, code = 3, col = "red")

axis(1, labels = c("22->2","2->6","6->10", "10->14", "14->18", "18->22"), at = 1:6)
legend("topright", legend = c("total", "appearances", "disappearances"), col = c("black", "blue", "red"), pch = 16, lwd = 1)
?legend

summary(aov(mendoTurnover$total~mendoTurnover$group))
summary(aov(mendoTurnover$disappear~mendoTurnover$group))
summary(aov(mendoTurnover$appear~mendoTurnover$group))



plot(mendoTurnover[,1], type = "b", ylim = c(0,1))
lines(mendoAppear[,1] , col = "blue")
lines(mendoDisappear[,1] , col = "red")

#Format a compiled data frame
KNZ_turnover$metric<-"total"
names(KNZ_turnover)[1]="turnover"

KNZ_appearance$metric<-"appearance"
names(KNZ_appearance)[1]="turnover"

KNZ_disappearance$metric<-"disappearance"
names(KNZ_disappearance)[1]="turnover"

KNZ_allturnover<-rbind(KNZ_turnover, KNZ_appearance, KNZ_disappearance)

#Create the graph
turn.graph <- ggplot(KNZ_allturnover, aes(x=year, y=turnover, color=metric)) + 
  geom_line(size = 1) + 
  facet_wrap(~replicate) + 
  theme_bw() + 
  theme(legend.position="bottom")

### Panel C: mean rank shifts ###

#Run the rank shift code
KNZ_rankshift <- rank_shift(df=collins08, 
                            time.var = "year", 
                            species.var = "species",
                            abundance.var = "abundance", 
                            replicate.var = "replicate")

#Select the final time point from the returned time.var_pair
KNZ_rankshift$year <- as.numeric(substr(KNZ_rankshift$year_pair, 6,9))

# Create the graph
rankshift.graph <- ggplot(KNZ_rankshift, aes(year, MRS)) + 
  geom_line(size = 1) + 
  facet_wrap(~replicate) +
  theme_bw() 


### Panel D: rate change ###

# Run the rate change code
comm.res <- rate_change_interval(collins08,   
                                 time.var= "year",    
                                 species.var= "species",  
                                 abundance.var= "abundance", 
                                 replicate.var = "replicate")

# Create the graph
rate.graph<-ggplot(comm.res, aes(interval, distance, group = replicate)) + 
  geom_point() + 
  facet_wrap(~replicate) + 
  stat_smooth(method = "lm", se = F, size = 1) +
  theme_bw() 

### Put the panel graph together ###
fig2 <- grid.arrange(rich.graph +
                       labs(x="Year", y=expression(paste("Richness (no / 10", m^2,")"))) +
                       theme(strip.text.x = element_text(size = 14),
                             strip.background = element_blank()) +
                       theme( plot.margin=unit(c(0,1,0,0), "cm")),
                     
                     turn.graph + 
                       labs(x="Year", y="Turnover") +
                       theme(legend.position="none") +
                       theme(strip.background = element_blank(),
                             strip.text.x = element_blank()) +
                       scale_color_manual( values = cbbPalette) + 
                       theme( plot.margin=unit(c(0,1,0,0), "cm")),
                     
                     rankshift.graph + 
                       labs(x="Year", y="Mean rank shift") +
                       theme(strip.background = element_blank(),
                             strip.text.x = element_blank()) +
                       theme( plot.margin=unit(c(0,1,0,0), "cm")), 
                     
                     rate.graph +  
                       labs(x="Time interval", y="Euclidean distance") +
                       theme(strip.background = element_blank(),
                             strip.text.x = element_blank()) +
                       theme( plot.margin=unit(c(0,1,0,0), "cm")), 
                     
                     ncol=1)


            
#make community matrices
        #subsets by malaise and pitfall
          dM <- d[d$Malaise.Pit=="M" , ]
          dP <- d[d$Malaise.Pit=="P" , ]
    
        #abundance within samples by family
          abund.FxS <- as.data.frame.matrix(table(d$timeChar, d$ord.fam))
          M.abund.FxS <- as.data.frame.matrix(table(dM[ ,"timeChar"], dM[, "ord.fam"]))
          P.abund.FxS <- as.data.frame.matrix(table(dP[ ,"timeChar"], dP[, "ord.fam"]))
          
        #abundance within samples by trophic position
          abund.TxS <-as.data.frame.matrix(table(d$timeChar,d$trophic.assign))
          M.abund.TxS <-as.data.frame.matrix(table(dM$timeChar,dM$trophic.assign))
          P.abund.TxS <-as.data.frame.matrix(table(dP$timeChar,dP$trophic.assign))
          
        #biomass within samples by family
          vol.FxS <- as.data.frame.matrix(t(with(d, tapply(volume, list(ord.fam, timeChar), sum))))
          vol.FxS[is.na(vol.FxS)] <- 0
        
        #biomass within samples by trophic position
          vol.TxS <- as.data.frame.matrix(t(with(d, tapply(volume, list(trophic.assign, timeChar), sum))))
          vol.TxS[is.na(vol.TxS)] <- 0


#multivariate abundance models using mvabund package
#describe which families and trophic positions vary by time of day
          #by family, abundance
                    mAbund <- manyglm(as.matrix(abund.FxS) ~ as.factor(clim$dateTimeEnd$hour), family = "negative.binomial")
                    mAbund.sig <- anova(mAbund, p.uni = "adjusted")
                    
                  #time of day is significant (.001) as a factor determining community composition
                  #check residuals vs fitted
                    plot.manyglm(mAbund)
                  #check mean variance relationship
                    meanvar.plot(abund.FxS[,] ~ as.factor(clim$dateTimeEnd$hour))
                  
                  #how many taxa varied in abundance significantly by time of day
                    sum(mAbund.sig$uni.p[2,]<=.05)
                    which(mAbund.sig$uni.p[2,]<=.05)
                  #some notable taxa that were abundant, but did not vary by time of day - ants, cecidomyiids
                  
                  #a closer look at taxa driving changes over the course of the day
                  #get top 15 likelihood statistics
                    top.fifteen.uni <- names(sort(-mAbund$two.loglike, decreasing=TRUE)[1:15])
                    top.fifteen.lik <- sort(-mAbund$two.loglike, decreasing=TRUE)[1:15]
                    top.fifteen.values <- top.fifteen.uni
          
          
          #by trophic, abundance
                    mTroph <- manyglm(as.matrix(abund.TxS) ~ as.factor(clim$dateTimeEnd$hour), family = "negative.binomial")
                    mTroph.sig <- anova(mTroph, p.uni = "adjusted")
                  
                  #time of day is significant (.001) as a factor determining community composition
                  #check residuals vs fitted
                    plot.manyglm(mTroph)
                  #check mean variance relationship
                    meanvar.plot(abund.TxS[,] ~ as.factor(clim$dateTimeEnd$hour))
                  
                  #how many taxa varied in abundance significantly by time of day
                    sum(mTroph.sig$uni.p[2,]<=.05)
                    which(mTroph.sig$uni.p[2,]<=.05)

#with abiotic covariates                    
           #by trophic, abundance
                    mTroph <- manyglm(as.matrix(abund.TxS) ~ clim$light.sun + clim$mean.T + clim$wind.max + clim$perc.clouds, family = "negative.binomial")
                    mTroph.sig <- anova(mTroph, p.uni = "adjusted")
                    summary(mTroph)
                    
                    #time of day is significant (.001) as a factor determining community composition
                    #check residuals vs fitted
                    plot.manyglm(mTroph)
                    #check mean variance relationship
                    meanvar.plot(abund.TxS[,] ~ as.factor(clim$dateTimeEnd$hour))
                    
                    mTroph.sig$uni.p
                    
                    
#a function to animate each samnple in a NMDS plot
ordination.animation <- function(comm.matrix, file.ext){
              
                set.seed(1231231)
                ord.all= metaMDS(comm=as.matrix(comm.matrix),distance="bray")
                xy.coord <- as.data.frame(scores(ord.all))
                xy.coord$hour <- strptime(rownames(xy.coord), format = "%Y-%m-%d %H:%M:%S")$hour
                cols <- c(rep(c( "purple", "black",  "gray50", "red", "yellow", "green"), 5))
            
            #for gif
                    setwd("animations")
                            #frame 31, with convex hulls around crepuscular, diurnal, nocturnal
                            name<-paste("ord.plot",file.ext,31,"png",sep=".")
                            png(name, width = 600, height = 600)
                            par(oma = c(1,1,0,0))
                            plot(x=xy.coord[,1], y=xy.coord[,2],col=cols,pch=16,cex=2.5,cex.lab=1.5,cex.axis=1.5,
                                 xlim = c(min(xy.coord[,1])*1.3 , max(xy.coord[,1])*1.3),
                                 ylim = c(min(xy.coord[,2])*1.3 , max(xy.coord[,2])*1.5),
                                 xlab = "NMDS 1", ylab = "NMDS 2", main = paste(file.ext, " by collection time", sep = ""))
                              legend("topleft",legend=c("10pm-2am","2am-6am","6am-10am","10am-2pm","2pm-6pm","6pm-10pm"),
                                     pch=rep(16,6),col=c("black","gray50","red","yellow","green","purple"),cex=1.5)
                          
                            #outer hulls
                              crep.points <- xy.coord[which(xy.coord$hour==22 | xy.coord$hour==10),]
                              polygon(crep.points[chull(crep.points),])
                              text(mean(crep.points[,1]), mean(crep.points[,2])-.5, "crepuscular", cex = 1.5)
                              
                              day.points <- xy.coord[which(xy.coord$hour==14 | xy.coord$hour==18),]
                              polygon(day.points[chull(day.points),])
                              text(mean(day.points[,1]), mean(day.points[,2])-.5, "diurnal", cex = 1.5)
                              
                              night.points <- xy.coord[which(xy.coord$hour==2 | xy.coord$hour==6),]
                              polygon(night.points[chull(night.points),])
                              text(mean(night.points[,1]), mean(night.points[,2])-.5, "nocturnal", cex = 1.5)
                          dev.off()
                    #sample level loop to produce figues
                    for(i in 1:30){
                    	name<-paste("ord.plot",file.ext,i,"png",sep=".")
                          	png(name, width = 600, height = 600)
                          	par(oma = c(1,1,0,0))
                            	plot(xy.coord[ , c("NMDS1", "NMDS2")],type="n",col=cols,pch=16,cex=1.5,cex.lab=1.5,cex.axis=1.5,
                            	xlim =c(min(xy.coord[,1])*1.3 , max(xy.coord[,1])*1.3),
                            	ylim = c(min(xy.coord[,2])*1.3 , max(xy.coord[,2])*1.5),
                            	xlab = "NMDS 1", ylab = "NMDS 2", main = paste(file.ext, " by collection time", sep = "")
                            	)
                          points(x=xy.coord[1:i,1], y=xy.coord[1:i,2],col=cols[1:i],pch=16,cex=2.5) 
                          legend("topleft",legend=c("10pm-2am","2am-6am","6am-10am","10am-2pm","2pm-6pm","6pm-10pm"),pch=rep(16,6),col=c("black","gray50","red","yellow","green","purple"),cex=1.5)
                              
                              #consecutive data points loop to plot all lines between previous samples   
                                  for (h in 1:i ){
                                  	segments(x0 = xy.coord[ifelse(h >1, h-1, h) , 1 ], y0= xy.coord[ifelse(h >1, h-1, h) , 2], x1 = xy.coord[h, 1], y1 = xy.coord[h, 2]   )	
                                  }
                                  arrows(x0 = xy.coord[ifelse(i >1, i-1, i) , 1 ], y0= xy.coord[ifelse(i >1, i-1, i) , 2], x1 = xy.coord[i, 1], y1 = xy.coord[i, 2]   , lwd=6 )
                                  dev.off()
                                  }
                                  setwd("..") }
            
#apply the NMDS animation function
      #family level abundance, malaise and pitfall bulked, then separatly
          ordination.animation(comm.matrix = abund.FxS, file.ext = "Family abundance")
          ordination.animation(comm.matrix = P.abund.FxS, file.ext = "Pitfall family abundance")
          ordination.animation(comm.matrix = M.abund.FxS, file.ext = "Malaise family abundance")
      
      #trophic abundance, malaise and pitfall bulked, then separatly
          ordination.animation(comm.matrix = abund.TxS, file.ext = "Trophic abundance")
          ordination.animation(comm.matrix = P.abund.TxS, file.ext = "Pitfall trophic abundance")
          ordination.animation(comm.matrix = M.abund.TxS, file.ext = "Malaise trophic abundance")
      
      #by volume, family and trophic level
          ordination.animation(comm.matrix = vol.FxS, file.ext = "Family volume")
          ordination.animation(comm.matrix = vol.FxS, file.ext = "Trophic volume")

          #PERMANOVA for family abundance by abiotic variables
    adonis(abund.FxS ~ clim$light.sun + clim$mean.T + clim$wind.max + clim$perc.clouds)

#PERMANOVA for trophic abundance  by abiotic variables
    adonis(abund.TxS ~ clim$light.sun + clim$mean.T + clim$wind.max + clim$perc.clouds)

#PERMANOVA for family volume  by abiotic variables
    adonis(vol.FxS ~ clim$light.sun + clim$mean.T + clim$wind.max + clim$perc.clouds)
    
#PERMANOVA for trophic volume  by abiotic variables
    adonis(vol.TxS ~ clim$light.sun + clim$mean.T + clim$wind.max + clim$perc.clouds)        

    permTable <- function(matrix) {
      modelResult <- adonis(matrix ~ as.factor(clim$dateTimeEnd$hour))
      pValue <- modelResult$aov.tab["Pr(>F)"][1,]
      R2value <- modelResult$aov.tab["R2"][1,]
      Fstat <- modelResult$aov.tab$F.Model[1]
      data.frame("pValue" = pValue ,"R2value" = R2value, "Fstat" = Fstat )
    }
    ?adonis

PERMsummaries  <- as.data.frame(t( sapply(     list("individuals in families" = abund.FxS, 
             "individuals in feeding strategy" = abund.TxS,
             "volume in families" = vol.FxS, 
             "volume in feeding strategies" = vol.TxS)
                        , permTable)))

PERMsummaries


    #plot average body length by order
      #mean and se within orders
      lenOrd <- data.frame (mean = with(d, tapply(Length, Order, mean)),se= with(d, tapply(Length, Order, se)))
      lenOrd <- lenOrd[order(lenOrd$mean, decreasing=TRUE) , ]
      lenOrd$order <- factor(rownames(lenOrd), levels = rownames(lenOrd)) 
      
      p1 <- ggplot(lenOrd)+aes(x= order, y = mean) + geom_bar(stat = "identity") 
      p1 <-  p1+ geom_errorbar(aes(ymin = mean-se, ymax=mean+se), width = 0.2)+theme(axis.text.x = element_text(angle = 45, vjust= 1,hjust=1))
      p1 <- p1 + theme(axis.text=element_text(size=14,face="bold")) + labs(y = "mean body length (mm)", x = "Order")
      p1 <- p1+ theme(axis.title = element_text(size=18,face="bold")) + scale_y_continuous( expand = c(0, 0), limits = c(0,20)) 
      p1
      
#calculate mean and se of mean sizes, volumes, and sample sizes within each order and time block

      
      summary.table <- with (dExpand, aggregate(x = cbind(Length, volume),
                                  by = list(Order = Order, endTime = End.Time), 
                                  FUN= function(x) {c(mean = mean(x, na.rm = TRUE), se = se(x), n = length(x))}
                                  ))
      summary.table[,"Length"][,"mean"][  is.nan(summary.table[,"Length"][,"mean"]) ] <- NA

#reorder by order then by time
      summary.table <- summary.table[ order( summary.table[,"Order"] , summary.table[ , "endTime"] ) , ]

  #body length   
      #body volume
      plots <- list()  
      #using a loop, fill in each of the 17 plots
      for(i in unique(d$Order) ) {
        d.temp <- summary.table[summary.table$Order == i , ]
        
        p1 <- ggplot(d.temp)+aes(x= endTime, y = Length[,"mean"]) +
          geom_bar(stat = "identity") + 
          ggtitle(i) + xlab("") + ylab("mean length (mm)")
        plots[[i]] <- p1
        
      }
      grid.arrange(grobs = plots)
      head(summary.table)
  #body volume
            plots <- list()  
            #using a loop, fill in each of the 17 plots
            for(i in unique(d$Order) ) {
              d.temp <- summary.table[summary.table$Order == i , ]
              
              p1 <- ggplot(d.temp)+aes(x= endTime, y = volume[,"mean"]) +
                geom_bar(stat = "identity") + 
                ggtitle(i) + xlab("") + ylab("mean volume (mm^3)")
              plots[[i]] <- p1
              
            }
            grid.arrange(grobs = plots)
            
  #abundance
         plots <- list()  
        #using a loop, fill in each of the 17 plots
        for(i in unique(d$Order) ) {
          d.temp <- summary.table[summary.table$Order == i , ]
          
          p1 <- ggplot(d.temp)+aes(x= endTime, y = Length[,"n"]) +
            geom_bar(stat = "identity") + 
            ggtitle(i) + xlab("") + ylab("abundance")
          plots[[i]] <- p1
          
        }
         grid.arrange(grobs = plots)
         
#summarize by day, night, crepuscular    
  day.night.summ <- with (dExpand, aggregate(x = cbind(Length, volume),
                          by = list(Order = Order, comm.type = comm.type), 
                          FUN= function(x) {c(mean = mean(x, na.rm = TRUE), se = se(x), n = length(x))}
        ))

  
#night and day BODY LENGTH, labeled with sample size on x axis
#using a loop, fill in each of the 17 plots
      plots <- list()
      for(i in unique(day.night.summ[,"Order"])){
      	
                    d.temp <- day.night.summ[ day.night.summ[ ,"Order" ]==i , ]
                    p1 <- ggplot(d.temp)+aes(x= comm.type, y = Length[,"mean"]) +
                      geom_bar(stat = "identity", fill = myColors[!is.na(d.temp[, "Length"][ , "mean"])]) + 
                      ggtitle(i) + xlab("") + ylab("mean length (mm)")
                    p1 <-  p1+ geom_errorbar(aes(ymin = Length[ , "mean"]-Length[ , "se"], ymax=Length[ , "mean"]+Length[ , "se"]), width = 0.2)
                    plots[[i]] <- p1
              }
              grid.arrange(grobs = plots)


#STATISTICS for day/night body length question

      #subset to only day and night, subset orders occuring in both
      d.dn<-subset(d,d$comm.type=="day"|d$comm.type=="night")
      orders.for.analysis <- names( which( rowSums( table(d.dn$Order, d.dn$comm.type) > 0 ) == 2 ) )

      #function to pull out p-value from lm object
      lmp <- function (modelobject) {
          if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
          f <- summary(modelobject)$fstatistic
          p <- pf(f[1],f[2],f[3],lower.tail=F)
          attributes(p) <- NULL
          return(p)
      }
      
      #loop to do the analysis for each order, print order, parameter estimate, and p-value
        coef <- list()
          for (i in orders.for.analysis){
          	d.temp <- d.dn[d.dn$Order == i , ]
          	m.temp <- lm(d.temp$Length ~ d.temp$comm.type)
          	coef[[i]][["p.value"]] <-  lmp(m.temp)
          	coef[[i]][["Intercept"]] <- m.temp$coefficients[1]
          	coef[[i]][["comm.typenight"]] <- m.temp$coefficients[2]
          }


        
#plot Hyms and Dips seperately to look at how different families contribute to the total counts
        
    #plot hymenopteran families by time
    Hym.fam.count <- with (dH<-d[ d$Order=="Hymenoptera" , ] ,  
                  aggregate(x = ord.fam, by = list(Family = Family,End.Time = End.Time), FUN = length ))
    
    ggplot(Hym.fam.count , aes(x=End.Time , y = x ,fill=Family)) + geom_bar(stat="identity") + labs(title="Hymenoptera abundance")+
      theme(axis.text=element_text(size=14))+labs(fill="Family",x=("Time"),y=("Abundance"))+guides(fill=guide_legend(reverse=TRUE))

    #plot hymenopteran families by time
      Dip.fam.count <- with (dH<-d[ d$Order=="Diptera" , ] ,
                  aggregate(x = ord.fam, by = list(Family = Family,End.Time = End.Time), FUN = length ))
      
      ggplot(Dip.fam.count , aes(x=End.Time , y = x ,fill=Family)) + geom_bar(stat="identity") + labs(title="Diptera abundance")+
        theme(axis.text=element_text(size=14))+labs(fill="Family",x=("Time"),y=("Abundance"))+guides(fill=guide_legend(reverse=TRUE))


        
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
setwd("figures")
png(file = "family accumulation")
plot(specaccum(abund.FxS), 
    xlab = "", ylab = "",
    main = "Family accumulation",cex = 2)
graphics.off()
setwd("..")

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

setwd("figures")
png(file = "parisitoid_map", res = 1000)
?png
tree.map.trophic("Parasitoid")
dev.off()

png(file = "detritivore_map", width = 1300, height = 1300, font = "Times")
tree.map.trophic("Detritivore")
dev.off()

png(file = "predator_map", width = 1300, height = 1300)
tree.map.trophic("Predator")
dev.off()

png(file = "herbivore_map", width = 1300, height = 1300)
tree.map.trophic("Herbivore")
dev.off()


for.tree<-with(d, aggregate(ord.fam, by=list(End.Time , trophic.assign , ord.fam), FUN=length))

colnames(for.tree) <- c("End.Time" , "trophic.assign" ,"Family","Time Communities")


dev.new(width=50, height=40)
treemap(for.tree, index = c("End.Time","trophic.assign","Family"),vSize="Time Communities",type="index", fontsize.labels=c(100,60,14),fontsize.title=80,bg.labels=0,overlap.labels=1,force.print=TRUE,ymod.labels=c(0,.5,0))

for.tree<-with(d, aggregate(volume, by=list(End.Time , trophic.assign , ord.fam), FUN=sum))

colnames(for.tree) <- c("End.Time" , "trophic.assign" ,"Family","Time Communities")

dev.new(width=50, height=40)
treemap(for.tree, index = c("End.Time","trophic.assign","Family"),vSize="Time Communities",type="index", fontsize.labels=c(100,60,14),fontsize.title=80,bg.labels=0,overlap.labels=1,force.print=TRUE,ymod.labels=c(0,.5,0))





