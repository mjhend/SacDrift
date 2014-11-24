setwd("C:/Users/mark.henderson/Desktop/Sacramento River drift data/R analysis/Lookup tables")

#Aquatic biology associates hierarchy and regression values
#ran this to compare my code with the results on the long output from Aq Biol Associates spreadsheet

#AqBiolLkup <- read.csv(file="AqBiolHier-LM.csv", stringsAsFactors=F,header=T)
#danger <- read.csv(file="AqBiolRecord-Danger.csv", stringsAsFactors=F, header=T)
#danger$Date<-as.Date(danger$Date, "%m/%d/%Y")

#dangerBio <- merge(danger, AqBiolLkup, by=c("Taxon", "Stage"))
#dangerBio <- dangerBio[,c("Taxon", "Stage", "Site", "Date", "Replicate", "Length", "Abundance", "a", "b")]
#dangerBio$biomass <- with(dangerBio, Abundance*a*Length^b)
#dangerBio <- aggregate(biomass~Taxon+Stage+Site+Date+Replicate, data=dangerBio, FUN=sum)
#dangerBio <- dangerBio[with(dangerBio,order(Taxon, Stage, Date, Replicate)),]

#Henderson data and linear regression values
deployData <- read.csv(file="deployData.csv", stringsAsFactors=F, header=T)
taxaList <- read.csv(file="taxonStageFeed.csv", stringsAsFactors=F, header=T)
LMreg <- read.csv(file="LMreg_hierarchy.csv", stringsAsFactors=F, header=T)

taxonData <- merge(taxaList, LMreg, by.x=c("Taxon","Stage"), by.y=c("Taxon","Life.stage"), all=T)
missingReg <- taxonData[is.na(taxonData$Class),] #identify the taxons and life stages that don't have associated length-weight regression parameters 
#write.table(test, file="hierQAQC", sep=",", col.names=T)

setwd("C:/Users/mark.henderson/Desktop/Sacramento River drift data/R analysis")
rawData <- read.csv(file="raw data 2012-2013.csv", stringsAsFactors=F, header=T)
rawData$Date<-as.Date(rawData$Date, "%m/%d/%Y")

allData <- merge(rawData, taxonData, by=c("Taxon", "Stage"))
#write.table(allData, file="driftData", sep=",", col.names=T)
allData$fullCount <- with(allData, Count*(1/prop.sampled))
allData$biomass <- with(allData, fullCount*a*Length.mm^b)
zoop <- c("Amphipoda", "Copepoda", "Cladocera", "Decapoda", "Ostracoda", "mysida")
allData$group[allData$Class=="Insecta"]<-"insect"
allData$group[allData$Analysis.taxon %in% zoop] <- "zoop"
allData$group[is.na(allData$group)]<-"other"
allData$group<-as.factor(allData$group)

#Summarize the biomass and count data by replicate
rep.summary <- merge(aggregate(biomass~Analysis.taxon+Site+Date+Replicate, data=allData, FUN=sum),
                      aggregate(fullCount~Analysis.taxon+Site+Date+Replicate, data=allData, FUN=sum))
meanLen.rep<- aggregate(Length.mm~Analysis.taxon+Site+Date+Replicate, data=allData, FUN=mean)
rep.summary<-merge(rep.summary, meanLen.rep, by=c("Analysis.taxon","Site","Date","Replicate"))

#summarize the biomass and count data by site and date

site.summary <- merge(aggregate(biomass~Analysis.taxon+Site+Date, data=allData, FUN=sum),
                aggregate(fullCount~Analysis.taxon+Site+Date, data=allData, FUN=sum))
meanLen.site<- aggregate(Length.mm~Analysis.taxon+Site+Date, data=allData, FUN=mean)
site.summary<-merge(site.summary, meanLen.site, by=c("Analysis.taxon","Site","Date"))

#look to see if there is any trend in the stages through time. Haven't looked at stages by analysis taxon
#but I should do that for the major taxa in the samples (Chironomid, baetid, etc.)
stageBio <- aggregate(biomass~Stage+Site+Date, data=allData, FUN=sum)
stageCount <- aggregate(Count*(1/prop.sampled)~Stage+Site+Date, data=allData, FUN=sum)

stageSummary <- merge(stageBio, stageCount, by=c("Stage","Site","Date"))
test<-stageSummary[with(stageSummary,order(Site,Date,Stage)),]
#write.table(test, file="stageTrend", sep=",", col.names=T)

#From the AqBiolDataAnalysis program - make some data summaries by different groups and origins to plot results
cols<-c("Analysis.taxon","group", "Origin")
lookup<-allData[!duplicated(allData$Analysis.taxon),cols]

### site summary ##
site.summary <- merge(site.summary,lookup)

site.abund<-aggregate(fullCount~Site+Date, data=site.summary, sum)
site.abund.grp<-aggregate(fullCount~Site+Date+group, data=site.summary, sum)

names(site.abund)[match("fullCount",names(site.abund))]<-"abund.tot"
names(site.abund.grp)[match("fullCount",names(site.abund.grp))]<-"abund.grp"

site.summary<-merge(site.summary,site.abund, by=c("Site", "Date"))
site.summary<-merge(site.summary,site.abund.grp, by=c("Site", "Date", "group"))

site.summary$taxon.prop<-with(site.summary, fullCount/abund.tot)
site.summary$grp.prop<-with(site.summary, abund.grp/abund.tot)

### replicate summary ##
drift.summary<-merge(rep.summary,lookup)

abund.tot<-aggregate(fullCount~Site+Date+Replicate, data=drift.summary, sum)
names(abund.tot)[match("fullCount",names(abund.tot))]<-"abund.tot"
drift.summary<-merge(drift.summary,abund.tot, by=c("Site", "Date", "Replicate"))
drift.summary$abund.prop<-with(drift.summary, fullCount/abund.tot)

#summary plots for total abundance
library(ggplot2)
qplot(x=Date,y=fullCount,data=abund.tot,colour=Site, ylim=c(0,5000))
boxplot(fullCount~Date,data=abund.tot, subset=Site=="Verona")
boxplot(fullCount~Date,data=abund.tot, subset=Site=="Colusa")
boxplot(fullCount~Date,data=abund.tot, subset=Site=="Hamilton City")
boxplot(fullCount~Date,data=abund.tot, subset=Site=="Mooney Island")
boxplot(fullCount~Date,data=abund.tot, subset=Site=="Danger Island")

###  Insect abundance  ###
insect.abund<-aggregate(fullCount~Site+Date+Replicate+Insect, data=drift.summary, sum)
insect.abund<-merge(insect.abund, abund.tot)
insect.abund<-insect.abund[which(insect.abund$Insect=='insect'),]
insect.abund$insect.prop<-insect.abund$fullCount/insect.abund$abund.tot

#summary plots for insect abundance
library(ggplot2)
qplot(x=Date,y=insect.prop,data=insect.abund,colour=Site)+geom_line()

### Aquatic abundance  ###
aquatic.abund<-aggregate(fullCount~Site+Date+Replicate+Origin, data=drift.summary, sum)
aquatic.abund<-merge(aquatic.abund, abund.tot)
aquatic.abund<-aquatic.abund[which(aquatic.abund$Origin=='Aquatic'),]
aquatic.abund$aquatic.prop<-aquatic.abund$fullCount/aquatic.abund$abund.tot

#summary plots for aquatic abundance
library(ggplot2)
qplot(x=Date,y=aquatic.prop,data=aquatic.abund,colour=Site)+geom_line()

### dominant taxon abundance  ###

keep <- c("Site", "Date", "Analysis.taxon", "group", "fullCount", "abund.prop")
dom.abund <- site.summary[with(site.summary,order(Site, Date, -abund.prop)),]

#summarize by the maximum value in each replicate
max.taxon <- aggregate(abund.prop~Site+Date, data=site.summary, FUN=max)
max.taxon <- merge(dom.abund,max.taxon)

#summarize by the taxon that have more than 10 percent of the catch
dom.abund <- dom.abund[which(dom.abund$abund.prop>0.1),keep]

write.table(max.taxon, file="maxTaxon", sep=",", col.names=T)

### dominant group abundance  ###

dom.abund.grp <- site.summary[with(site.summary,order(Site, Date, -grp.prop, -taxon.prop)),]
max.grp <- dom.abund.grp[!duplicated(dom.abund.grp[c("Site","Date")]),]

#summarize by the taxon that have more than 10 percent of the catch
dom.abund <- dom.abund[which(dom.abund$abund.prop>0.1),keep]

write.table(max.grp, file="maxGroup", sep=",", col.names=T)

######### Repeat analysis for biomass  ##################

biomass.tot<-aggregate(Biomass~Site+month, data=drift.summary, sum)
names(biomass.tot)[3]<-"biomass.tot"
drift.summary<-merge(drift.summary,biomass.tot, by=c("Site", "month"))
drift.summary$biomass.prop<-with(drift.summary, Biomass/biomass.tot)

insect.biomass<-aggregate(Biomass~Site+month+Insect, data=drift.summary, sum)
insect.biomass<-merge(insect.biomass, biomass.tot)
insect.biomass<-insect.biomass[which(insect.biomass$Insect=='insect'),]
insect.biomass$insect.prop<-insect.biomass$Biomass/insect.biomass$biomass.tot

aquatic.biomass<-aggregate(Biomass~Site+month+Origin, data=drift.summary, sum)
aquatic.biomass<-merge(aquatic.biomass, biomass.tot)
aquatic.biomass<-aquatic.biomass[which(aquatic.biomass$Origin=='Aquatic'),]
aquatic.biomass$aquatic.prop<-aquatic.biomass$Biomass/aquatic.biomass$biomass.tot

library(plyr)
dataSummary <-join_all(list(insect.abund, insect.biomass, aquatic.abund, aquatic.biomass), by=c("Site", "month"))

write.table(dataSummary, file="dataSummary", sep=",", col.names=T)

sample.cols <- c("Site", "month")
keep <- c("Site", "month", "summary.taxa", "Biomass", "biomass.prop")
dom.biomass <- drift.summary[with(drift.summary,order(Site, month, -biomass.prop)),]
dom.biomass <- dom.biomass[which(dom.biomass$biomass.prop>0.1), keep]

totAbund <- aggregate(Abundance~summary.taxa, data=drift.summary, sum)
totAbund <- totAbund[order(-totAbund$Abundance),]

totBiomass <- aggregate(Biomass~summary.taxa, data=drift.summary, sum)
totBiomass <- totBiomass[order(-totBiomass$Biomass),]

keep <- c("Chironomidae", "Baetidae", "Simuliidae", "Hydropsychidae", "Glossosomatidae", 
          "Annelida: Oligochaeta", "Crustacea: Copepoda", "Crustacea: Amphipoda", "Crustacea: Cladocera", "Pisces")
salmonPrey <- subset(drift.summary, subset=summary.taxa %in% keep)
write.table(salmonPrey, file="Chinook Diet Items", sep=",", col.names=T)
