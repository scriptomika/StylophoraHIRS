library(ape); library(phytools); library(vegan)

		# With big tree (w/ all clusters representative and novo oligotypes)...
		bigtre<-read.tree("RAxML_bestTree.allseqs")
		plot(bigtre,cex=0.3, no.margin=T)
		plot(ladderize(bigtre),cex=0.25, no.margin=T)
		nodelabels(cex=0.5)
		# ... pull out list of cluster IV/chlorophyllide reductases clade for now..
		drop<-extract.clade(tre, 316); tre$tip.label[!(tre$tip.label %in% drop$tip.label)]

#BLAST 'good' (non-chlorophyllides/ferrodoxins) against nr BLAST (culturable only)
#re-align & re-tree

tre<-ladderize(read.tree("RAxML_bipartitions.bootslabel"))
bs<-read.tree("allboots")
supp<-prop.clades(tre, bs)/length(bs)
tre$node.label<-supp
plot(tre,cex=0.3, no.margin=T)
nodelabels(cex=0.5)
tre<-reroot(tre,201)
plot(tre,cex=0.5, no.margin=T)


p<-numeric(length(tre$node.label))
p[tre$node.label >= 0.9] <- 1.5
p[tre$node.label < 0.9 & tre$node.label >= 0.8] <- 1
p[tre$node.label < 0.8 & tre$node.label >= 0.7] <- 0.75
p[tre$node.label < 0.7] <- NA
nodelabels(pch=16, cex = p, col="black")


#sample metadata
meta<-read.table("MEDmetafile.txt",row.names=1,header=T)
meta$Light<-factor(gsub("_light.*","",meta$Light,perl=T))

#read in count data
dat<-read.table("MATRIX-COUNT.txt",header=T,row.names=1)
colnames(dat)<-gsub("X00000","nifH_novo_",colnames(dat),perl=T)


#read in CDHIT clustering to combine counts for OTUs with 95% similarity
cd<-read.table("cdhit_clusters.txt",header=T,stringsAsFactors=F)
orddat<-dat[,cd$nifHtypeID]
aggregate(t(orddat)[,1],by=list(cd$clusterRepID),FUN=sum)

agg<-aggregate(t(orddat),by=list(cd$clusterRepID),FUN=sum); rownames(agg)<-agg$Group.1;
aggdat<-t(agg[,-1])

#summary table and figure for oligotyping
summary<-cbind(rowSums(aggdat[,!(colnames(aggdat)%in% c("nifH_novo_3411","nifH_novo_3412"))]),rowSums(aggdat[,c("nifH_novo_3411","nifH_novo_3412")]),rowSums(dat[,!(colnames(dat) %in% cd$nifHtypeID)]),rowSums(dat))
colnames(summary)<-c("cluster I","cluster III","non nifH","total reads")
#par(mfrow=c(2,1))
barplot(t(summary[,1:3]),legend=T,las=3,cex.names=0.5,ylab="Number of MiSeq reads")
barplot(t(summary[,1:3]/summary[,4]),las=3,cex.names=0.5,ylab="Relative proportion of MiSeq reads")

#does seawater have less than coral? proportion wise
x<-rowSums(summary[,1:2]/summary[,4]); barplot(x)
mean(x[meta$host!="none"])
mean(x[meta$host=="none"])
wilcox.test(x[meta$host!="none"],x[meta$host=="none"])


dat<-data.frame(aggdat)
datsq<-sqrt(dat)
rowSums(datsq)
datsq<-datsq[rowSums(datsq)>1,]

rownames(datsq)[which(rownames(datsq) %in% rownames(meta)==F)] # all libraries match metadata
datsq <-datsq[which(rownames(datsq) %in% rownames(meta)==T),]
rownames(meta)[which(rownames(meta)%in% rownames(datsq)  ==F)] 	
meta<-meta[which(rownames(meta)%in% rownames(datsq)  ==T),]

#transform to relative abundance
relabun<-datsq/rowSums(datsq,na.rm=T)
#write.csv(relabun,file="cdhit.nifH.relabun.csv")
adonis2(relabun ~ host, data=meta,permutations=999)

         # Df SumOfSqs     F Pr(>F)  
# host      1   0.7509 2.443  0.012 *
# Residual 29   8.9135               
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

hostdepthlight<-factor(meta$host:meta$Depth:meta$Light)
adonis2(relabun ~ hostdepthlight, permutations=999)


adonis2(relabun ~ host+Depth+Light+ExptRun, data=meta,permutations=999)
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999

# adonis2(formula = relabun ~ host + Depth + Light + ExptRun, data = meta, permutations = 999)
         # Df SumOfSqs      F Pr(>F)   
# host      1   0.7509 2.6414  0.007 **
# Depth     1   0.2252 0.7923  0.660   
# Light     1   0.9420 3.3136  0.002 **
# ExptRun   2   0.6394 1.1247  0.312   
# Residual 25   7.1068                   

adonis2(relabun ~ host+Light, data=meta,permutations=999)

myLightpvals<-c()
mydepthpvals<-c()
myhostpvals<-c()
tukeyp<-c()
for (c in 1:ncol(relabun)){
	aov<-aov(relabun[,c] ~meta$host +meta$Light+ meta$Depth + meta$ExptRun)
	tukeyp[c]<-(TukeyHSD(aov)[[1]][[4]])
	rawphost<-summary(aov(relabun[,c] ~meta$host +meta$Light+ meta$Depth + meta$ExptRun))[[1]][["Pr(>F)"]][[1]]
	rawpLight<-summary(aov(relabun[,c] ~meta$host +meta$Light+ meta$Depth + meta$ExptRun))[[1]][["Pr(>F)"]][[2]]
	rawpdepth<-summary(aov(relabun[,c] ~meta$host +meta$Light+ meta$Depth + meta$ExptRun))[[1]][["Pr(>F)"]][[3]]
	myhostpvals <-c(myhostpvals, rawphost);myLightpvals <-c(myLightpvals, rawpLight);	mydepthpvals <-c(mydepthpvals, rawpdepth)

}
for (c in 1:ncol(relabun)){
	print(colnames(relabun)[c])
	print(summary(aov(relabun[,c] ~meta$host +meta$Light+ meta$Depth + meta$ExptRun))[[1]])

}



names(myhostpvals)<-names(relabun)
p.adjust(myhostpvals,method="bonferroni")
names(myLightpvals)<-names(relabun)
names(mydepthpvals)<-names(relabun)
mypvals<-data.frame(cbind(myhostpvals ,myLightpvals, mydepthpvals))
light_otus<-mypvals[mypvals$myLightpvals<0.05,]
depth_otus<-mypvals[mypvals$mydepthpvals<0.05,]
host_otus<-mypvals[mypvals$myhostpvals<0.05,]

#compare t-test to anova: habitat effect on phylotype abundance:
myt<-c();myp<-c(); df<-c()
for (c in 1:ncol(relabun)){
	
	myt[c]<-t.test(relabun[,c] ~meta$host)$statistic
	df[c]<-t.test(relabun[,c] ~meta$host)$parameter
	myp[c]<-t.test(relabun[,c] ~meta$host)$p.value
}
ttest<-data.frame(cbind(myt,df,myp))
rownames(ttest)<-colnames(relabun)
pcorrected <-p.adjust(ttest$myp,method="bonferroni")
cbind(ttest,pcorrected)



sig_otus<-unique(c(rownames(light_otus),rownames(depth_otus),rownames(host_otus)))
summary(aov(relabun[,"nifH_novo_0006"] ~ meta$host +meta$Light+ meta$Depth + meta$ExptRun))
 aggregate(relabun[,"nifH_novo_0006"],by=list(meta$host),FUN=mean)
# 1                  none 0.31241950
# 2 Stylophora_pistillata 0.02415345
 
summary(aov(relabun[,"nifH_novo_3411"] ~ meta$host +meta$Light+ meta$Depth + meta$ExptRun))
 aggregate(relabun[,"nifH_novo_3411"],by=list(meta$host),FUN=mean)
# 1                  none 0.15606315
# 2 Stylophora_pistillata 0.02815804


#repeat PERMANOVA without seawater
relacoral<-relabun[meta$host!="none",]; metacoral<-droplevels(meta[meta$host!="none",])
adonis2(relacoral ~ Depth+Light+ExptRun, data=metacoral,permutations=1000)
         # Df SumOfSqs      F Pr(>F)    
# Depth     1   0.2252 0.7090  0.757    
# Light     1   0.9420 2.9652  0.001 ***
# ExptRun   2   0.6394 1.0064  0.430    
# Residual 21   6.6711                  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

adonis2(relacoral ~ Light, data=metacoral,permutations=1000)
adonis2(relacoral ~ depthlight, data=metacoral,permutations=1000)

#pairwise comparison of treatments:

depthlight<-metacoral$Depth:metacoral$Light
#15: low v hi
adonis2(relacoral[depthlight=="Depth_15-18m:Low"|depthlight=="Depth_15-18m:High",]~ Light+ExptRun, data=metacoral[depthlight=="Depth_15-18m:Low"|depthlight=="Depth_15-18m:High",],permutations=1000)
mean(vegdist(relacoral[depthlight=="Depth_15-18m:High"|depthlight=="Depth_15-18m:Low",]))

#low: 15 v 5
adonis2(relacoral[depthlight=="Depth_5m:Low"|depthlight=="Depth_15-18m:Low",]~ Depth+ExptRun, data=metacoral[depthlight=="Depth_5m:Low"|depthlight=="Depth_15-18m:Low",],permutations=1000)
mean(vegdist(relacoral[depthlight=="Depth_5m:Low"|depthlight=="Depth_15-18m:Low",]))


#lowdeep v high 5
adonis2(relacoral[depthlight=="Depth_15-18m:Low"|depthlight=="Depth_5m:High",]~ Depth+Light+ExptRun, data=metacoral[depthlight=="Depth_15-18m:Low"|depthlight=="Depth_5m:High",],permutations=1000)
mean(vegdist(relacoral[depthlight=="Depth_5m:High"|depthlight=="Depth_15-18m:Low",]))



#highdeep v low5
adonis2(relacoral[depthlight=="Depth_5m:Low"|depthlight=="Depth_15-18m:High",]~ Depth+Light+ExptRun, data=metacoral[depthlight=="Depth_5m:Low"|depthlight=="Depth_15-18m:High",],permutations=1000)
mean(vegdist(relacoral[depthlight=="Depth_5m:Low"|depthlight=="Depth_15-18m:High",]))


#high15 v high 5
adonis2(relacoral[depthlight=="Depth_5m:High"|depthlight=="Depth_15-18m:High",]~ Depth+ExptRun, data=metacoral[depthlight=="Depth_5m:High"|depthlight=="Depth_15-18m:High",],permutations=1000)
mean(vegdist(relacoral[depthlight=="Depth_5m:High"|depthlight=="Depth_15-18m:High",]))


#5: low vs high
adonis2(relacoral[depthlight=="Depth_5m:Low"|depthlight=="Depth_5m:High",]~ Light+ExptRun, data=metacoral[depthlight=="Depth_5m:Low"|depthlight=="Depth_5m:High",],permutations=1000)

mean(vegdist(relacoral[depthlight=="Depth_5m:Low"|depthlight=="Depth_5m:High",]))



mean(vegdist(relacoral[depthlight=="Depth_15-18m:High",]))
mean(vegdist(relacoral[depthlight=="Depth_15-18m:Low",]))
mean(vegdist(relacoral[depthlight=="Depth_5m:High",]))
mean(vegdist(relacoral[depthlight=="Depth_5m:Low",]))


#repeat ANOVAs excluding seawater, to find signif nifH OTUs
myLightpvals<-c(); mydepthpvals<-c()

for (c in 1:ncol(relacoral)){
	rawpLight<-summary(aov(relacoral[,c] ~ Light+ Depth, data=metacoral))[[1]][["Pr(>F)"]][[1]]
	rawpdepth<-summary(aov(relacoral[,c] ~ Light+ Depth, data=metacoral))[[1]][["Pr(>F)"]][[2]]
	myLightpvals <-c(myLightpvals, rawpLight);	mydepthpvals <-c(mydepthpvals, rawpdepth)

}

names(myLightpvals)<-names(relabun)
names(mydepthpvals)<-names(relabun)
mypvals<-data.frame(cbind(p.adjust(myLightpvals,method='fdr'), p.adjust(mydepthpvals,method="fdr")))
light_otus<-mypvals[mypvals$myLightpvals<0.05,]
depth_otus<-mypvals[mypvals$mydepthpvals<0.05,]
sig_otus<-unique(c(sig_otus,c(rownames(light_otus),rownames(depth_otus),rownames(host_otus))))

aggregate(relacoral[,"nifH_novo_0153"],by=list(metacoral $Light),FUN=mean)
# 1    High 0.11396974
# 2     Low 0.01276101
aggregate(relacoral[,"nifH_novo_7096"],by=list(metacoral $Light),FUN=mean)
# 1    High 0.02789451
# 2     Low 0.20686222
aggregate(relacoral[,"nifH_novo_7457"],by=list(metacoral $Light),FUN=mean)
# 1    High 0.24121529
# 2     Low 0.04416233

summary(aov(relacoral[,"nifH_novo_7457"] ~ metacoral $Light+ metacoral $Depth + metacoral $ExptRun))

factors<-factor(meta$Light:meta$Depth)
anosim(relabun, factors)
# ANOSIM statistic R: 0.01696 
#      Significance:  0.314 
anosim(relabun, meta$host)
# ANOSIM statistic R: -0.2623 
      # Significance: 0.947 

 library(colorspace)
 cols<-rev(heat_hcl(12, h = c(0, -100), l = c(75, 40), c = c(40, 80), power = 1))
 cols<-diverge_hcl(12,c=100,l=c(50,90),power=1)
library(gplots)
plot(tre,cex=0.3, no.margin=T); nodelabels(pch=20,cex=supp)

lmat=rbind(4:3,2:1) #= 1-heatmap, 2-row dendrogram, 3-column dendrogram, and 4-key)

labels<-meta$host:meta$Light:meta$Depth

sample.order<-order(labels)

heatmap.2(t(as.matrix(relabun[sample.order,])), scale="none",trace="none", Colv=F,col=cols, dendrogram="row",margins=c(6,5),cexRow=0.5,key=T,lmat=lmat, lhei=c(1,30), lwid=c(1,10))

#to get the key
heatmap.2(t(as.matrix(relabun[sample.order,])), scale="none",trace="none", Colv=F,col=cols, dendrogram="row",margins=c(6,5),cexRow=0.5,key=T)

