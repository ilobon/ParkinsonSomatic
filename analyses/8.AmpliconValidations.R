options(stringsAsFactors=F)

library(tidyr)
library(ggplot2)
library(plyr)
library(scales)
library(ggpubr)
library(cowplot)

#Load candidates info and ultra deep amplicon seq counts
CandidatesInfo<-read.delim2("data/PD_SomaticCandidates_Info.txt")
readCountsAmplicons<-read.delim2("data/candidates_ampliconReadCounts.txt",header=F)

#### Amplicon coverage ####

#Order both by tier and position, remove positions from read counts
CandidatesInfo<-CandidatesInfo[order(CandidatesInfo$Tier,CandidatesInfo$Chr,CandidatesInfo$Pos),]
readCountsAmplicons<-readCountsAmplicons[match(CandidatesInfo$Pos,readCountsAmplicons$V2),-c(1,2)]

#Create tissue, inds, sample, nucleotide vectors
tissues<-c("B","C","E","N","S")
inds<-paste("DV",seq(1,10),sep="")
samples<-unlist(lapply(inds, function(x) paste(x,tissues,sep="")))
colnames(readCountsAmplicons)<-samples
nucleotides<-c("A","C","G","T")

#List per candidate, DF of sample vs nucleotide counts
readCountsList<-list()
for (row in 1:nrow(readCountsAmplicons)){
  df<-list()
  for (j in 1:ncol(readCountsAmplicons)){
    df[[j]]<-do.call(cbind.data.frame,lapply(strsplit(readCountsAmplicons[row,j],",")[[1]], function(x) 
      as.numeric(x)))
    colnames(df[[j]])<-nucleotides}
  readCountsList[[row]]<-do.call(rbind.data.frame,df)
  rownames(readCountsList[[row]])<-samples
}

#Calculate alt VAF and collapse in one DF
VAF<-do.call(cbind.data.frame,lapply(1:nrow(CandidatesInfo), function(variant)
  sapply(1:nrow(readCountsList[[variant]]), function(row) 
    readCountsList[[variant]][row,CandidatesInfo$Alt[variant]]/sum(readCountsList[[variant]][row,])*100)))
colnames(VAF)<-1:ncol(VAF)
VAF<-cbind.data.frame(gather(VAF,Variant,VAF),Sample=factor(samples,levels=samples,ordered=T))

#Add variant called for sample (Y/N) in exome data
calledInfo<-list()
for (variant in 1:nrow(CandidatesInfo)){
called<-rep("N",50)
called[which(samples%in%strsplit(CandidatesInfo$Samples[variant],",")[[1]])]<-"Y"
calledInfo[[variant]]<-called
}
VAF<-cbind.data.frame(VAF,Called=unlist(calledInfo))

#Add total coverage
VAF<-cbind.data.frame(VAF, Coverage=unlist(lapply(readCountsList, function(x) apply(x,1,sum))))

#Add individual and tissue 
VAF<-cbind.data.frame(VAF,Individual=sapply(VAF$Sample, function(x) substr(x,1,nchar(as.character(x))-1)),
                      Tissue=sapply(VAF$Sample, function(x) 
                        substr(x,nchar(as.character(x)),nchar(as.character(x)))))
VAF$Individual<-factor(VAF$Individual,levels=inds,ordered=T)

#Plot
ggplot(VAF)+geom_histogram(aes(log10(Coverage),fill=Called),bins=50)+
  scale_fill_manual(values=c("#bababa","#f46d43"))+
  facet_grid(Individual~Tissue)+
  theme_bw()

#### Amplicon data per variant ####

#Add rank of alt allele (second most supported, third...)
AltPosition<-unlist(lapply(1:length(readCountsList), function(i)
  sapply(1:nrow(readCountsList[[i]]), function(row)
    which(names(sort(readCountsList[[i]][row,],decreasing=T))==CandidatesInfo$Alt[i]))))
VAF<-cbind.data.frame(VAF,AltPos=AltPosition)

#VAF and coverage in same variable to plot in Y as positive and negative values
VAFgathered<-gather(VAF,Data,Y,-Variant,-Sample,-Called,-Individual,-Tissue,-AltPos)
VAFgathered$Data<-factor(VAFgathered$Data,levels=c("VAF","Coverage"),ordered=T)
VAFgathered$Y<-ifelse(VAFgathered$Data=="VAF",as.numeric(VAFgathered$Y),-(as.numeric(VAFgathered$Y)))

#Make a colour variable, Called Y/N for VAF and AltPos for Coverage rows 
VAFgathered<-cbind.data.frame(VAFgathered,ColourVar=ifelse(VAFgathered$Data=="VAF",
                                                           VAFgathered$Called,VAFgathered$AltPos))

#### Amplicon data per variant plots ####

#Bar colours
colors=c("1"="#762a83","2"="#006837","3"="#66bd63","4"="#a6d96a",N="#fee08b",Y="#f46d43")

#Prepare legends
dfLegend<-cbind.data.frame(Exome=c("Called","Called","Called","NotCalled"),
                         AltAllele=c("1st","2nd","3rd","4th"),
                         n=1:4,
                         cov=rep("MeanCov+1SD",4))

leg1<-cowplot::get_legend(ggplot(dfLegend)+geom_bar(aes(Exome,n,fill=Exome),stat="identity")+
                            scale_fill_manual(values=unname(colors[c(6,5)])))

leg2<-cowplot::get_legend(ggplot(dfLegend)+geom_bar(aes(AltAllele,n,fill=AltAllele),stat="identity")+
  scale_fill_manual(values=unname(colors[1:4]))+
  geom_line(aes(AltAllele,n,group=cov,color=cov),linetype="dashed")+
  scale_color_manual(values="#525252",name="")+
  theme_bw())

#Open file
#pdf("SuppFile_Amplicons.pdf",10,4)

#Horizontal coverage threshold line value
hlineDF<-cbind.data.frame(Data=factor("Coverage",levels=c("VAF","Coverage"),ordered=T),yint=0)

for (var in c(2:59)){
  print(var)
  #Number of breaks in Y axis
  nbreaks=4
  #Coverage
  maxCOVbreak=round_any(min(VAFgathered$Y[which(VAFgathered$Variant==var & VAFgathered$Data=="Coverage")]),10000)
  COVbreaks=seq(maxCOVbreak,0,
                round_any(abs(maxCOVbreak)/nbreaks,10^round_any(log10(abs(maxCOVbreak)/nbreaks),1,f=floor)))
  #VAF
  maxVAF=max(VAFgathered$Y[which(VAFgathered$Variant==var & VAFgathered$Data=="VAF")],na.rm=T)
  maxVAFbreak=round_any(maxVAF,10^round_any(log10(maxVAF),1,f=floor))
  VAFbreaks=seq(0,maxVAFbreak,round_any(maxVAFbreak/nbreaks,10^round_any(log10(maxVAFbreak/nbreaks),1,f=floor)))
  ybreaks=unique(c(COVbreaks,VAFbreaks))
  ylabs=unique(c(comma(abs(COVbreaks)),VAFbreaks))
  #yintercept
  hlineDF$yint=mean(VAFgathered$Y[which(VAFgathered$Variant==var & VAFgathered$Data=="Coverage")])+
    sd(VAFgathered$Y[which(VAFgathered$Variant==var & VAFgathered$Data=="Coverage")])
  
  plot<-ggplot(VAFgathered[which(VAFgathered$Variant==var),])+
    geom_bar(aes(Sample,Y,fill=ColourVar),stat="identity")+
    geom_hline(data=hlineDF,aes(yintercept=yint),color="#525252",linetype="dashed")+
    scale_y_continuous(expand=c(0,0),breaks=ybreaks,labels=ylabs)+
    scale_fill_manual(values=colors[sort(unique(VAFgathered$ColourVar[which(VAFgathered$Variant==var)]))])+
    ggtitle(paste(CandidatesInfo$Gene[var]," - Tier ",CandidatesInfo$Tier[var],sep=""))+
    ylab("")+
    facet_grid(Data~.,scales="free_y")+
    theme(axis.text.x=element_text(angle=90,size=9,vjust=0.5,hjust=1),
          panel.spacing.y=unit(0,"lines"),
          panel.background = element_blank(),
          strip.background = element_rect(color="#bababa"),
          axis.ticks.x=element_blank(),
          legend.position="none")
  
  print(plot_grid(plot,
            plot_grid(leg1,leg2,NULL,nrow=3,align = "v",rel_heights = c(0.35,0.4,0.15)),
            ncol=2,rel_widths = c(0.85,0.15)))
  
}
#dev.off()


#### Validation decisions #####

decision<-list()
validatedSamples<-list()

for (var in 1:59){
  print(var)
  
  #Called individual
  calledInd=unique(as.character(VAF$Individual[VAF$Variant==var & VAF$Called=="Y"]))
  #Coverage threshold for position (mean-sd)
  coverageLim=mean(VAF$Coverage[VAF$Variant==var])-sd(VAF$Coverage[VAF$Variant==var])
  #Background samples VAFs
  backgVAF=VAF$VAF[VAF$Variant==var & VAF$Individual!=calledInd & VAF$Coverage>coverageLim]
  
  ##Called tissues info
  #Alt ranking
  calledAltRanking=VAF$AltPos[VAF$Variant==var & VAF$Called=="Y"]
  #VAFs
  calledVAF=VAF$VAF[VAF$Variant==var & VAF$Called=="Y"]
  #Passing VAF and Alt rank
  calledPassVAFAltrank<-(calledVAF>mean(backgVAF)+2*sd(backgVAF) & calledAltRanking==2)
  calledTissuesPassVAFAltrank<-as.character(VAF$Sample[VAF$Variant==var & VAF$Called=="Y"][calledPassVAFAltrank])
  
  ##Not called tissues
  #Pass alt ranking and coverage limit
  notcalledPassCovAltrank<-as.character(VAF$Sample[VAF$Variant==var & VAF$Individual==calledInd & 
                                                     VAF$Called=="N" & VAF$Coverage>coverageLim & VAF$AltPos==2])
  #Pass VAF
  notcalledPassVAF<-VAF$VAF[VAF$Variant==var & VAF$Sample%in%notcalledPassCovAltrank]>mean(backgVAF)+2*sd(backgVAF)
  #All VAFs >30?
  ghet<-(length(notcalledPassVAF)>0 & all(VAF$VAF[VAF$Variant==var & VAF$Sample%in%notcalledPassCovAltrank]>30))
  
  #Germline heterozygous
  if(any(calledAltRanking==1) | ghet==T ){
    decision[[var]]<-"GHET"
    validatedSamples[[var]]<-NA
  }
  
  #False positives
  if( (length(calledTissuesPassVAFAltrank)==0 | all(is.na(calledVAF))) & 
      (all(VAF$AltPos[VAF$Variant==var & VAF$Individual==calledInd & VAF$Coverage>coverageLim]!=1))){
    decision[[var]]<-"FP"
    validatedSamples[[var]]<-NA
  }
  
  #Tissue exclusive
  if(length(calledTissuesPassVAFAltrank)>0 & length(notcalledPassVAF)==0){
    decision[[var]]<-"TissueExclusive"
    validatedSamples[[var]]<-as.character(VAF$Sample[VAF$Variant==var & VAF$Called=="Y"])
  }
  
  #Multi tissue
  if(length(calledTissuesPassVAFAltrank)>0 & sum(notcalledPassVAF)>0 & ghet==F){
    decision[[var]]<-"MultiTissue"
    validatedSamples[[var]]<-paste(sort(c(calledTissuesPassVAFAltrank,notcalledPassCovAltrank[notcalledPassVAF])),
                                   collapse=",")
  }
}

decision<-unlist(decision)
validatedSamples<-unlist(validatedSamples)


#### Save amplicon VAFs ####

ampliconValidatedVAFs<-list()
for (var in 1:59){
    dfCounts<-readCountsList[[var]][rownames(readCountsList[[var]])%in%paste(CandidatesInfo$Individual[var],
                                                                             tissues,sep=""),]
    ampliconValidatedVAFs[[var]]<-unname(dfCounts[,CandidatesInfo$Alt[var]]/apply(dfCounts,1,sum))
  }

ampliconValidatedVAFs<-do.call(rbind.data.frame,ampliconValidatedVAFs)
colnames(ampliconValidatedVAFs)<-tissues


CandidatesInfoValidations<-cbind.data.frame(CandidatesInfo,
                                    Decision=decision, 
                                    ValidatedSamples=validatedSamples,
                                    ampliconValidatedVAFs)

write.table(CandidatesInfoValidations,"data/CandidatesInfoValidations.txt",
            row.names=F,quote=F,sep="\t")



