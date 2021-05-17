options(stringsAsFactors = F)

library(ggplot2)
library(tidyr)
library(lsa)

`%notin%`<-Negate(`%in%`)

#### Input data ####

infoValidations<-read.delim2("data/CandidatesInfoValidations.txt")

#Validated in at least one CNS sample

valCNS<-which(sapply(infoValidations$ValidatedSamples, 
                     function(x) length(strsplit(x,",")[[1]]))>1 | 
                grepl("C|E|N|S",infoValidations$ValidatedSamples))

infoValidationsSubset<-infoValidations[valCNS,]


#### Format data ####

#Collapse substitution and context

nucleotidePairs<-cbind.data.frame(a=c("A","C","G","T"),b=c("T","G","C","A"))

collapseSubs<-function(ref,alt){
  if(ref%in%c("C","T")){
    collapsedSubstitution<-paste(ref,alt,sep=">")
  }else{
    revref<-nucleotidePairs$b[nucleotidePairs$a==ref]
    revalt<-nucleotidePairs$b[nucleotidePairs$a==alt]
    collapsedSubstitution<-paste(revref,revalt,sep=">")
  }
  return(collapsedSubstitution)
}

collapseContext<-function(context){
  if(strsplit(context,"")[[1]][2]%in%c("C","T")){
    collapsedContext<-context
  }else{
    collapsedContext<-paste(sapply(strsplit(context,"")[[1]], function(x) 
      nucleotidePairs$b[nucleotidePairs$a==x]),collapse="")
  }
  return(collapsedContext)
}

collapsedSubstitutions<-sapply(1:nrow(infoValidationsSubset), function(i) 
  collapseSubs(infoValidationsSubset$Ref[i],infoValidationsSubset$Alt[i]))

collapsedContext<-sapply(infoValidationsSubset$Context, function(x) collapseContext(x))

strand<-sapply(1:nrow(infoValidationsSubset), function(i) 
  if(strsplit(collapsedSubstitutions[i],"")[[1]][1]==infoValidationsSubset$Ref[i]){
    if(infoValidationsSubset$Strand[i]==1){"Transcribed"}else{"Untranscribed"}}else{
      if(infoValidationsSubset$Strand[i]==1){"Untranscribed"}else{"Transcribed"}})


#Create empty DF for spectrum
nucleotides<-c("A","C","G","T")
contextC<-data.frame(lapply(expand.grid(nucleotides,"C",nucleotides),as.character), stringsAsFactors=FALSE)
contextC<-sapply(1:nrow(contextC), function(i) paste0(as.character(contextC[i,]),collapse=""))
contextT<-data.frame(lapply(expand.grid(nucleotides,"T",nucleotides),as.character), stringsAsFactors=FALSE)
contextT<-sapply(1:nrow(contextT), function(i) paste0(as.character(contextT[i,]),collapse=""))
substitutions<-c("C>A","C>G","C>T","T>A","T>C","T>G")
spectrumDFempty<-rbind.data.frame(expand.grid(substitutions[1:3],contextC),expand.grid(substitutions[4:6],contextT))
spectrumDFempty<-data.frame(lapply(spectrumDFempty,as.character), stringsAsFactors=FALSE)
colnames(spectrumDFempty)<-c("Substitution","Context")

spectrumDF<-cbind.data.frame(spectrumDFempty, n=sapply(1:nrow(spectrumDFempty), function(i) 
  sum(collapsedSubstitutions==spectrumDFempty[i,1] & collapsedContext==spectrumDFempty[i,2])))

#### Spectrum plot ####

pal<-c("#14BCF0","#000000","#D43431","#999999","#A0CF62","#ECC5C4")

p<-ggplot(spectrumDF)+geom_bar(aes(Context,n/sum(n),fill=Substitution),stat="identity",position="dodge")+
  scale_fill_manual(values=pal)+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))+
  facet_grid(.~Substitution,scales="free_x")+
  xlab("")+ylab("")+
  theme_bw()+ theme(axis.text.x=element_text(angle=90,size=6,margin=margin(-2,0,0,0)),
                    axis.ticks.x=element_blank(),
                    panel.grid=element_blank(),
                    legend.position="none",
                    strip.text=element_text(colour='white'))

g<-ggplot_gtable(ggplot_build(p))
strips<-which(grepl('strip-', g$layout$name))

for (i in 1:length(strips)) {
  k<-which(grepl('rect', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
  g$grobs[[strips[i]]]$grobs[[1]]$children[[k]]$gp$fill <- pal[i]
  g$grobs[[strips[i]]]$grobs[[1]]$children[[k]]$gp$col <- pal[i]
}

mutationalSpectrum<-g


#### Correlation with COSMIC SBS ####

SBS_signatures_2019_05_22<-read.csv("data/sigProfiler_SBS_signatures_2019_05_22.csv")
sequencingArtSBS<-c("SBS27","SBS43","SBS45","SBS46","SBS47","SBS48","SBS49",
                    "SBS50","SBS51","SBS52","SBS53","SBS54","SBS55","SBS56",
                    "SBS57","SBS58","SBS59","SBS60")

SBSandCounts<-merge(SBS_signatures_2019_05_22,spectrumDF,by.x=c(1,2),by.y=c(1,2))

corr<-cbind.data.frame(SBS=colnames(SBSandCounts)[3:69],
                       Pcor=sapply(3:69, function(i) 
                         cor(SBSandCounts[,i], SBSandCounts$n)))

corr<-cbind.data.frame(corr,Set=rep("Other",nrow(corr)))
corr$Set[corr$SBS%in%sequencingArtSBS]<-"SequencingArtifact"
corr<-corr[order(corr$Pcor),]
corr$SBS<-factor(corr$SBS,levels=corr$SBS,ordered=T)


SBS15_text<-"Defective DNA mismatch repair"
SBS5_text<-"Most tumours (Alexandrov 2020); 
De novo muts (Rahbari 2016);
Somatic muts (Bae 2018)"
SBS6_text<-"Defective DNA mismatch repair; 
De novo mutations (Rahbari 2016)"
SBS1_text<-"Endogenous 5-methylcytosine 
deamination"

SBStext<-ls()[grep("_text",ls())]

textAnnotations<-cbind.data.frame(SBS=sapply(strsplit(SBStext,"_"), function(x) x[1]),
                                  Text=sapply(SBStext, function(x) eval(parse(text=x))))


corrPlot<-ggplot(corr[corr$Pcor>0.3,])+geom_bar(aes(SBS,Pcor),stat="identity")+
  geom_text(data=textAnnotations,aes(SBS,0.01,label=Text),size=3,color="white",hjust=0)+
  coord_flip()+
  xlab("")+ ylab("Pearson correlation")+
  scale_y_continuous(breaks = c(0,0.2,0.4))+
  theme_bw()+
  theme(axis.text.y=element_text(size=13),
        axis.text.x=element_text(size=12),
        axis.title.x=element_text(size=13),
        panel.grid.major.x=element_blank(),
        panel.grid.minor.x=element_blank(),
        panel.grid.major.y=element_blank(),
        panel.grid.minor.y=element_blank())



