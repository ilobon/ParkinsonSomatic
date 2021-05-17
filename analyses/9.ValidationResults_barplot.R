options(stringsAsFactors = F)
library(ggplot2)

InfoValidations<-read.delim2("data/CandidatesInfoValidations.txt")

#### Histogram validations vs exome ####

#Format data for plot
InfoValidations$Decision<-factor(InfoValidations$Decision,levels=unique(InfoValidations$Decision)[c(2,3,1,4)],ordered=T)
InfoValidations<-cbind.data.frame(InfoValidations,ExomeTissue=ifelse(InfoValidations$nTissues==1,InfoValidations$Tissues,"Multiple"))
InfoValidations$ExomeTissue<-factor(InfoValidations$ExomeTissue,levels=c("B","C","N","E","S","Multiple"),ordered=T)
levels(InfoValidations$ExomeTissue)<-c("Blood","Cerebellum","Neocortex","Striatum","S.nigra","Multiple")
#Order tissues by mean of calls per tissue and tier
levelsOrdered<-names(sort(apply(table(InfoValidations[,c("ExomeTissue","Tier")])[,-1],1,mean)))
#But multiple last
levelsOrdered<-levelsOrdered[c(which(levelsOrdered=="Multiple"),1:(which(levelsOrdered=="Multiple")-1),(which(levelsOrdered=="Multiple")+1):6)]
InfoValidations$ExomeTissue<-factor(InfoValidations$ExomeTissue,levels=levelsOrdered,ordered=T)
#Format decision labels
decisionLabels<-c("Pass (Multiple tissues)","Pass (Tissue exclusive)","False positive","Germline heterozygous")
names(decisionLabels)<-levels(InfoValidations$Decision)

tissueCols<-c("#bababa", "#04935c", "#8754a1", "#3a75c4", "#ff9b38","#e23b3d")

decisionTierPlot<-ggplot(InfoValidations[-which(InfoValidations$Tier==0),])+geom_bar(aes(Tier,fill=ExomeTissue),stat="count")+
  scale_fill_manual(values=tissueCols,name="Exome call")+
  ylab("Number of candidate sSNVs")+
  xlab("Tier")+
  facet_grid(.~Decision,labeller=labeller(Decision=decisionLabels))+
  theme_bw()+
  theme(strip.background=element_rect(colour="#bdbdbd",fill="#d9d9d9"),
        panel.grid.major.x=element_blank(),
        panel.grid.minor.x=element_blank(),
        panel.grid.major.y=element_blank(),
        panel.grid.minor.y=element_blank())



