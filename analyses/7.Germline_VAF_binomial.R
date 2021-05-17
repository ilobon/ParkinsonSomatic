options(stringsAsFactors = F)
library(ggplot2)
library(gridExtra)

#Load depth data from high-confidence germline heterozygous positions, defined as:
#-The five tissues have support for alt and ref
#-VAF between 0.45 and 0.55 in at least one tissue
#-In dbsnp
#-Total depth 20-100
#-DV2 was removed from this list
GermHetDP<-read.delim("data/HighConf_GermlineHet_depths.txt")


#### Binomial observed vs theoretical #####

#Bin total depth
HConfHetInfoBins<-cbind.data.frame(GermHetDP,
                                   binDP=cut(GermHetDP$DP,seq(20,100,20)))
levels(HConfHetInfoBins$binDP)<-c("20-40","40-60","60-80","80-100")

#Subsample bins to make them comparable
minN<-sort(table(HConfHetInfoBins$binDP))[1]
HConfHetInfoBinsSampled<-HConfHetInfoBins[c(sample(which(HConfHetInfoBins$binDP=="80-100"),minN),
                                            sample(which(HConfHetInfoBins$binDP=="20-40"),minN),
                                            sample(which(HConfHetInfoBins$binDP=="40-60"),minN),
                                            sample(which(HConfHetInfoBins$binDP=="60-80"),minN)),]

#Binomial background 
randomBinomial<-cbind.data.frame(VAF=rbinom(nrow(GermHetDP),100,0.5)/100)
randomBinomialSampled<-cbind.data.frame(VAF=randomBinomial$VAF[sample(1:nrow(randomBinomial),minN)])

#Main plot
p<-ggplot()+geom_histogram(data=randomBinomialSampled,aes(VAF),fill="grey",binwidth=0.01)+
  geom_histogram(data=HConfHetInfoBinsSampled,aes(AD/DP),fill="#3182bd",binwidth=0.01,alpha=0.7)+
  xlab("VAF")+ggtitle("High confidence heterozygous positions")+
  scale_x_continuous(limits=c(0,1),breaks=seq(0,1,0.1))+theme_bw()+
  facet_grid(binDP~.)

#Make legend with both colors
dummy_df<-cbind.data.frame(VAF=c(2,3,2,4),set=c("Binomial","Binomial","Empirical","Empirical"))

dummy_plot<-ggplot(dummy_df)+geom_histogram(aes(VAF,fill=set),binwidth=1)+
  scale_fill_manual(name="",values=c("grey","#3182bd"))

extract.legend<-function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  return(legend)} 

legend<-extract.legend(dummy_plot)

#Plot together
grid.arrange(p,legend,nrow=1,widths=c(0.8,0.2))



##### Power to detect germline het variants by n of tissues #####

#Binomial test pvalue for each variant and tissue
binomialTests<-lapply(unique(GermHetDP$tissue), function(t)
  sapply(round(GermHetDP$AD[GermHetDP$tissue==t]/GermHetDP$DP[GermHetDP$tissue==t]*100),
         function(x) binom.test(c(x,100-x))$p.value))

#Significant binomial with one tissue (for each tissue)
sb1<-sapply(1:5, function(i) sum(binomialTests[[i]]<0.05)/length(binomialTests[[i]]))*100
#Significant binomial with two tissues
sb2<-unlist(lapply(1:5, function(t1) unlist(lapply(c(1:5)[-t1], function(t2) 
  sum(binomialTests[[t1]]<0.05 & binomialTests[[t2]]<0.05)/length(binomialTests[[t1]])*100))))
#Significant binomial with 3 tissues
sb3<-unlist(lapply(1:5, function(t1) unlist(lapply(c(1:5)[-t1], function(t2) 
  unlist(lapply(c(1:5)[-c(t1,t2)], function(t3)
    sum(binomialTests[[t1]]<0.05 & binomialTests[[t2]]<0.05 & binomialTests[[t3]]<0.05)/length(binomialTests[[t1]])*100))))))
#Significant binomial with 4 tissues
sb4<-unlist(lapply(1:5, function(t1) unlist(lapply(c(1:5)[-t1], function(t2) 
  unlist(lapply(c(1:5)[-c(t1,t2)], function(t3) unlist(lapply(c(1:5)[-c(t1,t2,t3)], function(t4)
    sum(binomialTests[[t1]]<0.05 & binomialTests[[t2]]<0.05 & binomialTests[[t3]]<0.05 & binomialTests[[t4]]<0.05)/length(binomialTests[[t1]])*100))))))))

#Calculate median and sd 
medians<-c(median(sb1),median(sb2),median(sb3),median(sb4))
sds<-c(sd(sb1),sd(sb2),sd(sb3),sd(sb4))

#Build DF for plotting
propGermHetUndetected<-cbind.data.frame(n=1:4,Percentage=medians, SD=sds)
propGermHetUndetected$n<-factor(propGermHetUndetected$n,levels=propGermHetUndetected$n,ordered=T)

#Exponential model of median % undetected germ het vars vs n tissues sampled
model<-lm(log(medians) ~ c(1:4))
undetectedExponential<-exp(predict(model,list(1:4)))
propGermHetUndetectedModel<-cbind.data.frame(propGermHetUndetected,Exponential=undetectedExponential)

#Plot
ggplot(propGermHetUndetectedModel)+geom_bar(aes(n,Percentage,fill=n),stat="identity",width=0.95,alpha=0.9)+
  geom_point(aes(n,Percentage))+
  geom_errorbar(aes(x=n,ymin=Percentage-SD,ymax=Percentage+SD),width=.2,position=position_dodge(.9))+
  geom_point(aes(x=n,y=Exponential),color="#cb181d",size=1,shape=18)+
  geom_line(aes(x=n,y=Exponential,group=1),color="#cb181d")+
  xlab("Number of tissues")+ggtitle("Positions with significant binomial test")+
  scale_fill_manual(values=c("#6baed6","#4292c6","#2171b5","#084594"))+
  theme_bw()+
  theme(legend.position="none")





