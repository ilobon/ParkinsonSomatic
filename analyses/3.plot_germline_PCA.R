options(stringsAsFactors=F)
library(ggplot2)
library(gridExtra)

#### 1. Import and format data ####

#Import smartpca output
DV.p2 <- read.table("data/DV.p2.hardFiltered.eigenvecs", quote="\"", comment.char="")
#Add colnames
colnames(DV.p2)<-c("Sample",paste0("PC",1:10),"Individual")
#New column with tissue code
DV.p2<-cbind.data.frame(DV.p2,Tissue=unname(sapply(DV.p2$Sample, function(i) substr(i,nchar(i),nchar(i)))))
#Convert Individual to factor
DV.p2$Individual<-factor(DV.p2$Individual,levels=unique(DV.p2$Individual)[c(2:10,1)],ordered=T)

#### 2. PC1/PC2 and PC3/PC4 plots ####

#Plot DV2C on top so that is visible
pc12<-ggplot(DV.p2[DV.p2$Sample!="DV2C",])+geom_point(aes(PC1,PC2,color=Individual, shape=Tissue),size=8)+
  geom_point(data=DV.p2[DV.p2$Sample=="DV2C",],aes(PC1,PC2,shape=Tissue,color=Individual),size=8)+
  scale_color_brewer(palette = "Paired")+theme_bw()+
  scale_shape_manual(values=c(19,17,15,3,7),labels=c("Blood","Cerebellum","Striatum","Neocortex","S.Nigra"))+
  xlab("PC1 (11.16%)")+ylab("PC2 (10.27%)")+theme(legend.position = "none")

pc34<-ggplot(DV.p2[DV.p2$Sample!="DV2C",])+geom_point(aes(PC3,PC4,color=Individual, shape=Tissue),size=8)+
  geom_point(data=DV.p2[DV.p2$Sample=="DV2C",],aes(PC3,PC4,shape=Tissue,color=Individual),size=8)+
  scale_color_brewer(palette = "Paired")+theme_bw()+
  scale_shape_manual(values=c(19,17,15,3,7),labels=c("Blood","Cerebellum","Striatum","Neocortex","S.Nigra"))+
  xlab("PC3 (10.18%)")+ylab("PC4 (10.12%)")+theme(legend.position = "none")

#### 3. Plot smaller legend ####

#Make a similar plot with smaller points so that legend is smaller
plotSmallerLegend<-ggplot(DV.p2[DV.p2$Sample!="DV2C",])+geom_point(aes(PC3,PC4,color=Individual, shape=Tissue),size=5)+
  geom_point(data=DV.p2[DV.p2$Sample=="DV2C",],aes(PC3,PC4,shape=Tissue,color=Individual),size=5)+
  scale_color_brewer(palette = "Paired")+theme_bw()+
  scale_shape_manual(values=c(19,17,15,3,7),labels=c("Blood","Cerebellum","Striatum","Neocortex","S.Nigra"))+
  xlab("PC3 (10.18%)")+ylab("PC4 (10.12%)")

#Get the legend from that plot
extract.legend<-function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  return(legend)} 

legend<-extract.legend(plotSmallerLegend)

#### 4. Final plot ####

grid.arrange(pc12,pc34,legend,nrow=1,widths=c(0.85,0.85,0.15))

