options(stringsAsFactors = F)
library(tidyr)
library(UpSetR)
library(RColorBrewer)
library(pheatmap)

#### 1. Import and format data ####

#Import data
variantPositions<- read.delim("data/chr21_VariantPositions_FilterInfo.txt",
                              check.names=FALSE)

#On target / off target
onTarget<-variantPositions[variantPositions$onTarget==0,]
offTarget<-variantPositions[variantPositions$onTarget==1,]

#### UpsetR germline heterozygous ####

#Dataframe with relevant filters
ghInfo<-onTarget[,c("PDpanel","PON","Binomial","VAF","BinomialInd")]

#Format for plot
ghInfo$PDpanel<-as.character(ghInfo$PDpanel)
ghInfo$PDpanel[ghInfo$PDpanel==1]<-"NotPass"
ghInfo$PDpanel[ghInfo$PDpanel==0]<-"Pass"

ghInfo$PON<-as.character(ghInfo$PON)
ghInfo$PON[ghInfo$PON==1]<-"NotPass"
ghInfo$PON[ghInfo$PON==0]<-"Pass"

#UpsetR plot
ghPlot<-upset(ghInfo, 
              mainbar.y.label="", 
              sets.x.label="Filtered positions",
              query.legend="bottom",
              queries = list(
                list(query=elements,params=list("PDpanel","NotPass"),color="#e7d4e8",active=T,query.name="PD panel"),
                list(query=elements,params=list("PON","NotPass"),color="#762a83",query.name="PON")),
              text.scale=c(1, 1.2, 1.4, 0.9, 1.7, 1.4))


#### UpsetR CNVs ####

#Dataframe with relevant filters
cnvInfo<-onTarget[onTarget$Binomial==0 & onTarget$BinomialInd==0 & onTarget$VAF==0 & onTarget$VAFhc==0,
                  c("PDpanel","PON","1000Gstrictmask","DepthRange","SegmentalDups")]

#Format for plot
cnvInfo$PDpanel<-as.character(cnvInfo$PDpanel)
cnvInfo$PDpanel[cnvInfo$PDpanel==1]<-"NotPass"
cnvInfo$PDpanel[cnvInfo$PDpanel==0]<-"Pass"

cnvInfo$PON<-as.character(cnvInfo$PON)
cnvInfo$PON[cnvInfo$PON==1]<-"NotPass"
cnvInfo$PON[cnvInfo$PON==0]<-"Pass"

cnvInfo$`1000Gstrictmask`<-as.character(cnvInfo$`1000Gstrictmask`)
cnvInfo$`1000Gstrictmask`[cnvInfo$`1000Gstrictmask`==1]<-"NotPass"
cnvInfo$`1000Gstrictmask`[cnvInfo$`1000Gstrictmask`==0]<-"Pass"

#UpsetR plot
cnvPlot<-upset(cnvInfo,
      mainbar.y.label="",
      sets.x.label="Filtered positions",
      query.legend = "bottom",
      queries=list(
        list(query=elements,params=list("PDpanel","NotPass"),color="#d9f0d3",query.name="PD panel"),
        list(query=elements,params=list("1000Gstrictmask","NotPass"),color="#a6dba0",active=T,query.name="1000G strict mask"),
        list(query=elements,params=list("PON","NotPass"),color="#1b7837",query.name="PON")),
      text.scale=c(1, 1.2, 1.4, 0.9, 1.7, 1.4))



#### Supp figure CNV and other features #####

#Fail DepthRange SegDups
cnvInfoExt<-onTarget[onTarget$Binomial==0 & onTarget$BinomialInd==0 & onTarget$VAF==0 & onTarget$VAFhc==0 &
                       onTarget$DepthRange==1 & onTarget$SegmentalDups==1,
                     c("PDpanel","PON","1000Gstrictmask","DepthRange","SegmentalDups",
                       "Mappability","nVarsReadLength","Mismatches","Clipping")]

cnvInfoExt$PDpanel<-as.character(cnvInfoExt$PDpanel)
cnvInfoExt$PDpanel[cnvInfoExt$PDpanel==1]<-"NotPass"
cnvInfoExt$PDpanel[cnvInfoExt$PDpanel==0]<-"Pass"

cnvInfoExt$PON<-as.character(cnvInfoExt$PON)
cnvInfoExt$PON[cnvInfoExt$PON==1]<-"NotPass"
cnvInfoExt$PON[cnvInfoExt$PON==0]<-"Pass"

cnvInfoExt$`1000Gstrictmask`<-as.character(cnvInfoExt$`1000Gstrictmask`)
cnvInfoExt$`1000Gstrictmask`[cnvInfoExt$`1000Gstrictmask`==1]<-"NotPass"
cnvInfoExt$`1000Gstrictmask`[cnvInfoExt$`1000Gstrictmask`==0]<-"Pass"

upset(cnvInfoExt,
      mainbar.y.label="CNV filters and related features",
      sets.x.label="Positions filtered",
      query.legend = "bottom",
      queries=list(
        list(query=elements,params=list("PDpanel","NotPass"),color="#d9f0d3",query.name="PD panel"),
        list(query=elements,params=list("1000Gstrictmask","NotPass"),color="#a6dba0",active=T,query.name="1000G strict mask"),
        list(query=elements,params=list("PON","NotPass"),color="#1b7837",query.name="PON")),
      text.scale=c(1.2, 1, 1.4, 0.8, 1.7, 1.2))


#### Heatmap InOut ####

#Get ratio with pseudocount 
getRatioInOut<-function(df, mainFilterCol, secondaryFilterCol){
  failBoth_outOf_failMain<-((sum(df[,mainFilterCol]==1 & df[,secondaryFilterCol]==1)/
                               sum(df[,mainFilterCol]==1))+0.001)
  failSec_outOf_passMain<-((sum(onTarget[,mainFilterCol]==0 & onTarget[,secondaryFilterCol]==1)/
                              sum(onTarget[,mainFilterCol]==0))+0.001)
  return(failBoth_outOf_failMain/failSec_outOf_passMain)
}

inOut<-do.call(rbind.data.frame, lapply(5:35, function(i) 
  unlist(lapply(5:35, function(j) getRatioInOut(onTarget,i,j)))))

colnames(inOut)<-colnames(onTarget)[5:35]
rownames(inOut)<-colnames(onTarget)[5:35]

#Order by number of failed positions
sizeOrder<-order(apply(onTarget[,5:35],2,sum),decreasing = T)
inOut<-inOut[sizeOrder,sizeOrder]

#Convert to log2
invariableFilts<-which(sapply(1:ncol(inOut), function(i) all(inOut[,i]==1,na.rm=T)))
inOutVariable<-inOut[-invariableFilts,-invariableFilts]
inOutVariableLog2<-log2(inOutVariable)

#Remove uninformative filters
uninformative<-which(colnames(inOutVariableLog2)%in%c("StrandRatio","VAFhc","MinAltDepthSingleSample"))
inOutVariableLog2<-inOutVariableLog2[-uninformative,-uninformative]

#Annotation with filter size
annotation<-data.frame(FilteredVars=unname(apply(onTarget[,5:35],2,sum)))
rownames(annotation)<-colnames(onTarget[5:35])

#Color palette with 0 in the middle + 40 on each side
colors=c(colorRampPalette(rev(brewer.pal(n=7,name="RdBu")))(100)[c(1:40)],
         "white",
         colorRampPalette(rev(brewer.pal(n=7,name="RdBu")))(100)[c(61:100)])

#Make heatmap
heatmapFiltersPlot<-pheatmap(as.matrix(inOutVariableLog2),color=colors,breaks=seq(-10,9.75,0.25),annotation_row=annotation,
         cluster_rows = T,
         cluster_cols = F,
         border_color=FALSE,
         legend_breaks=c(-8,0,8),legend_labels = c("-8","0","8"))

#Change order of cols to be that of the clustered rows
inOutVariableLog2<-inOutVariableLog2[,heatmapFiltersPlot$tree_row$order]

#Repeat plot
heatmapFiltersPlot<-pheatmap(as.matrix(inOutVariableLog2),
                             color=colors,
                             breaks=seq(-10,9.75,0.25),
                             annotation_row=annotation,
                             cluster_rows = T,
                             cluster_cols = F,
                             border_color=FALSE,
                             legend_breaks=c(-8,0,8),
                             legend_labels = c("-8","0","8"),
                             fontsize=14)


