options(stringsAsFactors=F)

library(pheatmap)
library(RColorBrewer)
library(plyr)
library(ggrepel)

`%notin%`<-Negate(`%in%`)

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

#### Input data ####

ValidationResults<-read.delim2("data/CandidatesInfoValidations.txt")
age_onset<-read.delim("data/age_onset.txt")
infoValidations<-merge(ValidationResults,age_onset,by.x=10,by.y=1,all.x=T)

#Add collapsed substitutions
nucleotidePairs<-cbind.data.frame(a=c("A","C","G","T"),b=c("T","G","C","A"))
Substitutions<-sapply(1:nrow(infoValidations), function(i) 
  collapseSubs(infoValidations$Ref[i],infoValidations$Alt[i]))
infoValidations<-cbind.data.frame(infoValidations,Substitution=Substitutions)


#### VAFs heatmap ####

#Filter multi tissue variants
infoValidationsSubset<-infoValidations[which(infoValidations$Decision%in%c("MultiTissue")),]

#Convert VAFs to numeric
tissues<-c("B", "C", "E", "N", "S")
ampliconVAFcols<-which(colnames(infoValidationsSubset)%in%tissues)
for (i in ampliconVAFcols){infoValidationsSubset[,i]<-as.numeric(infoValidationsSubset[,i])}
#VAF matrix
VAFs<-t(as.matrix(infoValidationsSubset[c(ampliconVAFcols)]))

#Validated tissues per variant
valTissues<-lapply(1:nrow(infoValidationsSubset), function(i)
  unname(sapply(strsplit(infoValidationsSubset$ValidatedSamples[i],",")[[1]], 
                function(x) substr(x,nchar(x),nchar(x)))))

#Change amplicon VAF to 0 if not validated in that tissue
for (i in 1:ncol(VAFs)){
  VAFs[which(tissues%notin%valTissues[[i]]),i]<-0
}


#Order data by age then total VAF
goodOrder<-order(infoValidationsSubset$Age,sapply(1:ncol(VAFs), function(i) sum(VAFs[,i]==0)))
#colnames trace the row number 
VAFs<-VAFs[,goodOrder]

#Age annotation
annotation_age<-cbind.data.frame(AgeAtDeath=infoValidations$Age)

#purples
mycolors<-c('#efedf5','#dadaeb','#bcbddc','#9e9ac8','#807dba','#6a51a3','#54278f','#3f007d')
names(mycolors)<-sort(unique(annotation_age$AgeAtDeath))
mycolors<-list(AgeAtDeath=mycolors)

rownames(VAFs)<-c("Blood","Cerebellum","Striatum","Neocortex","S.nigra")

#Genes in colnames and black in 0 VAF

#colors: 64 + black for 0
colors<-c("#000000",colorRampPalette(rev(brewer.pal(n=11,name="RdYlBu")))(65))


VAFheatmap<-pheatmap(VAFs,border_color=FALSE,cluster_cols=F,
         color=colors,breaks=c(0,seq(min(VAFs[VAFs!=0])*0.9,0.32,0.005)),
         labels_col=infoValidationsSubset$Gene[goodOrder],
         annotation_col=annotation_age,
         annotation_colors=mycolors)






#### Age correlations ####

#Subset putative deleterious variants
infoValidationsDel<-infoValidationsSubset[which(infoValidationsSubset$Type%notin%c("Intron","Synonymous")),]

#Variant count per patient
varCounts<-count(infoValidationsDel$Individual)
varCountsAge<-merge(age_onset,varCounts,by.x = 1, by.y = 1)

#Linear model
X=varCountsAge$Age
Y=varCountsAge$freq
model<-lm(Y ~ X)
lm_total_age<-predict(model, list(X=X))
lm_total_age<-cbind.data.frame(varCountsAge,lmTotal=lm_total_age)


ageCorrelationPlot<-ggplot(lm_total_age)+
  geom_point(aes(Age,freq,shape=Sex),color="#bdbdbd",size=6)+
  geom_line(aes(Age,lmTotal),group=1,color="black",size=1)+
  geom_label_repel(aes(Age,freq,label=Patient.ID),
                   fontface = 'bold', color = 'black',
                   point.padding =0.4,
                   segment.alpha=0)+
  scale_x_continuous(breaks=seq(76,90,2))+
  scale_y_continuous(breaks=seq(min(lm_total_age$freq),max(lm_total_age$freq),1))+
  ylab("Deleterious brain somatic SNVs")+
  xlab("Age at death")+
  annotate("text",x=85,y=max(varCountsAge$freq),size=4.5,
           label=paste0("Corr=",round(cor(varCountsAge$freq,varCountsAge$Age),2)))+
  theme_bw()+
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=12),
        legend.title=element_text(size=12), 
        legend.text=element_text(size=11),
        panel.grid.major.x=element_blank(),
        panel.grid.minor.x=element_blank(),
        panel.grid.major.y=element_blank(),
        panel.grid.minor.y=element_blank())



