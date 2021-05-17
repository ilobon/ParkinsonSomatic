options(stringsAsFactors=F)
library(WebGestaltR)
library(ggplot2)

#### 1. Import and format data ####

delGenes<-read.table("data/perInd.deleteriousGenes.txt", quote="\"", comment.char="")
#V1 is colnumber in file, make factor and translate to individual ID
delGenes$V1<-factor(delGenes$V1,levels=unique(delGenes$V1),ordered=T)
levels(delGenes$V1)<-paste("DV",c(10,1,3:9),sep="")

#### 2. WebGestalt enrichment ####

enrichment<-WebGestaltR(enrichmentMethod="ORA",
                        organism="hsapiens",
                        enrichDatabase=c("geneontology_Cellular_Component",
                                         "geneontology_Molecular_Function"),
                        referenceSet="genome_protein-coding",
                        interestGene=delGenes$V2,
                        interestGeneType="genesymbol",
                        isOutput = F,
                        minNum=2)

#### 3. Retrieve number of genes per individual and enriched term ####

#For each enriched term, return the individual IDs with variants overlapping its genes
indsPerTerm=lapply(1:nrow(enrichment), function(i) 
  as.character(delGenes$V1[which(!is.na(
    match(delGenes$V2,strsplit(enrichment$userId[i],";")[[1]])))]))

#Create a dataframe with row for each term and individual, enrichment ratio and number of genes overlapping
termsInd<-do.call(rbind.data.frame,lapply(1:nrow(enrichment), function(i) 
  cbind.data.frame(description=rep(enrichment$description[i],length(unique(indsPerTerm[[i]]))),
                   er=rep(enrichment$enrichmentRatio[i],length(unique(indsPerTerm[[i]]))),
                   inds=names(table(indsPerTerm[[i]])),
                   nGenes=as.numeric(unname(table(indsPerTerm[[i]]))))))

#Make enrichment description a factor
termsInd$description<-factor(termsInd$description,levels=unique(termsInd$description[order(termsInd$er)]),ordered=T)

#Order inds by age at death
age_onset<-read.delim("data/age_onset.txt")
termsInd$inds<-factor(termsInd$inds,levels=age_onset$Patient.ID[order(age_onset$Age)],ordered=T)

#### 4. Final plot ####

ggplot(termsInd)+geom_point(aes(inds,description,size=nGenes,color=factor(nGenes)))+
  scale_color_manual(values=c("#33a02c","#1f78b4","#6a3d9a","#e31a1c"),name="N of genes")+
  geom_text(aes(10.5,description,label=signif(er,3),size=1.5))+
  geom_text(aes(10.5,8,label="Enrichment",size=1.5))+
  geom_text(aes(10.5,7.7,label=" ratio",size=1.5))+
  coord_cartesian(xlim=c(1,9),ylim=c(1,7),clip = 'off')+
  theme_bw()+
  ylab("")+
  xlab("")+
  theme(legend.position="bottom",plot.margin = unit(c(1,4,1,1), "lines"))+
  guides(size=FALSE)




