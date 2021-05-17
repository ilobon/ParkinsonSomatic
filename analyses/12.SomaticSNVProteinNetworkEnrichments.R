options(stringsAsFactors = F)

library(ggplot2)
library(biomaRt)
library(WebGestaltR)

`%notin%`<-Negate(`%in%`)

#Data
ValidationResults<-read.delim2("data/CandidatesInfoValidations.txt")

infoValidationsMultiDel<-ValidationResults[ValidationResults$Decision%in%c("MultiTissue") & 
                                           ValidationResults$Type%notin%c("Intron","Synonymous"),]

#### 1. Enrichment of deleterious brain gene set ####

#Get first in gene/gene or gene,gene
genenames<-infoValidationsMultiDel$Gene
genenames<-sapply(strsplit(genenames,"[[:punct:]]"), function(x) x[[1]][1])

#Enrichment original set
enrichmentBrainDel<-WebGestaltR(enrichmentMethod="ORA",
                                 organism="hsapiens",
                                 enrichDatabase=c("geneontology_Biological_Process",
                                                  "geneontology_Cellular_Component",
                                                  "geneontology_Molecular_Function"),
                                 referenceSet="genome_protein-coding",
                                 interestGene=genenames,
                                 interestGeneType="genesymbol",
                             isOutput = F,
                             sigMethod="top",topThr=10)

enrichmentBrainDel<-enrichmentBrainDel[order(enrichmentBrainDel$enrichmentRatio),]
enrichmentBrainDel$description<-factor(enrichmentBrainDel$description,levels=unique(enrichmentBrainDel$description),ordered=T)
enrichmentBrainDel$database<-factor(enrichmentBrainDel$database,levels=unique(enrichmentBrainDel$database),ordered=T)
levels(enrichmentBrainDel$database)<-c("BioP", "CellC", "MolF")

#My gene set enrichment - supplementary
ggplot(enrichmentBrainDel)+geom_bar(aes(description,enrichmentRatio,fill=FDR<0.05),stat="identity")+
  scale_fill_manual(values="#B0C5DE")+
  xlab("")+
  ylab("Enrichment ratio")+
  facet_grid(database~.,scales="free_y",space = "free")+
  theme_bw()+
  coord_flip()

#### 2. Enrichment of expanded gene set ####

#Convert to entrezID with Biomart
ensembl<-useMart(biomart="ensembl",dataset="hsapiens_gene_ensembl")
geneName_to_Entrez<-getBM(mart=ensembl,
      attributes=c("external_gene_name","entrezgene_id"),
      filters="external_gene_name",
      values=genenames)
entrez<-geneName_to_Entrez$entrezgene_id

#Prepare STRING data

#Import alias info from STRING #https://stringdb-static.org/download/protein.aliases.v11.0.txt.gz
aliasesString<-read.delim("9606.protein.aliases.v11.0.txt.gz", header=FALSE, comment.char="#")
#Import human data table #https://stringdb-static.org/download/protein.links.detailed.v11.0.txt.gz
humanString<-read.csv("9606.protein.links.detailed.v11.0.txt.gz", sep="")

#Translate entrez to STRING IDs
stringGenes<-aliasesString$V1[aliasesString$V2%in%entrez]

#Get gene name from the STRING id
getGeneName<-function(stringID){
  aliasesString[aliasesString$V1==stringID,2][
    grep("Ensembl_WikiGene", aliasesString$V3[which(aliasesString$V1==stringID)])]}

genenames<-unlist(lapply(stringGenes,getGeneName))

#STRING data for each gene
genesStringData<-lapply(1:length(stringGenes), function(i) humanString[humanString$protein1==stringGenes[i],])
genesStringDataDF<-do.call(rbind.data.frame,genesStringData)

## Extend network based on sum co-expression score

#Get scores for every protein co-expressed with geneset
protein2Coexpression<-lapply(unique(genesStringDataDF$protein2), function(x) 
  genesStringDataDF$coexpression[genesStringDataDF$protein2==x])

#Filter by sum score > 900
protein2CoexpPass<-unlist(lapply(unique(genesStringDataDF$protein2)
                                 [sapply(protein2Coexpression,sum)>900],getGeneName))

##Enrichment
enrichmentTop25<-WebGestaltR(enrichmentMethod="ORA",
                                 organism="hsapiens",
                                 enrichDatabase=c("geneontology_Biological_Process",
                                                  "geneontology_Cellular_Component",
                                                  "geneontology_Molecular_Function",
                                                  "phenotype_Human_Phenotype_Ontology",
                                                  "disease_GLAD4U"),
                                 referenceSet="genome_protein-coding",
                                 interestGene=c(genenames,protein2CoexpPass),
                                 interestGeneType="genesymbol",
                                 sigMethod="top",
                                 topThr=25,
                                 isOutput = F)

enrichmentTop25<-enrichmentTop25[order(enrichmentTop25$enrichmentRatio),]
enrichmentTop25$description<-factor(enrichmentTop25$description,levels=unique(enrichmentTop25$description),ordered=T)
enrichmentTop25$database<-factor(enrichmentTop25$database,levels=unique(enrichmentTop25$database),ordered=T)

levels(enrichmentTop25$database)<-c("Bio Process", "Cellular Component")

extendedEnrichmentPlot<-ggplot(enrichmentTop25)+geom_bar(aes(description,enrichmentRatio,fill=FDR<0.05),stat="identity")+
  scale_fill_manual(values="#4682B4")+
  xlab("")+
  ylab("Enrichment ratio")+
  facet_grid(database~.,scales="free_y",space = "free")+
  theme_bw()+
  coord_flip()+
  theme(axis.text=element_text(size=12),
        axis.title.x=element_text(size=14),
        strip.text=element_text(size=12),
        panel.grid.major.x=element_blank(),
        panel.grid.minor.x=element_blank(),
        panel.grid.major.y=element_blank(),
        panel.grid.minor.y=element_blank())



## All enriched terms

enrichmentExtended<-WebGestaltR(enrichmentMethod="ORA",
                                 organism="hsapiens",
                                 enrichDatabase=c("geneontology_Biological_Process",
                                                  "geneontology_Cellular_Component",
                                                  "geneontology_Molecular_Function",
                                                  "phenotype_Human_Phenotype_Ontology",
                                                  "disease_GLAD4U"),
                                 referenceSet="genome_protein-coding",
                                 interestGene=c(genenames,protein2CoexpPass),
                                 interestGeneType="genesymbol",
                                 isOutput = F)

databaseOrder<-unique(enrichmentExtended$database)[c(4,2,3,1,5)]
enrichmentTableOrdered<-enrichmentExtended[order(match(enrichmentExtended$database,databaseOrder),enrichmentExtended$FDR),]
enrichmentTableOrdered<-enrichmentTableOrdered[,c(11,2,4,5,6,7,8,9)]
enrichmentTableOrdered$expect<-signif(enrichmentTableOrdered$expect,3)
enrichmentTableOrdered$enrichmentRatio<-signif(enrichmentTableOrdered$enrichmentRatio,3)
enrichmentTableOrdered$pValue<-signif(enrichmentTableOrdered$pValue,3)
enrichmentTableOrdered$FDR<-signif(enrichmentTableOrdered$FDR,3)
write.table(enrichmentTableOrdered,"data/allEnrichmentsExtended.txt",row.names=F,quote=F,sep="\t")




