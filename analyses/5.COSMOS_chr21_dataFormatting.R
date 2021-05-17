options(stringsAsFactors = F)
library(tidyr)

#Import results from all filters for chr21
rawFilterData<- read.delim("data/chr21_allCOSMOSFiltersResults.txt")

#Separate position-related filters and repeat so each position per sample is in a row
sepDataFirstCol<-separate(rawFilterData[,1:5],5,into=strsplit(colnames(rawFilterData)[5],"\\.")[[1]],sep=";")
sepDataFirstColLong<-do.call(rbind.data.frame,lapply(1:50, function(i) sepDataFirstCol))
#Separate SAMPLE cols and put each sample in one row 
sepDataLastColLong<-do.call(rbind.data.frame,lapply(7:56, function(col)
  cbind.data.frame(separate(as.data.frame(rawFilterData[,col]),1,
                            into=strsplit(colnames(rawFilterData)[6],"\\.")[[1]],sep=";"),
                   Sample=rep(colnames(rawFilterData)[col],nrow(rawFilterData)))))
#Join both DFs
sepData<-cbind.data.frame(sepDataFirstColLong,sepDataLastColLong)
rm(sepDataFirstCol,sepDataFirstColLong,sepDataLastColLong)

#Convert True and None to 0 (PASS filter) and False to 1 (FAIL filter)
sepData01<-sepData[,5:35]
sepData01[sepData01=="N"]<-0
sepData01[sepData01=="T"]<-0
sepData01[sepData01=="F"]<-1
sepData01<-as.data.frame(sapply(sepData01, as.numeric))
sepData01<-cbind.data.frame(sepData[,1:4],sepData01,Sample=sepData[,36])
colnames(sepData01)[5]<-"1000Gstrictmask"

#Only variant positions (most called in other samples or alternative homozygous)
variantPositions<-sepData01[sepData01$is01==0,]

write.table(variantPositions,file="data/chr21_VariantPositions_FilterInfo.txt",row.names=F,quote=F,sep="\t")


