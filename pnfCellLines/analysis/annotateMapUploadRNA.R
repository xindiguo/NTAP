###get all RNA-Seq files from local directory where kallisto was run
##map ENS transcripts to gene id
##create new tab-delimited file with transcript AND gene id
##annotate with appropriate cell line and upload

require(synapseClient)
library(data.table)
synapseLogin()

samp.mapping=read.table(synGet('syn10234072')@filePath,sep=',',header=T, stringsAsFactors = FALSE)

##################### GRCH38 ######################
local.dir <- "./rnaSeq/grch38Aligned/kallisto_result/"
filedirs <- list.files(local.dir)
all.dat<-lapply(filedirs,function(x) as.data.frame(fread(paste(local.dir,x,'abundance.tsv',sep='/'),sep='\t')))
names(all.dat)<-filedirs

##first get all gene names, add to files
all.genes=all.dat[[1]][,1]

require(biomaRt)
ensembl=useMart("ensembl",dataset="hsapiens_gene_ensembl")
filters = listFilters(ensembl)
attributes = listAttributes(ensembl)

epep="ensembl_transcript_id"
egene='hgnc_symbol'
gene.mapping<-getBM(attributes=c(epep,egene),filters=c(epep),values=as.list(all.genes),mart=ensembl)

#now create/write new files
newfiles<-lapply(all.dat,function(x){
  HGNCSymbol=gene.mapping[match(x[,1],gene.mapping[,1]),2]
  newf=data.frame(HGNCSymbol,x)
  return(newf)
})
names(newfiles)=names(all.dat)

##write files and store on Synapse
file.res<-lapply(filedirs,function(x){
  filename=paste(x,'_RNASeq_Kallisto_grch38_quants.tsv',sep='')
  write.table(newfiles[[x]],file=filename,sep='\t',row.names=F)
  newf=File(path=filename,parentId='syn5579783')
  synSetAnnotations(newf) <- list(cellType = "cultured cell",
                                  species = "Human",
                                  study = "Cell Culture",
                                  dataType = "geneExpression",
                                  dataSubtype = "processed",
                                  fileFormat = "tsv",
                                  fundingAgency = "NTAP",
                                  isCellLine = "true",
                                  assay = "rnaSeq",
                                  sampleIdentifier = samp.mapping$sampleIdentifier[samp.mapping$label==x])
  newf=synStore(newf,executed='https://raw.githubusercontent.com/Sage-Bionetworks/NTAP/master/pnfCellLines/analysis/annotateMapUploadRNA.R')
})
##################### GENCODE #####################
#### tsv ####
local.dir <- "./rnaSeq/gencodeV24Aligned/kallisto_result/"
filedirs <- list.files(local.dir)
all.dat<-lapply(filedirs,function(x) as.data.frame(fread(paste(local.dir,x,'abundance.tsv',sep='/'),sep='\t')))
names(all.dat)<-filedirs

##first get all gene names, add to files
all.genes=all.dat[[1]][,1]
newfiles<-lapply(all.dat,function(x){
  dat<-x[,-1]
  ids<-t(sapply(x[,1],function(x) unlist(strsplit(x,split='|',fixed=T))))
  colnames(ids)<-c('EnsGene','EnsTrans','OttGene','OttTrans','HugoTrans','HugoSymbol','Length','Descrip')
  new.dat<-data.frame(ids,dat)
  return(new.dat)
})

##write files and store on Synapse
file.res<-lapply(filedirs,function(x){
  filename=paste(x,'_RNASeq_Kallisto_gencodev24_quants.tsv',sep='')
  write.table(newfiles[[x]],file=filename,sep='\t',row.names=F)
  newf=File(path=filename,parentId='syn5579785')
  synSetAnnotations(newf) <- list(cellType = "cultured cell",
                                  species = "Human",
                                  study = "Cell Culture",
                                  dataType = "geneExpression",
                                  dataSubtype = "processed",
                                  fileFormat = "tsv",
                                  fundingAgency = "NTAP",
                                  isCellLine = "true",
                                  assay = "rnaSeq",
                                  sampleIdentifier = samp.mapping$sampleIdentifier[samp.mapping$label==x])
  newf=synStore(newf,executed='https://raw.githubusercontent.com/Sage-Bionetworks/NTAP/master/pnfCellLines/analysis/annotateMapUploadRNA.R')
})

#### h5 ####
local.dir <- "./rnaSeq/gencodeV24Aligned/kallisto_result/"
filedirs <- list.files(local.dir)

##write files and store on Synapse
file.res<-lapply(filedirs,function(x){
  filename=paste(x,'_RNASeq_Kallisto_gencodev24_quants.h5',sep='')
  file.copy(paste(local.dir,x,'abundance.h5',sep='/'),filename)
  newf=File(path=filename,parentId='syn5579785')
  synSetAnnotations(newf) <- list(cellType = "cultured cell",
                                  species = "Human",
                                  study = "Cell Culture",
                                  dataType = "geneExpression",
                                  dataSubtype = "processed",
                                  fileFormat = "h5",
                                  fundingAgency = "NTAP",
                                  isCellLine = "true",
                                  assay = "rnaSeq",
                                  sampleIdentifier = samp.mapping$sampleIdentifier[samp.mapping$label==x])
  newf=synStore(newf,executed='https://raw.githubusercontent.com/Sage-Bionetworks/NTAP/master/pnfCellLines/analysis/annotateMapUploadRNA.R')
})
