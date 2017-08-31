source("../dataAccess/RNASeqData.R")

samp.mappings<-rnaSeqFiles[c("sampleIdentifier","rnaSeq","nf1Genotype")]

tpm=rnaKallistoMatrix()
counts=rnaKallistoMatrix(metric = 'est_counts')

##now begin to do basic plotting
library(pheatmap)
library(ggbiplot)
plotMostVariableTranscripts<-function(count.mat,num=100,metric='tpm'){
  samp.names=samp.mappings$sampleIdentifier[match(colnames(count.mat),samp.mappings$rnaSeq)]
  samp.gen=samp.mappings$nf1Genotype[match(colnames(count.mat),samp.mappings$rnaSeq)]
  names(samp.gen)<-samp.names
  colnames(count.mat)<-samp.names
  var.vals<-apply(count.mat,1,var,na.rm=T)
  print(num)
  var.mat<-count.mat[order(var.vals,decreasing=T)[1:num],]
  pheatmap(var.mat,cellwidth=10,cellheight=10,annotation_col=data.frame(Genotype=samp.gen),filename=paste('top',num,'mostVariableGenesBy',metric,'.png',sep=''))
}

getGTAssociatedTranscripts<-function(count.mat,num=100,metric='tpm'){
  samp.names=samp.mappings$sampleIdentifier[match(colnames(count.mat),samp.mappings$rnaSeq)]
  samp.gen=samp.mappings$nf1Genotype[match(colnames(count.mat),samp.mappings$rnaSeq)]
  names(samp.gen)<-samp.names
  colnames(count.mat)<-samp.names
  ##make it abinary test
  samp.gen[which(samp.gen=='+/-')]='+/+'
  ##now for each gene, evaluate significance of linear model predicting genotype
  lmvals<-apply(count.mat,1,function(x)
    summary(lm(Expr~Genotype,data=data.frame(Genotype=as.factor(samp.gen),Expr=x))))
  pvals<-sapply(lmvals,function(x){
    fs=x$fstatistic
    pf(fs[1],fs[2],fs[3],lower.tail=F)[1]
  })
  fdr.vals<-p.adjust(pvals,method='BH')
  #sig<-which(fdr.vals<0.15)
  #names(sig)<-sapply(names(sig),function(x) paste(unlist(strsplit(x,split='.',fixed=T))[1:2],collapse='.'))
  ##now plo those associated
  
  var.mat<-count.mat[order(pvals,decreasing=F)[1:num],]
  pheatmap(var.mat,cellwidth=10,cellheight=10,annotation_col=data.frame(Genotype=samp.gen),filename=paste('top',num,'mostGTCorrelatedGenesBy',metric,'.png',sep=''))
  
}

plotNF1Transcripts<-function(count.mat,metric='tpm'){
  samp.names=samp.mappings$sampleIdentifier[match(colnames(count.mat),samp.mappings$rnaSeq)]
  samp.gen=samp.mappings$nf1Genotype[match(colnames(count.mat),samp.mappings$rnaSeq)]
  names(samp.gen)<-samp.names
  orig.samp.gen<-samp.gen
  colnames(count.mat)<-samp.names
  nf1=grep("^NF1.ENST0",rownames(count.mat))
  print(paste("Found",length(nf1),'NF1 transcripts'))
  n.mat<-count.mat[nf1,]
  pheatmap(n.mat,cellwidth=10,cellheight=10,annotation_col=data.frame(Genotype=orig.samp.gen),filename=paste('NF1TranscriptsBy',metric,'.png',sep=''))
  
}


plotPCA<-function(count.mat,metric='tpm'){
  samp.names=samp.mappings$sampleIdentifier[match(colnames(count.mat),samp.mappings$rnaSeq)]
  samp.gen=samp.mappings$nf1Genotype[match(colnames(count.mat),samp.mappings$rnaSeq)]
  names(samp.gen)<-samp.names
  colnames(count.mat)<-samp.names
  
  zv<-which(apply(count.mat,1,var)==0)
  if(length(zv)>0)
    count.mat=count.mat[-zv,]
  pn<-prcomp(t(count.mat),center=T,scale=T)
  png(paste(metric,'valuesinRNASeqData.png',sep=''))
  p<-ggbiplot(pn,groups=samp.gen,var.axes=F)
  print(p)
  dev.off()
  
}

##plot PCA
plotPCA(log2(tpm+0.01),'tpm')
plotPCA(log2(counts+1),'est_counts')

##plotHeatmaps
getSum <- function(df,refCol){  
  result <- do.call(rbind,lapply(unique(refCol),FUN =function(x){
    temp <- df[df$genes == x,]
    temp$genes <- NULL
    final <- matrix(colSums(temp), nrow=1)
    row.names(final) <- x
    return(final)
  }))
  result <- as.data.frame(result)
  colnames(result) <- colnames(df)[-length(colnames(df))]
  return(result)
}

tpm.genes <- tpm
tpm.genes$genes <- sub("\\.ENST.+","",row.names(tpm.genes))
tpm.genes <- tpm.genes[tpm.genes$genes != "" & tpm.genes$genes != "NA",]
tpm.genes <- getSum(tpm.genes, tpm.genes$genes)
tpm.genes <- log2(tpm.genes+0.01)

counts.genes <- counts
counts.genes$genes <- sub("\\.ENST.+","",row.names(counts.genes))
counts.genes <- counts.genes[counts.genes$genes != "" & counts.genes$genes != "NA",]
counts.genes <- getSum(counts.genes, counts.genes$genes)
counts.genes <- log2(counts.genes+0.01)
# by genes
plotMostVariableTranscripts(tpm.genes,100,'tpm')
plotMostVariableTranscripts(counts.genes,100,'est_counts')

# by transcripts
plotMostVariableTranscripts(log2(tpm+1),100,'tpm')
plotMostVariableTranscripts(log2(counts+1),100,'est_counts')

getGTAssociatedTranscripts(tmp.new,100,'tpm')
getGTAssociatedTranscripts(log2(counts+0.01),100,'est_counts')

plotNF1Transcripts(log2(tpm+0.01),'tpm')
plotNF1Transcripts(log2(counts+0.01),'est_counts')

#### upload files
uploaded.folder<- 'syn5731446'
scripturl <- 'https://raw.githubusercontent.com/Sage-Bionetworks/NTAP/master/pnfCellLines/analysis/prelimRnaAnalysis.R'

uploadFile2Synapse <- function(file,parentId,usedScript,forceVersion=T){
  synStore(File(path = file,parentId = parentId),
           used=list(list(url=usedScript,wasExecuted=TRUE)), forceVersion=forceVersion)
}


