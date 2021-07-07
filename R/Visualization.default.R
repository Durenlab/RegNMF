
#library(GenomicRanges)
#library(Gviz)
#library(GenomicInteractions)

Visualization.default <- function(wholef,peakf,regf,chr,from,to,clusterlist,width=6,height=4){


filename=paste0(wholef,"features.tsv")
data=read.table(filename,sep='\t')
x=match(data$V3,"Peaks")
data=data[!is.na(x),]
x=match(data$V4,chr)
data=data[!is.na(x),]
gr=GenomicRanges::GRanges(seqnames = data$V4,ranges = IRanges(start = data$V5,end = data$V6))
gtrack <- GenomeAxisTrack()

atrackw<-AnnotationTrack(gr,name="whole")

Gene <- BiomartGeneRegionTrack(genome = "hg38",chromosome = chr,start=from,end=to,name="Gene",collapseTranscripts = "meta",transcriptAnnotation = "symbol")


cluster=list();
atrack=list();
for (i in clusterlist) {
  file=paste0(peakf,i,"_peaks.narrowPeak")
  data=read.table(file)
  x=match(data$V1,chr)
  data=data[!is.na(x),]
  cluster[[i]]=GRanges(seqnames = data$V1,ranges = IRanges(start = data$V2,end = data$V3))
  atrack[[i]]= AnnotationTrack(cluster[[i]])
}



ftmp=factor(Gene@range@elementMetadata@listData$symbol)
flevel=levels(ftmp)
RL=rep(0,length(levels(ftmp)))
ftmp=match(ftmp,levels(ftmp))
for (i in 1:length(ftmp)) {
  start=Gene@range@ranges@start[i]
  if(RL[ftmp[i]]==0){
    RL[ftmp[i]]=start
  }
  RL[ftmp[i]]=min(start,RL[ftmp[i]])
}
#ftmp=RL[ftmp]

REcluster=list();
RETG=list();
interaction_track=list();

for (i in clusterlist) {
  file=paste0(regf,i,".bed")
  data=read.table(file)
  x=match(data$V1,chr)
  data=data[!is.na(x),]
  x=match(data$V4,flevel)
  data=data[!is.na(x),]
  if(nrow(data)<=0){
    next;
  }
  print(i)
  x=as.integer(na.omit(x))
  data$V5=RL[x]
  REcluster[[i]]=GRanges(seqnames = data$V1,ranges = IRanges(start = data$V2,end = data$V3))
  RETG[[i]]=GRanges(seqnames = rep(chr,length(data$V1)) ,IRanges(start = data$V5,end=data$V5+100))
  a=Pairs(REcluster[[i]],RETG[[i]])
  b=makeGInteractionsFromGRangesPairs(a)
  c=GenomicInteractions(b)
  name=paste0("R",i)
  interaction_track[[i]] = InteractionTrack(c, name = name, chromosome = chr)
}



####Make Plot
Plist=list()
Plist=append(Plist,gtrack)
for (i in clusterlist) {
  if(class(interaction_track[[i]])=="NULL") {next};
  Plist=append(Plist,interaction_track[[i]])
  Plist=append(Plist,atrack[[i]])
}
Plist=append(Plist,atrackw)
Plist=append(Plist,Gene)

pdf("./result.pdf",width = width,height = height)
plotTracks(Plist,from=from,to=to,plot.outside=FALSE)
dev.off()


}
