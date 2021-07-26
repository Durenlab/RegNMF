callpeak.default<-function(outfolder,fragment,barcord_cluster_whole,oldRegFolder,macs2path,awkpath,cluster,clusterL){

barcord_clusterF=paste0(outfolder,"barcord_cluster/")

#cmd=paste0("mkdir -p ",barcord_clusterF,";rm -f ",barcord_clusterF,"*")
#system(command=cmd)

if(!dir.exists(barcord_clusterF)){
    dir.create(barcord_clusterF,recursive = TRUE )
} else if (length(dir(barcord_clusterF))) {
   file.remove(paste0(barcord_clusterF,dir(barcord_clusterF)));
}

RE_clusterF=paste0(outfolder,"RE_cluster/")
#cmd=paste0("mkdir -p ",RE_clusterF,";rm -f ",RE_clusterF,"*")
#system(command=cmd)


if(!dir.exists(RE_clusterF)){
    dir.create(RE_clusterF,recursive = TRUE )
} else if (length(dir(RE_clusterF))) {
   file.remove(paste0(RE_clusterF,dir(RE_clusterF)));
}

peak_clusterF=paste0(outfolder,"peak_cluster/")
#cmd=paste0("mkdir -p ",peak_clusterF,";rm -f ",peak_clusterF,"*")
#system(command=cmd)
if(!dir.exists(peak_clusterF)){
    dir.create(peak_clusterF,recursive = TRUE )
} else if (length(dir(peak_clusterF))) {
   file.remove(paste0(peak_clusterF,dir(peak_clusterF)));
}

#cmd=paste0("line=`wc -l ", barcord_cluster_whole ,"|cut -f 1 -d ' ' `
#            echo $line
#           awk -v l=$line -v folder=",barcord_clusterF," \'{if(NR<=l){col[$1]=$2} else{if(length(col[$4])!=0){print >>folder col[$4]\".bed\"}}}\' ",barcord_cluster_whole," ",fragment)

cmd=paste0(awkpath," -v l=",clusterL," -v folder=",barcord_clusterF," \'{if(NR<=l){col[$1]=$2} else{if(length(col[$4])!=0){print >>folder col[$4]\".bed\"}}}\' ",barcord_cluster_whole," ",fragment)

system(command=cmd)



#cmd=paste0("for i in `ls ",barcord_clusterF,"`;
#           do i=`echo $i|tr -d \'.bed\'`;",
#           macs2path," callpeak -t ",barcord_clusterF,"$i.bed -g hs -f BED --nomodel --shift -100 --extsize 200 -n ",peak_clusterF,"$i;"
#           ,bedtoolspath," intersect -a ",peak_clusterF,"${i}_peaks.narrowPeak -b ",oldRegFolder,"Reg_cluster${i}.bed -wa -wb | cut -f 1,2,3,14,15 > ",RE_clusterF,"${i}.bed
#           done")




for(i in cluster){
    cmd=paste0(macs2path," callpeak -t ",barcord_clusterF,i,".bed -g hs -f BED --nomodel --shift -100 --extsize 200 -n ",peak_clusterF,i)
    system(command=cmd)
    genome <- Seqinfo(genome = NA_character_)
    peakfile=paste0(peak_clusterF,i,"_peaks.narrowPeak")
    oldREfile=paste0(oldRegFolder,"Reg_cluster",i,".bed")
    REfile=paste0(RE_clusterF,i,".bed")
    gr_a <- import(peakfile, genome = genome)
    gr_b <- import(oldREfile, genome = genome)
    pairs <- findOverlapPairs(gr_a, gr_b, ignore.strand = TRUE)
    df=data.frame(chr=pairs@first@seqnames,start=pairs@first@ranges@start,end=(pairs@first@ranges@start+pairs@first@ranges@width),TG=pairs@second@elementMetadata@listData$name,score=pairs@second@elementMetadata@listData$score)
    write.table(df,REfile,sep="\t",col.names = F,row.names = F,quote = FALSE)
}


return(c(RE_clusterF=RE_clusterF,peak_clusterF=peak_clusterF))

}
