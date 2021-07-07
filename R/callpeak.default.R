callpeak.default<-function(outfolder,fragment,barcord_cluster_whole,oldRegFolder,macs2path,bedtoolspath){

#outfolder="Rout/"
barcord_clusterF=paste0(outfolder,"barcord_cluster")
cmd=paste0("mkdir -p ",barcord_clusterF,";rm -f ",barcord_clusterF,"/*")
system(command=cmd)

RE_clusterF=paste0(outfolder,"RE_cluster")
cmd=paste0("mkdir -p ",RE_clusterF,";rm -f ",RE_clusterF,"/*")
system(command=cmd)

peak_clusterF=paste0(outfolder,"peak_cluster")
cmd=paste0("mkdir -p ",peak_clusterF,";rm -f ",peak_clusterF,"/*")
system(command=cmd)



#fragment="frag_test.tsv"
#barcord_cluster_whole="outdata/testout/renamed_barcord.txt"
cmd=paste0("line=`wc -l ", barcord_cluster_whole ,"|cut -f 1 -d ' ' `
            echo $line
           awk -v l=$line -v folder=",barcord_clusterF," \'{if(NR<=l){col[$1]=$2} else{if(length(col[$4])!=0){print >>folder\"/\"col[$4]\".bed\"}}}\' ",barcord_cluster_whole," ",fragment)
system(command=cmd)


#macs2path="~/.local/bin/macs2"
#bedtoolspath="/opt/ohpc/pub/Software/bedtools2-2.29.2/bin/bedtools"

cmd=paste0("for i in `ls ",barcord_clusterF,"`;
           do i=`echo $i|tr -d \'.bed\'`;",
           macs2path," callpeak -t ",barcord_clusterF,"/$i.bed -g hs -f BED --nomodel --shift -100 --extsize 200 -n ",peak_clusterF,"/$i;"
           ,bedtoolspath," intersect -a ",peak_clusterF,"/${i}_peaks.narrowPeak -b ",oldRegFolder,"Reg_cluster${i}.bed -wa -wb | cut -f 1,2,3,14,15 > ",RE_clusterF,"/${i}.bed
           done")

system(command=cmd)

return(c(RE_clusterF=RE_clusterF,peak_clusterF=peak_clusterF))

}
