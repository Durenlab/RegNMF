\name{callpeak}
\alias{callpeak}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
call peak
}
\description{
Call for each cluster
}
\usage{
callpeak(outfolder,fragment,barcord_cluster_whole,oldRegFolder,macs2path,awkpath,cluster,clusterL)
}

\arguments{
  \item{outfolder}{result of RegNMF(W123$H)}
  \item{fragment}{fragment path}
  \item{barcord_cluster_whole}{A file contain barcodes and cluster infomation}
  \item{oldRegFolder}{A folder contain predicted regulation for each cluster}
  \item{macs2path}{macs2 path}
  \item{cluster}{Vector of cluster names}
  \item{clusterL}{Number of cluster}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{peak_clusterF}{Folder contain peak infomation for each cluster}
  \item{RE_clusterF}{Folder contain regulation(peak - gene), filtered by peak in each cluster}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Cham
}
