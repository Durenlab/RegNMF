\name{RgeNMF}
\alias{RgeNMF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Reg Nonnegative Matrix Factorization
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
  RgeNMF(E,
         O,
         Symbol,
         PeakName,
         Symbol_location,
         Peak_location,
         feature_cut_perc=0.01,
         corr_cut_k=100000,
         cell_num_smooth=sqrt(dim(E)[2]),
         core=8)
}



%- maybe also 'usage' for other objects documented here.
\arguments{
\item{wholef=in_foldername}{Filtered feature barcode matrix file}
\item{fragment}{ATAC Per fragment information file}
\item{barcode_cluster_whole}{A tsv file contents barcodes with their cluster information}
\item{oldRegFolder}{A folder contents regulations in each cluster as tsv file, will be filtered by peak of each cluster}
\item{macs2path}{Path of macs2}
\item{bedtoolspath}{Path of bedtools}
\item{chr}{Which chromatin you want to see in the result(ex. "chr16").}
\item{from,to}{Which location you want to seein the result.}
\item{width,height}{Figure size of result.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\author{
Cham, fchang@clemson.edu
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
